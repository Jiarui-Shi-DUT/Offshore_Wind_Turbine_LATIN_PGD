# -*- coding: utf-8 -*-
"""
Two-dimensional Euler-Bernoulli beam-column utilities for a tapered
offshore-wind-turbine tower.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from fem.fiber_section import (
    create_annular_fiber_section,
    evaluate_linear_elastic_section,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _validated_integer(value: int, name: str) -> int:
    """Return a validated integer parameter."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


def _validated_finite_scalar(value: float, name: str) -> float:
    """Return a finite floating-point scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error

    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _validated_positive_scalar(value: float, name: str) -> float:
    """Return a finite positive floating-point scalar."""
    result = _validated_finite_scalar(value, name)
    if result <= 0.0:
        raise ValueError(name + " must be positive.")
    return result


@dataclass(frozen=True)
class LinearTaperedTowerGeometry:
    """
    Linearly tapered circular annular tower geometry.

    Outer diameter and wall thickness vary linearly with the tower-axis
    coordinate z measured upward from the tower base.
    """

    height: float
    base_outer_diameter: float
    top_outer_diameter: float
    base_thickness: float
    top_thickness: float

    def __post_init__(self) -> None:
        height = _validated_positive_scalar(self.height, "height")
        base_outer_diameter = _validated_positive_scalar(
            self.base_outer_diameter,
            "base_outer_diameter",
        )
        top_outer_diameter = _validated_positive_scalar(
            self.top_outer_diameter,
            "top_outer_diameter",
        )
        base_thickness = _validated_positive_scalar(
            self.base_thickness,
            "base_thickness",
        )
        top_thickness = _validated_positive_scalar(
            self.top_thickness,
            "top_thickness",
        )

        if 2.0 * base_thickness >= base_outer_diameter:
            raise ValueError(
                "base_thickness must be smaller than the base radius."
            )
        if 2.0 * top_thickness >= top_outer_diameter:
            raise ValueError(
                "top_thickness must be smaller than the top radius."
            )

        object.__setattr__(self, "height", height)
        object.__setattr__(
            self,
            "base_outer_diameter",
            base_outer_diameter,
        )
        object.__setattr__(
            self,
            "top_outer_diameter",
            top_outer_diameter,
        )
        object.__setattr__(
            self,
            "base_thickness",
            base_thickness,
        )
        object.__setattr__(
            self,
            "top_thickness",
            top_thickness,
        )

    def _normalized_height(self, height_coordinate: float) -> float:
        """Validate z and return z / H."""
        z = _validated_finite_scalar(
            height_coordinate,
            "height_coordinate",
        )
        tolerance = 1.0e-12 * max(1.0, self.height)

        if z < -tolerance or z > self.height + tolerance:
            raise ValueError(
                "height_coordinate must lie within the tower height."
            )

        z = float(np.clip(z, 0.0, self.height))
        return z / self.height

    def outer_diameter_at(self, height_coordinate: float) -> float:
        """Return the linearly interpolated outer diameter."""
        ratio = self._normalized_height(height_coordinate)
        return float(
            self.base_outer_diameter
            + ratio
            * (
                self.top_outer_diameter
                - self.base_outer_diameter
            )
        )

    def thickness_at(self, height_coordinate: float) -> float:
        """Return the linearly interpolated wall thickness."""
        ratio = self._normalized_height(height_coordinate)
        return float(
            self.base_thickness
            + ratio
            * (
                self.top_thickness
                - self.base_thickness
            )
        )

    def inner_diameter_at(self, height_coordinate: float) -> float:
        """Return the local inner diameter."""
        outer_diameter = self.outer_diameter_at(height_coordinate)
        thickness = self.thickness_at(height_coordinate)
        return float(outer_diameter - 2.0 * thickness)

    def exact_area_at(self, height_coordinate: float) -> float:
        """Return the analytical local annular-section area."""
        outer_diameter = self.outer_diameter_at(height_coordinate)
        inner_diameter = self.inner_diameter_at(height_coordinate)
        return float(
            0.25
            * np.pi
            * (
                outer_diameter ** 2
                - inner_diameter ** 2
            )
        )

    def exact_second_moment_at(
        self,
        height_coordinate: float,
    ) -> float:
        """Return the analytical local second moment of area."""
        outer_diameter = self.outer_diameter_at(height_coordinate)
        inner_diameter = self.inner_diameter_at(height_coordinate)
        return float(
            np.pi
            / 64.0
            * (
                outer_diameter ** 4
                - inner_diameter ** 4
            )
        )


@dataclass(frozen=True)
class BeamMesh2D:
    """Two-node planar beam mesh with consecutive connectivity."""

    coordinates: FloatArray
    connectivity: IntArray
    tower_axis_coordinates: FloatArray

    def __post_init__(self) -> None:
        coordinates = np.asarray(
            self.coordinates,
            dtype=np.float64,
        )
        connectivity = np.asarray(
            self.connectivity,
            dtype=np.int64,
        )
        tower_axis_coordinates = np.asarray(
            self.tower_axis_coordinates,
            dtype=np.float64,
        )

        if coordinates.ndim != 2 or coordinates.shape[1] != 2:
            raise ValueError(
                "coordinates must have shape (n_nodes, 2)."
            )
        if coordinates.shape[0] < 2:
            raise ValueError(
                "The mesh must contain at least two nodes."
            )
        if np.any(~np.isfinite(coordinates)):
            raise ValueError("coordinates must be finite.")

        if tower_axis_coordinates.ndim != 1:
            raise ValueError(
                "tower_axis_coordinates must be one-dimensional."
            )
        if tower_axis_coordinates.size != coordinates.shape[0]:
            raise ValueError(
                "tower_axis_coordinates must contain one value per node."
            )
        if np.any(~np.isfinite(tower_axis_coordinates)):
            raise ValueError(
                "tower_axis_coordinates must be finite."
            )
        if np.any(np.diff(tower_axis_coordinates) <= 0.0):
            raise ValueError(
                "tower_axis_coordinates must be strictly increasing."
            )

        if connectivity.ndim != 2 or connectivity.shape[1] != 2:
            raise ValueError(
                "connectivity must have shape (n_elements, 2)."
            )
        if connectivity.shape[0] != coordinates.shape[0] - 1:
            raise ValueError(
                "A continuous beam mesh must have n_nodes - 1 elements."
            )

        expected = np.column_stack(
            (
                np.arange(coordinates.shape[0] - 1, dtype=np.int64),
                np.arange(1, coordinates.shape[0], dtype=np.int64),
            )
        )
        if not np.array_equal(connectivity, expected):
            raise ValueError(
                "Elements must connect consecutive nodes in ascending order."
            )

        element_vectors = (
            coordinates[connectivity[:, 1]]
            - coordinates[connectivity[:, 0]]
        )
        element_lengths = np.linalg.norm(
            element_vectors,
            axis=1,
        )
        if np.any(element_lengths <= 0.0):
            raise ValueError("All beam elements must have positive length.")

        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "connectivity", connectivity)
        object.__setattr__(
            self,
            "tower_axis_coordinates",
            tower_axis_coordinates,
        )

    @property
    def n_nodes(self) -> int:
        """Return the number of mesh nodes."""
        return int(self.coordinates.shape[0])

    @property
    def n_elements(self) -> int:
        """Return the number of beam elements."""
        return int(self.connectivity.shape[0])

    @property
    def n_dof(self) -> int:
        """Return the total number of planar frame degrees of freedom."""
        return 3 * self.n_nodes

    @property
    def element_lengths(self) -> FloatArray:
        """Return Euclidean element lengths."""
        start = self.coordinates[self.connectivity[:, 0]]
        end = self.coordinates[self.connectivity[:, 1]]
        return np.linalg.norm(end - start, axis=1)


@dataclass(frozen=True)
class ElasticBeamElementResult:
    """Stiffness and integration data for one elastic beam element."""

    length: float
    transformation: FloatArray
    local_stiffness: FloatArray
    global_stiffness: FloatArray
    gauss_coordinates: FloatArray
    gauss_weights: FloatArray
    gauss_heights: FloatArray
    outer_diameters: FloatArray
    thicknesses: FloatArray
    section_areas: FloatArray
    section_second_moments: FloatArray

    def __post_init__(self) -> None:
        length = _validated_positive_scalar(self.length, "length")
        transformation = np.asarray(
            self.transformation,
            dtype=np.float64,
        )
        local_stiffness = np.asarray(
            self.local_stiffness,
            dtype=np.float64,
        )
        global_stiffness = np.asarray(
            self.global_stiffness,
            dtype=np.float64,
        )

        if transformation.shape != (6, 6):
            raise ValueError(
                "transformation must have shape (6, 6)."
            )
        if local_stiffness.shape != (6, 6):
            raise ValueError(
                "local_stiffness must have shape (6, 6)."
            )
        if global_stiffness.shape != (6, 6):
            raise ValueError(
                "global_stiffness must have shape (6, 6)."
            )

        vectors = (
            np.asarray(self.gauss_coordinates, dtype=np.float64),
            np.asarray(self.gauss_weights, dtype=np.float64),
            np.asarray(self.gauss_heights, dtype=np.float64),
            np.asarray(self.outer_diameters, dtype=np.float64),
            np.asarray(self.thicknesses, dtype=np.float64),
            np.asarray(self.section_areas, dtype=np.float64),
            np.asarray(
                self.section_second_moments,
                dtype=np.float64,
            ),
        )

        size = vectors[0].size
        for vector in vectors:
            if vector.ndim != 1 or vector.size != size:
                raise ValueError(
                    "All Gauss-point arrays must be one-dimensional "
                    "and have equal length."
                )
            if np.any(~np.isfinite(vector)):
                raise ValueError(
                    "All Gauss-point arrays must be finite."
                )

        if np.any(vectors[1] <= 0.0):
            raise ValueError("gauss_weights must be positive.")
        if np.any(vectors[3] <= 0.0):
            raise ValueError("outer_diameters must be positive.")
        if np.any(vectors[4] <= 0.0):
            raise ValueError("thicknesses must be positive.")
        if np.any(vectors[5] <= 0.0):
            raise ValueError("section_areas must be positive.")
        if np.any(vectors[6] <= 0.0):
            raise ValueError(
                "section_second_moments must be positive."
            )

        matrices = (
            transformation,
            local_stiffness,
            global_stiffness,
        )
        for matrix in matrices:
            if np.any(~np.isfinite(matrix)):
                raise ValueError("Element matrices must be finite.")

        object.__setattr__(self, "length", length)
        object.__setattr__(
            self,
            "transformation",
            transformation,
        )
        object.__setattr__(
            self,
            "local_stiffness",
            local_stiffness,
        )
        object.__setattr__(
            self,
            "global_stiffness",
            global_stiffness,
        )
        object.__setattr__(
            self,
            "gauss_coordinates",
            vectors[0],
        )
        object.__setattr__(self, "gauss_weights", vectors[1])
        object.__setattr__(self, "gauss_heights", vectors[2])
        object.__setattr__(self, "outer_diameters", vectors[3])
        object.__setattr__(self, "thicknesses", vectors[4])
        object.__setattr__(self, "section_areas", vectors[5])
        object.__setattr__(
            self,
            "section_second_moments",
            vectors[6],
        )


def create_uniform_vertical_tower_mesh(
    height: float,
    n_elements: int,
    horizontal_coordinate: float = 0.0,
) -> BeamMesh2D:
    """Create a straight vertical tower mesh with uniform element length."""
    height = _validated_positive_scalar(height, "height")
    n_elements = _validated_integer(n_elements, "n_elements")
    horizontal_coordinate = _validated_finite_scalar(
        horizontal_coordinate,
        "horizontal_coordinate",
    )

    if n_elements < 1:
        raise ValueError("n_elements must be at least 1.")

    tower_axis_coordinates = np.linspace(
        0.0,
        height,
        n_elements + 1,
        dtype=np.float64,
    )
    coordinates = np.column_stack(
        (
            np.full(
                n_elements + 1,
                horizontal_coordinate,
                dtype=np.float64,
            ),
            tower_axis_coordinates,
        )
    )
    connectivity = np.column_stack(
        (
            np.arange(n_elements, dtype=np.int64),
            np.arange(1, n_elements + 1, dtype=np.int64),
        )
    )

    return BeamMesh2D(
        coordinates=coordinates,
        connectivity=connectivity,
        tower_axis_coordinates=tower_axis_coordinates,
    )


def gauss_legendre_rule(
    n_points: int = 4,
) -> Tuple[FloatArray, FloatArray]:
    """Return Gauss-Legendre coordinates and weights on [-1, 1]."""
    n_points = _validated_integer(n_points, "n_points")
    if n_points < 1:
        raise ValueError("n_points must be at least 1.")

    coordinates, weights = np.polynomial.legendre.leggauss(n_points)
    return (
        np.asarray(coordinates, dtype=np.float64),
        np.asarray(weights, dtype=np.float64),
    )


def planar_frame_transformation(
    node_coordinates: FloatArray,
) -> FloatArray:
    """
    Return the 6-by-6 transformation from global to local DOFs.

    The local element x-axis points from node 1 to node 2.
    """
    coordinates = np.asarray(
        node_coordinates,
        dtype=np.float64,
    )
    if coordinates.shape != (2, 2):
        raise ValueError(
            "node_coordinates must have shape (2, 2)."
        )
    if np.any(~np.isfinite(coordinates)):
        raise ValueError("node_coordinates must be finite.")

    delta = coordinates[1] - coordinates[0]
    length = float(np.linalg.norm(delta))
    if length <= 0.0:
        raise ValueError(
            "Beam-element node coordinates must be distinct."
        )

    cosine = float(delta[0] / length)
    sine = float(delta[1] / length)

    rotation = np.array(
        [
            [cosine, sine, 0.0],
            [-sine, cosine, 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=np.float64,
    )

    transformation = np.zeros((6, 6), dtype=np.float64)
    transformation[:3, :3] = rotation
    transformation[3:, 3:] = rotation
    return transformation


def euler_bernoulli_strain_displacement(
    length: float,
    natural_coordinate: float,
) -> FloatArray:
    """
    Return the generalized-strain matrix B at one natural coordinate.

    Local DOF ordering is

        [u1, v1, theta1, u2, v2, theta2].

    Generalized strains are

        [axial_strain, curvature],

    with curvature equal to d2v/dx2.
    """
    length = _validated_positive_scalar(length, "length")
    xi = _validated_finite_scalar(
        natural_coordinate,
        "natural_coordinate",
    )
    if xi < -1.0 or xi > 1.0:
        raise ValueError(
            "natural_coordinate must lie in [-1, 1]."
        )

    s = 0.5 * (xi + 1.0)

    axial = np.array(
        [
            -1.0 / length,
            0.0,
            0.0,
            1.0 / length,
            0.0,
            0.0,
        ],
        dtype=np.float64,
    )

    curvature = np.array(
        [
            0.0,
            (-6.0 + 12.0 * s) / length ** 2,
            (-4.0 + 6.0 * s) / length,
            0.0,
            (6.0 - 12.0 * s) / length ** 2,
            (-2.0 + 6.0 * s) / length,
        ],
        dtype=np.float64,
    )

    return np.vstack((axial, curvature))


def compute_elastic_beam_element_stiffness(
    node_coordinates: FloatArray,
    tower_axis_start: float,
    tower_axis_end: float,
    tower_geometry: LinearTaperedTowerGeometry,
    elastic_modulus: float,
    n_gauss: int = 4,
    n_circumferential: int = 32,
    n_radial: int = 2,
) -> ElasticBeamElementResult:
    """
    Integrate one variable-section Euler-Bernoulli beam stiffness matrix.

    Local geometry and fiber-section properties are evaluated separately
    at every Gauss-Legendre integration point.
    """
    if not isinstance(tower_geometry, LinearTaperedTowerGeometry):
        raise TypeError(
            "tower_geometry must be a LinearTaperedTowerGeometry."
        )

    coordinates = np.asarray(
        node_coordinates,
        dtype=np.float64,
    )
    if coordinates.shape != (2, 2):
        raise ValueError(
            "node_coordinates must have shape (2, 2)."
        )
    if np.any(~np.isfinite(coordinates)):
        raise ValueError("node_coordinates must be finite.")

    tower_axis_start = _validated_finite_scalar(
        tower_axis_start,
        "tower_axis_start",
    )
    tower_axis_end = _validated_finite_scalar(
        tower_axis_end,
        "tower_axis_end",
    )
    if tower_axis_end <= tower_axis_start:
        raise ValueError(
            "tower_axis_end must be greater than tower_axis_start."
        )

    tolerance = 1.0e-12 * max(1.0, tower_geometry.height)
    if tower_axis_start < -tolerance:
        raise ValueError(
            "tower_axis_start must lie within the tower."
        )
    if tower_axis_end > tower_geometry.height + tolerance:
        raise ValueError(
            "tower_axis_end must lie within the tower."
        )

    elastic_modulus = _validated_positive_scalar(
        elastic_modulus,
        "elastic_modulus",
    )
    n_gauss = _validated_integer(n_gauss, "n_gauss")
    n_circumferential = _validated_integer(
        n_circumferential,
        "n_circumferential",
    )
    n_radial = _validated_integer(n_radial, "n_radial")

    transformation = planar_frame_transformation(coordinates)
    element_length = float(
        np.linalg.norm(coordinates[1] - coordinates[0])
    )
    gauss_coordinates, gauss_weights = gauss_legendre_rule(
        n_gauss
    )

    gauss_heights = (
        0.5
        * (
            tower_axis_start
            + tower_axis_end
        )
        + 0.5
        * (
            tower_axis_end
            - tower_axis_start
        )
        * gauss_coordinates
    )

    outer_diameters = np.empty(n_gauss, dtype=np.float64)
    thicknesses = np.empty(n_gauss, dtype=np.float64)
    section_areas = np.empty(n_gauss, dtype=np.float64)
    section_second_moments = np.empty(
        n_gauss,
        dtype=np.float64,
    )

    local_stiffness = np.zeros((6, 6), dtype=np.float64)
    jacobian = 0.5 * element_length

    for index in range(n_gauss):
        height_coordinate = float(gauss_heights[index])
        outer_diameter = tower_geometry.outer_diameter_at(
            height_coordinate
        )
        thickness = tower_geometry.thickness_at(
            height_coordinate
        )

        section = create_annular_fiber_section(
            outer_diameter=outer_diameter,
            thickness=thickness,
            n_circumferential=n_circumferential,
            n_radial=n_radial,
        )
        section_response = evaluate_linear_elastic_section(
            section=section,
            axial_strain=0.0,
            curvature=0.0,
            elastic_modulus=elastic_modulus,
        )
        strain_displacement = (
            euler_bernoulli_strain_displacement(
                length=element_length,
                natural_coordinate=float(
                    gauss_coordinates[index]
                ),
            )
        )

        local_stiffness += (
            strain_displacement.T
            @ section_response.tangent
            @ strain_displacement
            * gauss_weights[index]
            * jacobian
        )

        outer_diameters[index] = outer_diameter
        thicknesses[index] = thickness
        section_areas[index] = section.area
        section_second_moments[index] = (
            section.second_moment_x
        )

    local_stiffness = 0.5 * (
        local_stiffness + local_stiffness.T
    )
    global_stiffness = (
        transformation.T
        @ local_stiffness
        @ transformation
    )
    global_stiffness = 0.5 * (
        global_stiffness + global_stiffness.T
    )

    return ElasticBeamElementResult(
        length=element_length,
        transformation=transformation,
        local_stiffness=local_stiffness,
        global_stiffness=global_stiffness,
        gauss_coordinates=gauss_coordinates,
        gauss_weights=gauss_weights,
        gauss_heights=gauss_heights,
        outer_diameters=outer_diameters,
        thicknesses=thicknesses,
        section_areas=section_areas,
        section_second_moments=section_second_moments,
    )

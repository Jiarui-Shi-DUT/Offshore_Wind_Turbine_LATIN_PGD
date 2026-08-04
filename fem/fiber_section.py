# -*- coding: utf-8 -*-
"""
Fiber discretization and section-response utilities for thin-walled
circular annular sections.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class AnnularFiberSection:
    """
    Polar fiber discretization of a circular annular cross-section.

    Each fiber represents one annular sector. Its area and centroid are
    evaluated analytically rather than approximated by a rectangular cell.
    """

    outer_radius: float
    inner_radius: float
    n_circumferential: int
    n_radial: int
    areas: FloatArray
    x_coordinates: FloatArray
    y_coordinates: FloatArray
    radial_indices: IntArray
    circumferential_indices: IntArray

    def __post_init__(self) -> None:
        outer_radius = float(self.outer_radius)
        inner_radius = float(self.inner_radius)
        n_circumferential = int(self.n_circumferential)
        n_radial = int(self.n_radial)

        areas = np.asarray(self.areas, dtype=np.float64)
        x_coordinates = np.asarray(
            self.x_coordinates,
            dtype=np.float64,
        )
        y_coordinates = np.asarray(
            self.y_coordinates,
            dtype=np.float64,
        )
        radial_indices = np.asarray(
            self.radial_indices,
            dtype=np.int64,
        )
        circumferential_indices = np.asarray(
            self.circumferential_indices,
            dtype=np.int64,
        )

        if not np.isfinite(outer_radius) or outer_radius <= 0.0:
            raise ValueError("outer_radius must be finite and positive.")
        if not np.isfinite(inner_radius) or inner_radius < 0.0:
            raise ValueError(
                "inner_radius must be finite and non-negative."
            )
        if inner_radius >= outer_radius:
            raise ValueError(
                "inner_radius must be smaller than outer_radius."
            )

        if n_circumferential < 4:
            raise ValueError(
                "n_circumferential must be at least 4."
            )
        if n_circumferential % 4 != 0:
            raise ValueError(
                "n_circumferential must be divisible by 4."
            )
        if n_radial < 1:
            raise ValueError("n_radial must be at least 1.")

        expected_size = n_circumferential * n_radial
        arrays = (
            areas,
            x_coordinates,
            y_coordinates,
            radial_indices,
            circumferential_indices,
        )

        for array in arrays:
            if array.ndim != 1:
                raise ValueError("Fiber arrays must be one-dimensional.")
            if array.size != expected_size:
                raise ValueError(
                    "Fiber arrays must contain "
                    "n_circumferential * n_radial entries."
                )

        if np.any(~np.isfinite(areas)) or np.any(areas <= 0.0):
            raise ValueError(
                "All fiber areas must be finite and positive."
            )
        if np.any(~np.isfinite(x_coordinates)):
            raise ValueError(
                "All fiber x-coordinates must be finite."
            )
        if np.any(~np.isfinite(y_coordinates)):
            raise ValueError(
                "All fiber y-coordinates must be finite."
            )

        if np.any(radial_indices < 0):
            raise ValueError("radial_indices must be non-negative.")
        if np.any(radial_indices >= n_radial):
            raise ValueError(
                "radial_indices exceed the radial layer count."
            )
        if np.any(circumferential_indices < 0):
            raise ValueError(
                "circumferential_indices must be non-negative."
            )
        if np.any(
            circumferential_indices >= n_circumferential
        ):
            raise ValueError(
                "circumferential_indices exceed the sector count."
            )

        object.__setattr__(self, "outer_radius", outer_radius)
        object.__setattr__(self, "inner_radius", inner_radius)
        object.__setattr__(
            self,
            "n_circumferential",
            n_circumferential,
        )
        object.__setattr__(self, "n_radial", n_radial)
        object.__setattr__(self, "areas", areas)
        object.__setattr__(
            self,
            "x_coordinates",
            x_coordinates,
        )
        object.__setattr__(
            self,
            "y_coordinates",
            y_coordinates,
        )
        object.__setattr__(
            self,
            "radial_indices",
            radial_indices,
        )
        object.__setattr__(
            self,
            "circumferential_indices",
            circumferential_indices,
        )

    @property
    def n_fibers(self) -> int:
        """Return the total number of section fibers."""
        return int(self.areas.size)

    @property
    def area(self) -> float:
        """Return the area represented by all fibers."""
        return float(np.sum(self.areas))

    @property
    def exact_area(self) -> float:
        """Return the analytical annular-section area."""
        return float(
            np.pi
            * (
                self.outer_radius ** 2
                - self.inner_radius ** 2
            )
        )

    @property
    def first_moment_y(self) -> float:
        """Return the first moment of area with respect to y."""
        return float(np.dot(self.areas, self.y_coordinates))

    @property
    def centroid_x(self) -> float:
        """Return the fiber-discretized x-coordinate of the centroid."""
        return float(
            np.dot(self.areas, self.x_coordinates) / self.area
        )

    @property
    def centroid_y(self) -> float:
        """Return the fiber-discretized y-coordinate of the centroid."""
        return float(self.first_moment_y / self.area)

    @property
    def second_moment_x(self) -> float:
        """
        Return the fiber-discretized second moment about the x-axis.

        Fibers are represented as concentrated areas at their exact
        annular-sector centroids.
        """
        return float(
            np.dot(self.areas, self.y_coordinates ** 2)
        )

    @property
    def second_moment_y(self) -> float:
        """
        Return the fiber-discretized second moment about the y-axis.

        Fibers are represented as concentrated areas at their exact
        annular-sector centroids.
        """
        return float(
            np.dot(self.areas, self.x_coordinates ** 2)
        )

    @property
    def product_moment_xy(self) -> float:
        """Return the fiber-discretized product moment of area."""
        return float(
            np.dot(
                self.areas,
                self.x_coordinates * self.y_coordinates,
            )
        )

    @property
    def exact_second_moment(self) -> float:
        """
        Return the analytical second moment about either centroidal axis.
        """
        return float(
            0.25
            * np.pi
            * (
                self.outer_radius ** 4
                - self.inner_radius ** 4
            )
        )

    @property
    def coordinates(self) -> FloatArray:
        """Return fiber centroid coordinates with shape (n_fibers, 2)."""
        return np.column_stack(
            (
                self.x_coordinates,
                self.y_coordinates,
            )
        )


@dataclass(frozen=True)
class ElasticSectionResponse:
    """
    Linear-elastic response of one annular fiber section.

    The section kinematic convention is

        epsilon_f = epsilon_0 - kappa * y_f

    and the bending-moment convention is

        M = -sum(sigma_f * y_f * A_f).

    Consequently, positive curvature produces a positive bending moment.
    """

    axial_strain: float
    curvature: float
    elastic_modulus: float
    fiber_strains: FloatArray
    fiber_stresses: FloatArray
    axial_force: float
    bending_moment: float
    tangent: FloatArray

    def __post_init__(self) -> None:
        axial_strain = _validated_finite_scalar(
            self.axial_strain,
            "axial_strain",
        )
        curvature = _validated_finite_scalar(
            self.curvature,
            "curvature",
        )
        elastic_modulus = _validated_positive_scalar(
            self.elastic_modulus,
            "elastic_modulus",
        )

        fiber_strains = np.asarray(
            self.fiber_strains,
            dtype=np.float64,
        )
        fiber_stresses = np.asarray(
            self.fiber_stresses,
            dtype=np.float64,
        )
        tangent = np.asarray(self.tangent, dtype=np.float64)

        if fiber_strains.ndim != 1:
            raise ValueError(
                "fiber_strains must be one-dimensional."
            )
        if fiber_stresses.shape != fiber_strains.shape:
            raise ValueError(
                "fiber_stresses must match fiber_strains."
            )
        if np.any(~np.isfinite(fiber_strains)):
            raise ValueError("fiber_strains must be finite.")
        if np.any(~np.isfinite(fiber_stresses)):
            raise ValueError("fiber_stresses must be finite.")
        if tangent.shape != (2, 2):
            raise ValueError("tangent must have shape (2, 2).")
        if np.any(~np.isfinite(tangent)):
            raise ValueError("tangent must be finite.")

        axial_force = _validated_finite_scalar(
            self.axial_force,
            "axial_force",
        )
        bending_moment = _validated_finite_scalar(
            self.bending_moment,
            "bending_moment",
        )

        object.__setattr__(self, "axial_strain", axial_strain)
        object.__setattr__(self, "curvature", curvature)
        object.__setattr__(
            self,
            "elastic_modulus",
            elastic_modulus,
        )
        object.__setattr__(
            self,
            "fiber_strains",
            fiber_strains,
        )
        object.__setattr__(
            self,
            "fiber_stresses",
            fiber_stresses,
        )
        object.__setattr__(self, "axial_force", axial_force)
        object.__setattr__(
            self,
            "bending_moment",
            bending_moment,
        )
        object.__setattr__(self, "tangent", tangent)

    @property
    def resultants(self) -> FloatArray:
        """Return [axial_force, bending_moment]."""
        return np.array(
            [self.axial_force, self.bending_moment],
            dtype=np.float64,
        )


def _validated_integer(value: int, name: str) -> int:
    """Return a validated integer discretization parameter."""
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


def create_annular_fiber_section(
    outer_diameter: float,
    thickness: float,
    n_circumferential: int = 32,
    n_radial: int = 2,
) -> AnnularFiberSection:
    """
    Create an annular-sector fiber discretization.

    Circumferential sector centroids are placed at angles

        0, delta_theta, 2 * delta_theta, ...

    so that fibers are aligned with the positive and negative coordinate
    axes whenever n_circumferential is divisible by four. Radial layers
    have equal thickness.

    Parameters
    ----------
    outer_diameter
        Outside diameter of the circular section.
    thickness
        Wall thickness of the circular section.
    n_circumferential
        Number of equal circumferential sectors.
    n_radial
        Number of equal radial layers through the wall thickness.
    """
    outer_diameter = _validated_positive_scalar(
        outer_diameter,
        "outer_diameter",
    )
    thickness = _validated_positive_scalar(
        thickness,
        "thickness",
    )
    n_circumferential = _validated_integer(
        n_circumferential,
        "n_circumferential",
    )
    n_radial = _validated_integer(
        n_radial,
        "n_radial",
    )

    outer_radius = 0.5 * outer_diameter
    inner_radius = outer_radius - thickness

    if inner_radius <= 0.0:
        raise ValueError(
            "thickness must be smaller than the outer radius."
        )
    if n_circumferential < 4:
        raise ValueError(
            "n_circumferential must be at least 4."
        )
    if n_circumferential % 4 != 0:
        raise ValueError(
            "n_circumferential must be divisible by 4."
        )
    if n_radial < 1:
        raise ValueError("n_radial must be at least 1.")

    delta_theta = 2.0 * np.pi / n_circumferential

    radial_edges = np.linspace(
        inner_radius,
        outer_radius,
        n_radial + 1,
        dtype=np.float64,
    )

    radial_inner = np.repeat(
        radial_edges[:-1],
        n_circumferential,
    )
    radial_outer = np.repeat(
        radial_edges[1:],
        n_circumferential,
    )

    circumferential_indices = np.tile(
        np.arange(n_circumferential, dtype=np.int64),
        n_radial,
    )
    radial_indices = np.repeat(
        np.arange(n_radial, dtype=np.int64),
        n_circumferential,
    )

    centroid_angles = (
        circumferential_indices.astype(np.float64)
        * delta_theta
    )

    areas = (
        0.5
        * (
            radial_outer ** 2
            - radial_inner ** 2
        )
        * delta_theta
    )

    centroid_radii = (
        4.0
        * np.sin(0.5 * delta_theta)
        / (3.0 * delta_theta)
        * (
            radial_outer ** 3
            - radial_inner ** 3
        )
        / (
            radial_outer ** 2
            - radial_inner ** 2
        )
    )

    x_coordinates = centroid_radii * np.cos(
        centroid_angles
    )
    y_coordinates = centroid_radii * np.sin(
        centroid_angles
    )

    return AnnularFiberSection(
        outer_radius=outer_radius,
        inner_radius=inner_radius,
        n_circumferential=n_circumferential,
        n_radial=n_radial,
        areas=areas,
        x_coordinates=x_coordinates,
        y_coordinates=y_coordinates,
        radial_indices=radial_indices,
        circumferential_indices=circumferential_indices,
    )


def compute_fiber_strains(
    section: AnnularFiberSection,
    axial_strain: float,
    curvature: float,
) -> FloatArray:
    """
    Evaluate the Euler-Bernoulli axial strain at every fiber centroid.

    Positive curvature gives compression at fibers with positive y.
    """
    if not isinstance(section, AnnularFiberSection):
        raise TypeError(
            "section must be an AnnularFiberSection."
        )

    axial_strain = _validated_finite_scalar(
        axial_strain,
        "axial_strain",
    )
    curvature = _validated_finite_scalar(
        curvature,
        "curvature",
    )

    return (
        axial_strain
        - curvature * section.y_coordinates
    ).astype(np.float64, copy=False)


def integrate_fiber_stresses(
    section: AnnularFiberSection,
    fiber_stresses: FloatArray,
) -> FloatArray:
    """
    Integrate fiber stresses into [axial_force, bending_moment].

    The bending-moment sign convention is

        M = -sum(sigma_f * y_f * A_f).
    """
    if not isinstance(section, AnnularFiberSection):
        raise TypeError(
            "section must be an AnnularFiberSection."
        )

    stresses = np.asarray(fiber_stresses, dtype=np.float64)

    if stresses.ndim != 1:
        raise ValueError(
            "fiber_stresses must be one-dimensional."
        )
    if stresses.size != section.n_fibers:
        raise ValueError(
            "fiber_stresses must contain one value per fiber."
        )
    if np.any(~np.isfinite(stresses)):
        raise ValueError("fiber_stresses must be finite.")

    axial_force = float(np.dot(section.areas, stresses))
    bending_moment = float(
        -np.dot(
            section.areas * section.y_coordinates,
            stresses,
        )
    )

    return np.array(
        [axial_force, bending_moment],
        dtype=np.float64,
    )


def evaluate_linear_elastic_section(
    section: AnnularFiberSection,
    axial_strain: float,
    curvature: float,
    elastic_modulus: float,
) -> ElasticSectionResponse:
    """
    Evaluate fiber stresses, section resultants, and tangent stiffness.

    The generalized strain and resultant vectors are

        e = [epsilon_0, kappa]
        s = [N, M].

    The returned tangent satisfies ds = K_section de.
    """
    if not isinstance(section, AnnularFiberSection):
        raise TypeError(
            "section must be an AnnularFiberSection."
        )

    axial_strain = _validated_finite_scalar(
        axial_strain,
        "axial_strain",
    )
    curvature = _validated_finite_scalar(
        curvature,
        "curvature",
    )
    elastic_modulus = _validated_positive_scalar(
        elastic_modulus,
        "elastic_modulus",
    )

    fiber_strains = compute_fiber_strains(
        section=section,
        axial_strain=axial_strain,
        curvature=curvature,
    )
    fiber_stresses = elastic_modulus * fiber_strains
    resultants = integrate_fiber_stresses(
        section=section,
        fiber_stresses=fiber_stresses,
    )

    coupling = -elastic_modulus * section.first_moment_y
    tangent = np.array(
        [
            [
                elastic_modulus * section.area,
                coupling,
            ],
            [
                coupling,
                elastic_modulus * section.second_moment_x,
            ],
        ],
        dtype=np.float64,
    )

    return ElasticSectionResponse(
        axial_strain=axial_strain,
        curvature=curvature,
        elastic_modulus=elastic_modulus,
        fiber_strains=fiber_strains,
        fiber_stresses=fiber_stresses,
        axial_force=float(resultants[0]),
        bending_moment=float(resultants[1]),
        tangent=tangent,
    )

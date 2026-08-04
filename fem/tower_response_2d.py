# -*- coding: utf-8 -*-
"""
Elastic response recovery for a two-dimensional tapered tower.

This module maps the solved global nodal displacement vector back to:

    global element DOFs
        -> local element DOFs
        -> Gauss-point generalized strains
        -> annular-fiber strains and stresses
        -> section axial force and bending moment.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    euler_bernoulli_strain_displacement,
)
from fem.fiber_section import (
    create_annular_fiber_section,
    evaluate_linear_elastic_section,
)
from fem.tower_system_2d import (
    ElasticTowerAssembly,
    LinearStaticSolution,
)


FloatArray = NDArray[np.float64]


def _validated_positive_scalar(
    value: float,
    name: str,
) -> float:
    """Return a finite positive scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error

    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    if result <= 0.0:
        raise ValueError(name + " must be positive.")
    return result


def _validated_integer(
    value: int,
    name: str,
) -> int:
    """Return a validated integer."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


@dataclass(frozen=True)
class ElasticTowerResponse:
    """
    Recovered elastic response at all beam integration points.

    Array shapes
    ------------
    local_displacements
        (n_elements, 6)
    gauss_heights and all scalar Gauss-point response arrays
        (n_elements, n_gauss)
    fiber_strains and fiber_stresses
        (n_elements, n_gauss, n_fibers)
    """

    local_displacements: FloatArray
    gauss_heights: FloatArray
    axial_strains: FloatArray
    curvatures: FloatArray
    axial_forces: FloatArray
    bending_moments: FloatArray
    fiber_strains: FloatArray
    fiber_stresses: FloatArray

    def __post_init__(self) -> None:
        local_displacements = np.asarray(
            self.local_displacements,
            dtype=np.float64,
        )
        gauss_heights = np.asarray(
            self.gauss_heights,
            dtype=np.float64,
        )
        axial_strains = np.asarray(
            self.axial_strains,
            dtype=np.float64,
        )
        curvatures = np.asarray(
            self.curvatures,
            dtype=np.float64,
        )
        axial_forces = np.asarray(
            self.axial_forces,
            dtype=np.float64,
        )
        bending_moments = np.asarray(
            self.bending_moments,
            dtype=np.float64,
        )
        fiber_strains = np.asarray(
            self.fiber_strains,
            dtype=np.float64,
        )
        fiber_stresses = np.asarray(
            self.fiber_stresses,
            dtype=np.float64,
        )

        if (
            local_displacements.ndim != 2
            or local_displacements.shape[1] != 6
        ):
            raise ValueError(
                "local_displacements must have shape "
                "(n_elements, 6)."
            )

        scalar_arrays = (
            gauss_heights,
            axial_strains,
            curvatures,
            axial_forces,
            bending_moments,
        )
        scalar_shape = gauss_heights.shape

        if gauss_heights.ndim != 2:
            raise ValueError(
                "gauss_heights must have shape "
                "(n_elements, n_gauss)."
            )

        for array in scalar_arrays:
            if array.shape != scalar_shape:
                raise ValueError(
                    "All scalar Gauss-point arrays must "
                    "have equal shape."
                )
            if np.any(~np.isfinite(array)):
                raise ValueError(
                    "All scalar Gauss-point arrays must be finite."
                )

        if scalar_shape[0] != local_displacements.shape[0]:
            raise ValueError(
                "The response must contain the same number "
                "of elements in every array."
            )

        if fiber_strains.ndim != 3:
            raise ValueError(
                "fiber_strains must have shape "
                "(n_elements, n_gauss, n_fibers)."
            )
        if fiber_stresses.shape != fiber_strains.shape:
            raise ValueError(
                "fiber_stresses must match fiber_strains."
            )
        if fiber_strains.shape[:2] != scalar_shape:
            raise ValueError(
                "Fiber response arrays must match the "
                "element and Gauss-point dimensions."
            )
        if fiber_strains.shape[2] < 1:
            raise ValueError(
                "At least one section fiber is required."
            )
        if np.any(~np.isfinite(fiber_strains)):
            raise ValueError("fiber_strains must be finite.")
        if np.any(~np.isfinite(fiber_stresses)):
            raise ValueError("fiber_stresses must be finite.")

        object.__setattr__(
            self,
            "local_displacements",
            local_displacements,
        )
        object.__setattr__(
            self,
            "gauss_heights",
            gauss_heights,
        )
        object.__setattr__(
            self,
            "axial_strains",
            axial_strains,
        )
        object.__setattr__(
            self,
            "curvatures",
            curvatures,
        )
        object.__setattr__(
            self,
            "axial_forces",
            axial_forces,
        )
        object.__setattr__(
            self,
            "bending_moments",
            bending_moments,
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

    @property
    def n_elements(self) -> int:
        """Return the number of beam elements."""
        return int(self.local_displacements.shape[0])

    @property
    def n_gauss(self) -> int:
        """Return the number of Gauss points per beam element."""
        return int(self.gauss_heights.shape[1])

    @property
    def n_fibers(self) -> int:
        """Return the number of fibers per section."""
        return int(self.fiber_strains.shape[2])

    @property
    def minimum_fiber_strains(self) -> FloatArray:
        """Return the minimum discrete-fiber strain at each Gauss point."""
        return np.min(self.fiber_strains, axis=2)

    @property
    def maximum_fiber_strains(self) -> FloatArray:
        """Return the maximum discrete-fiber strain at each Gauss point."""
        return np.max(self.fiber_strains, axis=2)

    @property
    def minimum_fiber_stresses(self) -> FloatArray:
        """Return the minimum discrete-fiber stress at each Gauss point."""
        return np.min(self.fiber_stresses, axis=2)

    @property
    def maximum_fiber_stresses(self) -> FloatArray:
        """Return the maximum discrete-fiber stress at each Gauss point."""
        return np.max(self.fiber_stresses, axis=2)

    @property
    def flattened_gauss_heights(self) -> FloatArray:
        """Return all Gauss-point heights as one ascending vector."""
        return self.gauss_heights.reshape(-1)

    @property
    def flattened_bending_moments(self) -> FloatArray:
        """Return all Gauss-point bending moments as one vector."""
        return self.bending_moments.reshape(-1)


def recover_elastic_tower_response(
    assembly: ElasticTowerAssembly,
    solution: LinearStaticSolution,
    tower_geometry: LinearTaperedTowerGeometry,
    elastic_modulus: float,
    n_circumferential: int = 32,
    n_radial: int = 2,
) -> ElasticTowerResponse:
    """
    Recover Gauss-point and fiber response from a solved tower model.

    Sign conventions
    ----------------
    Local generalized strains are

        [epsilon_0, kappa],

    and the fiber strain is

        epsilon_f = epsilon_0 - kappa * y_f.

    The section bending moment is

        M = -sum(sigma_f * y_f * A_f),

    so positive curvature produces positive bending moment.
    """
    if not isinstance(assembly, ElasticTowerAssembly):
        raise TypeError(
            "assembly must be an ElasticTowerAssembly."
        )
    if not isinstance(solution, LinearStaticSolution):
        raise TypeError(
            "solution must be a LinearStaticSolution."
        )
    if not isinstance(
        tower_geometry,
        LinearTaperedTowerGeometry,
    ):
        raise TypeError(
            "tower_geometry must be a "
            "LinearTaperedTowerGeometry."
        )

    elastic_modulus = _validated_positive_scalar(
        elastic_modulus,
        "elastic_modulus",
    )
    n_circumferential = _validated_integer(
        n_circumferential,
        "n_circumferential",
    )
    n_radial = _validated_integer(
        n_radial,
        "n_radial",
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

    mesh = assembly.mesh
    if solution.displacements.shape != (mesh.n_dof,):
        raise ValueError(
            "solution displacement size is inconsistent "
            "with the assembled mesh."
        )

    if len(assembly.element_results) == 0:
        raise ValueError(
            "assembly must contain at least one element."
        )

    n_gauss = int(
        assembly.element_results[0].gauss_coordinates.size
    )
    if n_gauss < 1:
        raise ValueError(
            "Each element must contain at least one Gauss point."
        )

    for element_result in assembly.element_results:
        if element_result.gauss_coordinates.size != n_gauss:
            raise ValueError(
                "All elements must use the same number "
                "of Gauss points."
            )

    n_elements = mesh.n_elements
    n_fibers = n_circumferential * n_radial

    local_displacements = np.empty(
        (n_elements, 6),
        dtype=np.float64,
    )
    gauss_heights = np.empty(
        (n_elements, n_gauss),
        dtype=np.float64,
    )
    axial_strains = np.empty_like(gauss_heights)
    curvatures = np.empty_like(gauss_heights)
    axial_forces = np.empty_like(gauss_heights)
    bending_moments = np.empty_like(gauss_heights)
    fiber_strains = np.empty(
        (n_elements, n_gauss, n_fibers),
        dtype=np.float64,
    )
    fiber_stresses = np.empty_like(fiber_strains)

    for element_index in range(n_elements):
        element_result = assembly.element_results[element_index]
        dofs = assembly.element_dofs[element_index]

        element_global_displacements = (
            solution.displacements[dofs]
        )
        element_local_displacements = (
            element_result.transformation
            @ element_global_displacements
        )
        local_displacements[element_index] = (
            element_local_displacements
        )

        for gauss_index in range(n_gauss):
            height_coordinate = float(
                element_result.gauss_heights[gauss_index]
            )
            natural_coordinate = float(
                element_result.gauss_coordinates[gauss_index]
            )

            strain_displacement = (
                euler_bernoulli_strain_displacement(
                    length=element_result.length,
                    natural_coordinate=natural_coordinate,
                )
            )
            generalized_strains = (
                strain_displacement
                @ element_local_displacements
            )
            axial_strain = float(generalized_strains[0])
            curvature = float(generalized_strains[1])

            outer_diameter = (
                tower_geometry.outer_diameter_at(
                    height_coordinate
                )
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
            section_response = (
                evaluate_linear_elastic_section(
                    section=section,
                    axial_strain=axial_strain,
                    curvature=curvature,
                    elastic_modulus=elastic_modulus,
                )
            )

            gauss_heights[
                element_index,
                gauss_index,
            ] = height_coordinate
            axial_strains[
                element_index,
                gauss_index,
            ] = axial_strain
            curvatures[
                element_index,
                gauss_index,
            ] = curvature
            axial_forces[
                element_index,
                gauss_index,
            ] = section_response.axial_force
            bending_moments[
                element_index,
                gauss_index,
            ] = section_response.bending_moment
            fiber_strains[
                element_index,
                gauss_index,
                :,
            ] = section_response.fiber_strains
            fiber_stresses[
                element_index,
                gauss_index,
                :,
            ] = section_response.fiber_stresses

    return ElasticTowerResponse(
        local_displacements=local_displacements,
        gauss_heights=gauss_heights,
        axial_strains=axial_strains,
        curvatures=curvatures,
        axial_forces=axial_forces,
        bending_moments=bending_moments,
        fiber_strains=fiber_strains,
        fiber_stresses=fiber_stresses,
    )

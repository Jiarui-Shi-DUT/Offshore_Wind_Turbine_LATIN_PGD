# -*- coding: utf-8 -*-
"""
Material-point metric and reference structural projection for the tower
LATIN-PGD solver.

The canonical material-point coordinate q is owned by MaterialPointLayout.
For one beam-Gauss-fiber point,

    v_q = J_e * w_g * A_egf,

and the geometric material-point metric is conceptually

    M = diag(v_q).

The dense metric matrix is never formed.  The reference structural projection
uses the fiber compatibility operator H and reference elastic modulus C0:

    E = H (H^T M C0 H)^(-1) H^T M C0.

For a source strain r, the returned compatible strain and stress correction are

    epsilon = E r,
    sigma = C0 (epsilon - r),

so the free-DOF equilibrium residual satisfies H^T M sigma = 0.

This module deliberately does not own LATIN states, PGD bases, constitutive
history variables, or nonlinear solver transactions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from latin.tower_state import MaterialPointLayout


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
IntegratedValue = Union[float, FloatArray]


def _readonly_float_array(
    values: FloatArray,
    name: str,
    ndim: int,
) -> FloatArray:
    """Return a detached finite float64 array with write access disabled."""
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != ndim:
        raise ValueError(
            name + " must have " + str(ndim) + " dimension(s)."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain only finite values.")
    array.setflags(write=False)
    return array


def _readonly_int_array(
    values: Sequence[int],
    name: str,
) -> IntArray:
    """Return a detached one-dimensional int64 array."""
    array = np.array(values, dtype=np.int64, copy=True)
    if array.ndim != 1:
        raise ValueError(name + " must be one-dimensional.")
    array.setflags(write=False)
    return array


def _validated_vector_last_axis(
    values: FloatArray,
    n_material_points: int,
    name: str,
) -> FloatArray:
    """Return a finite float array whose last axis is the material axis."""
    array = np.asarray(values, dtype=np.float64)
    if array.ndim < 1:
        raise ValueError(name + " must have at least one dimension.")
    if array.shape[-1] != n_material_points:
        raise ValueError(
            name
            + " must have last-axis size "
            + str(n_material_points)
            + "."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain only finite values.")
    return array


@dataclass(frozen=True)
class MaterialPointMetric:
    """Immutable geometric integration weights on canonical material points."""

    weights: FloatArray

    def __post_init__(self) -> None:
        weights = _readonly_float_array(
            self.weights,
            "weights",
            ndim=1,
        )
        if weights.size < 1:
            raise ValueError(
                "weights must contain at least one material point."
            )
        if np.any(weights <= 0.0):
            raise ValueError(
                "All material-point metric weights must be positive."
            )
        object.__setattr__(self, "weights", weights)

    @property
    def n_material_points(self) -> int:
        """Return the number of material-point weights."""
        return int(self.weights.size)

    @property
    def total_measure(self) -> float:
        """Return the total material volume measure represented by the metric."""
        return float(np.sum(self.weights))

    def apply(self, values: FloatArray) -> FloatArray:
        """Apply M without forming a dense matrix."""
        array = _validated_vector_last_axis(
            values,
            self.n_material_points,
            "values",
        )
        result = np.asarray(
            array * self.weights,
            dtype=np.float64,
        )
        return result.copy()

    def integrate(self, values: FloatArray) -> IntegratedValue:
        """Integrate over the final canonical material-point axis."""
        array = _validated_vector_last_axis(
            values,
            self.n_material_points,
            "values",
        )
        result = np.sum(array * self.weights, axis=-1)
        if np.ndim(result) == 0:
            return float(result)
        return np.asarray(result, dtype=np.float64)

    def inner_product(
        self,
        left: FloatArray,
        right: FloatArray,
    ) -> float:
        """Return left^T M right for two spatial material-point vectors."""
        left_array = np.asarray(left, dtype=np.float64)
        right_array = np.asarray(right, dtype=np.float64)
        expected_shape = (self.n_material_points,)

        if left_array.shape != expected_shape:
            raise ValueError(
                "left must have shape " + str(expected_shape) + "."
            )
        if right_array.shape != expected_shape:
            raise ValueError(
                "right must have shape " + str(expected_shape) + "."
            )
        if np.any(~np.isfinite(left_array)):
            raise ValueError("left must contain only finite values.")
        if np.any(~np.isfinite(right_array)):
            raise ValueError("right must contain only finite values.")

        return float(
            np.dot(
                self.weights * left_array,
                right_array,
            )
        )

    def norm(self, values: FloatArray) -> float:
        """Return the material-point M norm of one spatial vector."""
        squared = self.inner_product(values, values)
        if squared < 0.0:
            raise FloatingPointError(
                "The material-point metric produced a negative norm square."
            )
        return float(np.sqrt(max(0.0, squared)))


@dataclass(frozen=True)
class EquilibriumProjectionTower:
    """Ephemeral result of one reference tower equilibrium projection."""

    source_strain: FloatArray
    compatible_strain: FloatArray
    stress: FloatArray
    displacement_free: FloatArray

    def __post_init__(self) -> None:
        source = np.array(
            self.source_strain,
            dtype=np.float64,
            copy=True,
        )
        compatible = np.array(
            self.compatible_strain,
            dtype=np.float64,
            copy=True,
        )
        stress = np.array(
            self.stress,
            dtype=np.float64,
            copy=True,
        )
        displacement = np.array(
            self.displacement_free,
            dtype=np.float64,
            copy=True,
        )

        if source.ndim not in (1, 2):
            raise ValueError(
                "source_strain must be spatial or time-by-material."
            )
        if compatible.shape != source.shape:
            raise ValueError(
                "compatible_strain must match source_strain."
            )
        if stress.shape != source.shape:
            raise ValueError("stress must match source_strain.")

        expected_displacement_ndim = source.ndim
        if displacement.ndim != expected_displacement_ndim:
            raise ValueError(
                "displacement_free must have the same rank as source_strain."
            )
        if source.ndim == 2:
            if displacement.shape[0] != source.shape[0]:
                raise ValueError(
                    "History displacement_free must match the time dimension."
                )

        for name, array in (
            ("source_strain", source),
            ("compatible_strain", compatible),
            ("stress", stress),
            ("displacement_free", displacement),
        ):
            if np.any(~np.isfinite(array)):
                raise ValueError(name + " must contain only finite values.")
            array.setflags(write=False)

        object.__setattr__(self, "source_strain", source)
        object.__setattr__(self, "compatible_strain", compatible)
        object.__setattr__(self, "stress", stress)
        object.__setattr__(
            self,
            "displacement_free",
            displacement,
        )

    @property
    def is_history(self) -> bool:
        """Return whether this result stores a time-by-material projection."""
        return bool(self.source_strain.ndim == 2)


@dataclass(frozen=True)
class TowerEquilibriumOperator:
    """
    Immutable reference elastic projection on the free structural DOFs.

    compatibility_matrix has shape (n_material_points, n_free_dofs) and maps
    free structural displacements to fiber axial strains.
    """

    layout: MaterialPointLayout
    metric: MaterialPointMetric
    reference_modulus: FloatArray
    compatibility_matrix: FloatArray
    free_dofs: IntArray
    n_dof: int

    _weighted_modulus: FloatArray = field(
        init=False,
        repr=False,
        compare=False,
    )
    _reduced_stiffness: FloatArray = field(
        init=False,
        repr=False,
        compare=False,
    )
    _cholesky: FloatArray = field(
        init=False,
        repr=False,
        compare=False,
    )

    def __post_init__(self) -> None:
        if not isinstance(self.layout, MaterialPointLayout):
            raise TypeError(
                "layout must be a MaterialPointLayout."
            )
        if not isinstance(self.metric, MaterialPointMetric):
            raise TypeError(
                "metric must be a MaterialPointMetric."
            )
        if (
            self.metric.n_material_points
            != self.layout.n_material_points
        ):
            raise ValueError(
                "layout and metric must have the same material-point count."
            )

        n_dof = int(self.n_dof)
        if n_dof < 1:
            raise ValueError("n_dof must be at least 1.")

        free_dofs = _readonly_int_array(
            self.free_dofs,
            "free_dofs",
        )
        if free_dofs.size < 1:
            raise ValueError(
                "At least one free structural DOF is required."
            )
        if np.any(free_dofs < 0) or np.any(free_dofs >= n_dof):
            raise ValueError(
                "free_dofs contain an out-of-range DOF index."
            )
        if np.unique(free_dofs).size != free_dofs.size:
            raise ValueError(
                "free_dofs must not contain duplicates."
            )

        reference_modulus = _readonly_float_array(
            self.reference_modulus,
            "reference_modulus",
            ndim=1,
        )
        expected_material_shape = (
            self.layout.n_material_points,
        )
        if reference_modulus.shape != expected_material_shape:
            raise ValueError(
                "reference_modulus must have shape "
                + str(expected_material_shape)
                + "."
            )
        if np.any(reference_modulus <= 0.0):
            raise ValueError(
                "reference_modulus must be strictly positive."
            )

        compatibility = _readonly_float_array(
            self.compatibility_matrix,
            "compatibility_matrix",
            ndim=2,
        )
        expected_compatibility_shape = (
            self.layout.n_material_points,
            free_dofs.size,
        )
        if compatibility.shape != expected_compatibility_shape:
            raise ValueError(
                "compatibility_matrix must have shape "
                + str(expected_compatibility_shape)
                + "."
            )

        weighted_modulus = (
            self.metric.weights * reference_modulus
        )
        reduced_stiffness = (
            compatibility.T
            @ (
                weighted_modulus[:, np.newaxis]
                * compatibility
            )
        )
        reduced_stiffness = 0.5 * (
            reduced_stiffness + reduced_stiffness.T
        )

        if np.any(~np.isfinite(reduced_stiffness)):
            raise FloatingPointError(
                "Reference reduced stiffness contains non-finite values."
            )

        try:
            cholesky = np.linalg.cholesky(reduced_stiffness)
        except np.linalg.LinAlgError as error:
            raise np.linalg.LinAlgError(
                "Reference tower stiffness on the free DOFs is not "
                "positive definite. Check boundary conditions, H, metric, "
                "and reference modulus."
            ) from error

        weighted_modulus = np.array(
            weighted_modulus,
            dtype=np.float64,
            copy=True,
        )
        reduced_stiffness = np.array(
            reduced_stiffness,
            dtype=np.float64,
            copy=True,
        )
        cholesky = np.array(
            cholesky,
            dtype=np.float64,
            copy=True,
        )
        weighted_modulus.setflags(write=False)
        reduced_stiffness.setflags(write=False)
        cholesky.setflags(write=False)

        object.__setattr__(self, "reference_modulus", reference_modulus)
        object.__setattr__(
            self,
            "compatibility_matrix",
            compatibility,
        )
        object.__setattr__(self, "free_dofs", free_dofs)
        object.__setattr__(self, "n_dof", n_dof)
        object.__setattr__(
            self,
            "_weighted_modulus",
            weighted_modulus,
        )
        object.__setattr__(
            self,
            "_reduced_stiffness",
            reduced_stiffness,
        )
        object.__setattr__(
            self,
            "_cholesky",
            cholesky,
        )

    @property
    def n_material_points(self) -> int:
        """Return the number of canonical material points."""
        return int(self.layout.n_material_points)

    @property
    def n_free_dofs(self) -> int:
        """Return the number of free structural DOFs."""
        return int(self.free_dofs.size)

    @property
    def reduced_stiffness(self) -> FloatArray:
        """Return a defensive copy of H^T M C0 H."""
        return self._reduced_stiffness.copy()

    def _solve_reduced(
        self,
        right_hand_side: FloatArray,
    ) -> FloatArray:
        """Solve the positive-definite reduced reference system."""
        rhs = np.asarray(right_hand_side, dtype=np.float64)
        y = np.linalg.solve(self._cholesky, rhs)
        return np.linalg.solve(self._cholesky.T, y)

    def apply_spatial(
        self,
        source_strain: FloatArray,
    ) -> EquilibriumProjectionTower:
        """Apply the reference structural projection to one spatial field."""
        source = np.asarray(source_strain, dtype=np.float64)
        expected_shape = (self.n_material_points,)
        if source.shape != expected_shape:
            raise ValueError(
                "source_strain must have shape "
                + str(expected_shape)
                + "."
            )
        if np.any(~np.isfinite(source)):
            raise ValueError(
                "source_strain must contain only finite values."
            )

        rhs = self.compatibility_matrix.T @ (
            self._weighted_modulus * source
        )
        displacement_free = self._solve_reduced(rhs)
        compatible = (
            self.compatibility_matrix @ displacement_free
        )
        stress = self.reference_modulus * (
            compatible - source
        )

        return EquilibriumProjectionTower(
            source_strain=source,
            compatible_strain=compatible,
            stress=stress,
            displacement_free=displacement_free,
        )

    def apply_history(
        self,
        source_strain: FloatArray,
    ) -> EquilibriumProjectionTower:
        """Apply the same spatial projection to every time row."""
        source = np.asarray(source_strain, dtype=np.float64)
        if source.ndim != 2:
            raise ValueError(
                "source_strain must have shape "
                "(n_time, n_material_points)."
            )
        if source.shape[0] < 1:
            raise ValueError(
                "source_strain must contain at least one time row."
            )
        if source.shape[1] != self.n_material_points:
            raise ValueError(
                "source_strain must have shape "
                "(n_time, n_material_points)."
            )
        if np.any(~np.isfinite(source)):
            raise ValueError(
                "source_strain must contain only finite values."
            )

        rhs = (
            source * self._weighted_modulus
        ) @ self.compatibility_matrix
        displacement_free = self._solve_reduced(
            rhs.T
        ).T
        compatible = (
            displacement_free
            @ self.compatibility_matrix.T
        )
        stress = (
            compatible - source
        ) * self.reference_modulus

        return EquilibriumProjectionTower(
            source_strain=source,
            compatible_strain=compatible,
            stress=stress,
            displacement_free=displacement_free,
        )

    def equilibrium_residual(
        self,
        stress: FloatArray,
    ) -> FloatArray:
        """
        Return H^T M sigma on free DOFs.

        For a history input, the returned shape is (n_time, n_free_dofs).
        """
        stress_array = np.asarray(stress, dtype=np.float64)

        if stress_array.ndim == 1:
            if stress_array.shape != (self.n_material_points,):
                raise ValueError(
                    "A spatial stress must have shape "
                    "(n_material_points,)."
                )
            if np.any(~np.isfinite(stress_array)):
                raise ValueError(
                    "stress must contain only finite values."
                )
            return self.compatibility_matrix.T @ (
                self.metric.weights * stress_array
            )

        if stress_array.ndim == 2:
            if stress_array.shape[1] != self.n_material_points:
                raise ValueError(
                    "A stress history must have shape "
                    "(n_time, n_material_points)."
                )
            if np.any(~np.isfinite(stress_array)):
                raise ValueError(
                    "stress must contain only finite values."
                )
            return (
                stress_array * self.metric.weights
            ) @ self.compatibility_matrix

        raise ValueError(
            "stress must be a spatial vector or time-by-material array."
        )


def _free_dofs_from_fixed(
    n_dof: int,
    fixed_dofs: Sequence[int],
) -> IntArray:
    """Return the sorted complement of fixed structural DOFs."""
    fixed = np.array(fixed_dofs, dtype=np.int64, copy=True)
    if fixed.ndim != 1:
        raise ValueError("fixed_dofs must be one-dimensional.")
    if np.any(fixed < 0) or np.any(fixed >= n_dof):
        raise ValueError(
            "fixed_dofs contain an out-of-range DOF index."
        )
    if np.unique(fixed).size != fixed.size:
        raise ValueError(
            "fixed_dofs must not contain duplicates."
        )

    mask = np.ones(n_dof, dtype=bool)
    mask[fixed] = False
    return np.flatnonzero(mask).astype(np.int64)


def build_tower_equilibrium_operator(
    system: Any,
    layout: Optional[MaterialPointLayout] = None,
    fixed_dofs: Optional[Sequence[int]] = None,
) -> TowerEquilibriumOperator:
    """
    Build the reference material-point metric and fiber compatibility operator
    from the current viscoplastic fiber-beam tower discretisation.

    The adapter intentionally relies on the frozen tower-system interface
    instead of taking ownership of the stateful FEM system.  Required members
    are:

        system.mesh.n_elements
        system.mesh.n_dof
        system.element_dofs
        system.elements
        system.material.E

    Each element supplies:

        jacobian
        gauss_weights
        b_matrices
        transformation
        sections[g].section.areas
        sections[g].section.y_coordinates

    The current tower v1 uses uniform n_gauss and n_fibers per Gauss section.
    """
    if not hasattr(system, "mesh"):
        raise TypeError("system must provide a mesh.")
    if not hasattr(system, "elements"):
        raise TypeError("system must provide beam elements.")
    if not hasattr(system, "element_dofs"):
        raise TypeError("system must provide element_dofs.")
    if not hasattr(system, "material"):
        raise TypeError("system must provide material data.")

    n_elements = int(system.mesh.n_elements)
    n_dof = int(system.mesh.n_dof)
    elements = tuple(system.elements)

    if len(elements) != n_elements:
        raise ValueError(
            "system.elements must contain one element per mesh element."
        )
    if n_elements < 1:
        raise ValueError(
            "The tower discretisation must contain at least one element."
        )

    first_element = elements[0]
    n_gauss = int(len(first_element.sections))
    if n_gauss < 1:
        raise ValueError(
            "Each tower element must contain at least one Gauss section."
        )
    first_section = first_element.sections[0].section
    n_fibers = int(first_section.n_fibers)
    if n_fibers < 1:
        raise ValueError(
            "Each Gauss section must contain at least one fiber."
        )

    if layout is None:
        layout = MaterialPointLayout(
            n_elements=n_elements,
            n_gauss=n_gauss,
            n_fibers=n_fibers,
        )
    elif not isinstance(layout, MaterialPointLayout):
        raise TypeError(
            "layout must be a MaterialPointLayout or None."
        )

    expected_layout = (
        n_elements,
        n_gauss,
        n_fibers,
    )
    actual_layout = (
        layout.n_elements,
        layout.n_gauss,
        layout.n_fibers,
    )
    if actual_layout != expected_layout:
        raise ValueError(
            "layout does not match the current tower discretisation."
        )

    if fixed_dofs is None:
        if n_dof < 4:
            raise ValueError(
                "A cantilever tower must contain free DOFs beyond the base."
            )
        fixed = np.array([0, 1, 2], dtype=np.int64)
    else:
        fixed = np.array(
            fixed_dofs,
            dtype=np.int64,
            copy=True,
        )
    free_dofs = _free_dofs_from_fixed(
        n_dof,
        fixed,
    )
    if free_dofs.size < 1:
        raise ValueError(
            "The tower projection requires at least one free DOF."
        )

    element_dofs = np.asarray(
        system.element_dofs,
        dtype=np.int64,
    )
    if element_dofs.shape != (n_elements, 6):
        raise ValueError(
            "system.element_dofs must have shape (n_elements, 6)."
        )
    if np.any(element_dofs < 0) or np.any(element_dofs >= n_dof):
        raise ValueError(
            "system.element_dofs contain an out-of-range DOF."
        )

    free_lookup = np.full(
        n_dof,
        -1,
        dtype=np.int64,
    )
    free_lookup[free_dofs] = np.arange(
        free_dofs.size,
        dtype=np.int64,
    )

    weights = np.empty(
        layout.n_material_points,
        dtype=np.float64,
    )
    compatibility = np.zeros(
        (
            layout.n_material_points,
            free_dofs.size,
        ),
        dtype=np.float64,
    )

    for element_index, element in enumerate(elements):
        sections = tuple(element.sections)
        gauss_weights = np.asarray(
            element.gauss_weights,
            dtype=np.float64,
        )
        b_matrices = np.asarray(
            element.b_matrices,
            dtype=np.float64,
        )
        transformation = np.asarray(
            element.transformation,
            dtype=np.float64,
        )
        jacobian = float(element.jacobian)

        if len(sections) != n_gauss:
            raise ValueError(
                "All tower elements must use the same n_gauss in tower v1."
            )
        if gauss_weights.shape != (n_gauss,):
            raise ValueError(
                "element.gauss_weights must have shape (n_gauss,)."
            )
        if b_matrices.shape != (n_gauss, 2, 6):
            raise ValueError(
                "element.b_matrices must have shape (n_gauss, 2, 6)."
            )
        if transformation.shape != (6, 6):
            raise ValueError(
                "element.transformation must have shape (6, 6)."
            )
        if (
            not np.isfinite(jacobian)
            or jacobian <= 0.0
        ):
            raise ValueError(
                "element.jacobian must be finite and positive."
            )
        if (
            np.any(~np.isfinite(gauss_weights))
            or np.any(gauss_weights <= 0.0)
        ):
            raise ValueError(
                "element.gauss_weights must be finite and positive."
            )

        dofs = element_dofs[element_index]

        for gauss_index, section_wrapper in enumerate(sections):
            section = section_wrapper.section
            areas = np.asarray(
                section.areas,
                dtype=np.float64,
            )
            y_coordinates = np.asarray(
                section.y_coordinates,
                dtype=np.float64,
            )

            if int(section.n_fibers) != n_fibers:
                raise ValueError(
                    "All Gauss sections must use the same fiber count "
                    "in tower v1."
                )
            if areas.shape != (n_fibers,):
                raise ValueError(
                    "section.areas must have shape (n_fibers,)."
                )
            if y_coordinates.shape != (n_fibers,):
                raise ValueError(
                    "section.y_coordinates must have shape (n_fibers,)."
                )
            if (
                np.any(~np.isfinite(areas))
                or np.any(areas <= 0.0)
            ):
                raise ValueError(
                    "All fiber areas must be finite and positive."
                )
            if np.any(~np.isfinite(y_coordinates)):
                raise ValueError(
                    "All fiber y-coordinates must be finite."
                )

            generalized_to_global = (
                b_matrices[gauss_index]
                @ transformation
            )

            for fiber_index in range(n_fibers):
                q = layout.flatten(
                    element_index,
                    gauss_index,
                    fiber_index,
                )
                weights[q] = (
                    jacobian
                    * gauss_weights[gauss_index]
                    * areas[fiber_index]
                )

                fiber_row_local_global = (
                    np.array(
                        [
                            1.0,
                            -float(y_coordinates[fiber_index]),
                        ],
                        dtype=np.float64,
                    )
                    @ generalized_to_global
                )

                for local_dof, global_dof in enumerate(dofs):
                    free_column = int(
                        free_lookup[int(global_dof)]
                    )
                    if free_column >= 0:
                        compatibility[q, free_column] += (
                            fiber_row_local_global[local_dof]
                        )

    reference_modulus_scalar = float(system.material.E)
    if (
        not np.isfinite(reference_modulus_scalar)
        or reference_modulus_scalar <= 0.0
    ):
        raise ValueError(
            "system.material.E must be finite and positive."
        )
    reference_modulus = np.full(
        layout.n_material_points,
        reference_modulus_scalar,
        dtype=np.float64,
    )

    metric = MaterialPointMetric(weights=weights)

    return TowerEquilibriumOperator(
        layout=layout,
        metric=metric,
        reference_modulus=reference_modulus,
        compatibility_matrix=compatibility,
        free_dofs=free_dofs,
        n_dof=n_dof,
    )

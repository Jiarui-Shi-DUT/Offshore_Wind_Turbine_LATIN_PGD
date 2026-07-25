# -*- coding: utf-8 -*-
"""
Residual-driven enrichment of the one-dimensional LATIN-PGD basis.

After the temporal functions of the existing PGD basis have been updated,
the remaining mechanical residual is reduced by adding one separated pair

    delta_eps_p_dot(x, t) = lambda_dot(t) * eps_p_bar(x),
    delta_sigma(x, t)     = lambda(t) * sigma_bar(x),

with

    sigma_bar = C(E - I) eps_p_bar.

The new spatial and temporal functions are obtained by an alternating
fixed-point procedure:

1. keep the temporal functions fixed and solve a weighted spatial problem;
2. keep the spatial function fixed and update the temporal functions;
3. repeat until the separated correction changes sufficiently little.

The spatial function is orthogonalized against the existing PGD basis and
normalized with the finite-element volume inner product.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from latin.equilibrium_operator import apply_equilibrium_operator
from latin.pgd_basis import PGDBasis1D, PGDMode1D
from latin.pgd_time_update import update_pgd_time_functions
from latin.search_directions import DescentSearchDirections
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class PGDEnrichmentResult:
    """One residual-driven PGD pair and its fixed-point diagnostics."""

    mode: PGDMode1D
    residual: FloatArray
    residual_norm_before: float
    residual_norm_after: float
    relative_residual: float
    fixed_point_history: FloatArray
    iterations: int
    converged: bool
    accepted: bool

    @property
    def residual_reduction(self) -> float:
        """Fraction of the initial residual removed by the new pair."""
        if self.residual_norm_before <= np.finfo(float).eps:
            return 0.0
        return float(
            1.0
            - self.residual_norm_after / self.residual_norm_before
        )


def _validate_inputs(
    basis: PGDBasis1D,
    time: FloatArray,
    residual: FloatArray,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    fixed_point_tolerance: float,
    max_fixed_point_iterations: int,
    minimum_spatial_norm: float,
    rcond: float,
) -> Tuple[FloatArray, FloatArray]:
    """Validate and normalise the enrichment inputs."""
    time_array = np.asarray(time, dtype=np.float64)
    residual_array = np.asarray(residual, dtype=np.float64)

    if time_array.ndim != 1:
        raise ValueError("time must be one-dimensional.")
    if time_array.size != basis.n_time:
        raise ValueError(
            "time and basis must contain the same number of points."
        )
    if np.any(np.diff(time_array) <= 0.0):
        raise ValueError("time must be strictly increasing.")
    if not np.all(np.isfinite(time_array)):
        raise ValueError("time contains non-finite values.")

    if residual_array.shape != basis.field_shape:
        raise ValueError(
            "residual must have shape (n_time, n_elements)."
        )
    if not np.all(np.isfinite(residual_array)):
        raise ValueError("residual contains non-finite values.")

    if directions.field_shape != basis.field_shape:
        raise ValueError(
            "directions and basis must use the same space-time grid."
        )
    if mesh.n_elements != basis.n_elements:
        raise ValueError(
            "mesh and basis must contain the same number of elements."
        )
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if fixed_point_tolerance <= 0.0:
        raise ValueError(
            "fixed_point_tolerance must be positive."
        )
    if max_fixed_point_iterations < 1:
        raise ValueError(
            "max_fixed_point_iterations must be at least one."
        )
    if minimum_spatial_norm <= 0.0:
        raise ValueError("minimum_spatial_norm must be positive.")
    if rcond <= 0.0:
        raise ValueError("rcond must be positive.")

    return time_array, residual_array


def _trapezoidal_weights(time: FloatArray) -> FloatArray:
    """Return nodal weights equivalent to trapezoidal time integration."""
    weights = np.zeros(time.size, dtype=np.float64)
    time_steps = np.diff(time)

    weights[0] = 0.5 * time_steps[0]
    weights[-1] = 0.5 * time_steps[-1]

    if time.size > 2:
        weights[1:-1] = 0.5 * (
            time_steps[:-1] + time_steps[1:]
        )

    return weights


def _weighted_space_time_inner_product(
    first: FloatArray,
    second: FloatArray,
    time: FloatArray,
    directions: DescentSearchDirections,
    element_volumes: FloatArray,
) -> float:
    """Evaluate the H_sigma^{-1}-weighted space-time inner product."""
    density = first * second / directions.H_sigma
    space_integral = density @ element_volumes
    value = float(np.trapz(space_integral, x=time))

    if not np.isfinite(value):
        raise FloatingPointError(
            "The PGD weighted inner product is non-finite."
        )

    return value


def _weighted_space_time_norm(
    field: FloatArray,
    time: FloatArray,
    directions: DescentSearchDirections,
    element_volumes: FloatArray,
) -> float:
    """Evaluate the H_sigma^{-1}-weighted space-time norm."""
    density = field**2 / directions.H_sigma
    space_integral = density @ element_volumes
    squared_norm = float(np.trapz(space_integral, x=time))

    if squared_norm < -1.0e-12:
        raise FloatingPointError(
            "The squared PGD residual norm became negative."
        )

    return float(np.sqrt(max(squared_norm, 0.0)))


def _equilibrium_stress_matrix(
    mesh: BarMesh1D,
    materials: Sequence[MaterialParameters],
) -> FloatArray:
    """
    Return the matrix mapping a spatial source strain to equilibrated stress.

    For the one-dimensional bar,

        sigma_bar = A_sigma eps_p_bar,

    and every row of A_sigma is identical because the axial correction stress
    is constant along the bar.
    """
    element_lengths = mesh.element_lengths
    elastic_moduli = np.asarray(
        [material.E for material in materials],
        dtype=np.float64,
    )

    if np.any(elastic_moduli <= 0.0):
        raise ValueError("All Young's moduli must be positive.")

    compliance_length = float(
        np.sum(element_lengths / elastic_moduli)
    )
    if compliance_length <= 0.0:
        raise FloatingPointError(
            "The assembled bar compliance is invalid."
        )

    stress_row = -element_lengths / compliance_length
    return np.tile(
        stress_row[np.newaxis, :],
        (mesh.n_elements, 1),
    )


def _orthogonalize_and_normalise(
    spatial_function: FloatArray,
    basis: PGDBasis1D,
    element_volumes: FloatArray,
    minimum_spatial_norm: float,
    rcond: float,
) -> Tuple[FloatArray, float]:
    """
    Remove the existing spatial-basis component and apply volume normalization.
    """
    vector = np.asarray(
        spatial_function,
        dtype=np.float64,
    ).copy()

    existing_basis = basis.spatial_plastic_strain_matrix()

    if existing_basis.shape[1] > 0:
        weighted_basis = (
            element_volumes[:, np.newaxis] * existing_basis
        )
        gram_matrix = existing_basis.T @ weighted_basis
        right_hand_side = (
            existing_basis.T
            @ (element_volumes * vector)
        )
        coefficients = np.linalg.lstsq(
            gram_matrix,
            right_hand_side,
            rcond=rcond,
        )[0]
        vector -= existing_basis @ coefficients

    norm = float(
        np.sqrt(
            max(
                np.dot(
                    element_volumes * vector,
                    vector,
                ),
                0.0,
            )
        )
    )

    if not np.isfinite(norm) or norm <= minimum_spatial_norm:
        raise FloatingPointError(
            "The candidate PGD spatial function is linearly dependent "
            "or numerically negligible."
        )

    return vector / norm, norm


def _initial_spatial_function(
    residual: FloatArray,
    directions: DescentSearchDirections,
    basis: PGDBasis1D,
    element_volumes: FloatArray,
    minimum_spatial_norm: float,
    rcond: float,
) -> FloatArray:
    """Build a deterministic initial spatial direction from the largest row."""
    row_energy = np.sum(
        residual**2
        * element_volumes[np.newaxis, :]
        / directions.H_sigma,
        axis=1,
    )
    seed_step = int(np.argmax(row_energy))

    if row_energy[seed_step] <= np.finfo(float).eps:
        raise FloatingPointError(
            "The residual is too small to generate a new PGD pair."
        )

    spatial_seed = -residual[seed_step, :]

    spatial_seed, _ = _orthogonalize_and_normalise(
        spatial_function=spatial_seed,
        basis=basis,
        element_volumes=element_volumes,
        minimum_spatial_norm=minimum_spatial_norm,
        rcond=rcond,
    )
    return spatial_seed


def _solve_spatial_function(
    temporal_amplitude: FloatArray,
    temporal_rate: FloatArray,
    residual: FloatArray,
    time: FloatArray,
    directions: DescentSearchDirections,
    stress_matrix: FloatArray,
    basis: PGDBasis1D,
    element_volumes: FloatArray,
    minimum_spatial_norm: float,
    rcond: float,
) -> FloatArray:
    """
    Solve the fixed-temporal weighted least-squares spatial problem.
    """
    n_elements = basis.n_elements
    time_weights = _trapezoidal_weights(time)

    normal_matrix = np.zeros(
        (n_elements, n_elements),
        dtype=np.float64,
    )
    right_hand_side = np.zeros(
        n_elements,
        dtype=np.float64,
    )
    identity = np.eye(n_elements, dtype=np.float64)

    for step in range(time.size):
        H_sigma = directions.H_sigma[step, :]
        operator = (
            temporal_rate[step] * identity
            - temporal_amplitude[step]
            * H_sigma[:, np.newaxis]
            * stress_matrix
        )
        weights = (
            time_weights[step]
            * element_volumes
            / H_sigma
        )

        weighted_operator = weights[:, np.newaxis] * operator
        normal_matrix += operator.T @ weighted_operator
        right_hand_side -= (
            operator.T @ (weights * residual[step, :])
        )

    spatial_function = np.linalg.lstsq(
        normal_matrix,
        right_hand_side,
        rcond=rcond,
    )[0]

    if not np.all(np.isfinite(spatial_function)):
        raise FloatingPointError(
            "The PGD spatial fixed-point solution is non-finite."
        )

    spatial_function, _ = _orthogonalize_and_normalise(
        spatial_function=spatial_function,
        basis=basis,
        element_volumes=element_volumes,
        minimum_spatial_norm=minimum_spatial_norm,
        rcond=rcond,
    )
    return spatial_function


def _mode_correction_residual(
    residual: FloatArray,
    mode: PGDMode1D,
    directions: DescentSearchDirections,
) -> FloatArray:
    """Add one separated correction to the current mechanical residual."""
    return (
        residual
        + mode.plastic_strain_rate_correction()
        - directions.H_sigma * mode.stress_correction()
    )


def enrich_pgd_basis_once(
    basis: PGDBasis1D,
    time: FloatArray,
    residual: FloatArray,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    *,
    iteration_added: int = 0,
    fixed_point_tolerance: float = 1.0e-6,
    max_fixed_point_iterations: int = 30,
    minimum_spatial_norm: float = 1.0e-14,
    acceptance_tolerance: float = 1.0e-8,
    rcond: float = 1.0e-12,
) -> PGDEnrichmentResult:
    """
    Generate one new PGD pair from the residual of the current reduced solution.

    The returned mode is not appended automatically.  This allows the caller
    to inspect ``accepted`` and the residual reduction before modifying the
    persistent PGD basis.
    """
    time_array, residual_array = _validate_inputs(
        basis=basis,
        time=time,
        residual=residual,
        directions=directions,
        mesh=mesh,
        area=area,
        materials=materials,
        fixed_point_tolerance=fixed_point_tolerance,
        max_fixed_point_iterations=max_fixed_point_iterations,
        minimum_spatial_norm=minimum_spatial_norm,
        rcond=rcond,
    )

    if acceptance_tolerance < 0.0:
        raise ValueError(
            "acceptance_tolerance must be non-negative."
        )

    element_volumes = area * mesh.element_lengths
    stress_matrix = _equilibrium_stress_matrix(
        mesh=mesh,
        materials=materials,
    )

    residual_norm_before = _weighted_space_time_norm(
        field=residual_array,
        time=time_array,
        directions=directions,
        element_volumes=element_volumes,
    )
    if residual_norm_before <= np.finfo(float).eps:
        raise FloatingPointError(
            "The residual is already zero; enrichment is unnecessary."
        )

    spatial_function = _initial_spatial_function(
        residual=residual_array,
        directions=directions,
        basis=basis,
        element_volumes=element_volumes,
        minimum_spatial_norm=minimum_spatial_norm,
        rcond=rcond,
    )

    previous_correction = None
    fixed_point_values = []
    converged = False
    final_mode = None

    for fixed_point_iteration in range(
        1,
        max_fixed_point_iterations + 1,
    ):
        projection = apply_equilibrium_operator(
            mesh=mesh,
            area=area,
            materials=materials,
            source_strain=spatial_function,
        )
        spatial_stress = projection.stress[0, :]

        temporary_basis = PGDBasis1D(
            n_elements=basis.n_elements,
            n_time=basis.n_time,
            modes=[
                PGDMode1D(
                    spatial_plastic_strain=spatial_function,
                    spatial_stress=spatial_stress,
                    temporal_amplitude=np.zeros(
                        basis.n_time,
                        dtype=np.float64,
                    ),
                    temporal_rate=np.zeros(
                        basis.n_time,
                        dtype=np.float64,
                    ),
                    iteration_added=iteration_added,
                )
            ],
        )

        time_result = update_pgd_time_functions(
            basis=temporary_basis,
            time=time_array,
            directions=directions,
            forcing=-residual_array,
            mesh=mesh,
            area=area,
            rcond=rcond,
        )
        temporal_mode = time_result.basis.modes[0]

        spatial_function = _solve_spatial_function(
            temporal_amplitude=temporal_mode.temporal_amplitude,
            temporal_rate=temporal_mode.temporal_rate,
            residual=residual_array,
            time=time_array,
            directions=directions,
            stress_matrix=stress_matrix,
            basis=basis,
            element_volumes=element_volumes,
            minimum_spatial_norm=minimum_spatial_norm,
            rcond=rcond,
        )

        projection = apply_equilibrium_operator(
            mesh=mesh,
            area=area,
            materials=materials,
            source_strain=spatial_function,
        )
        final_mode = PGDMode1D(
            spatial_plastic_strain=spatial_function,
            spatial_stress=projection.stress[0, :],
            temporal_amplitude=temporal_mode.temporal_amplitude,
            temporal_rate=temporal_mode.temporal_rate,
            iteration_added=iteration_added,
        )

        correction = (
            final_mode.plastic_strain_rate_correction()
            - directions.H_sigma
            * final_mode.stress_correction()
        )

        if previous_correction is None:
            fixed_point_change = np.inf
        else:
            change_norm = _weighted_space_time_norm(
                field=correction - previous_correction,
                time=time_array,
                directions=directions,
                element_volumes=element_volumes,
            )
            correction_norm = _weighted_space_time_norm(
                field=correction,
                time=time_array,
                directions=directions,
                element_volumes=element_volumes,
            )
            fixed_point_change = float(
                change_norm
                / max(
                    correction_norm,
                    float(np.finfo(np.float64).eps),
                )
            )

        fixed_point_values.append(fixed_point_change)
        previous_correction = correction.copy()

        if (
            fixed_point_iteration > 1
            and fixed_point_change <= fixed_point_tolerance
        ):
            converged = True
            break

    if final_mode is None:
        raise RuntimeError(
            "The PGD enrichment terminated before producing a mode."
        )

    # Recompute the temporal function for the final spatial vector so that
    # the returned pair is mutually consistent.
    final_basis = PGDBasis1D(
        n_elements=basis.n_elements,
        n_time=basis.n_time,
        modes=[
            PGDMode1D(
                spatial_plastic_strain=final_mode.spatial_plastic_strain,
                spatial_stress=final_mode.spatial_stress,
                temporal_amplitude=np.zeros(
                    basis.n_time,
                    dtype=np.float64,
                ),
                temporal_rate=np.zeros(
                    basis.n_time,
                    dtype=np.float64,
                ),
                iteration_added=iteration_added,
            )
        ],
    )
    final_time_result = update_pgd_time_functions(
        basis=final_basis,
        time=time_array,
        directions=directions,
        forcing=-residual_array,
        mesh=mesh,
        area=area,
        rcond=rcond,
    )
    final_mode = final_time_result.basis.modes[0]

    # The sequential backward-Euler temporal update does not minimise the
    # complete space-time residual in one coupled solve.  Apply one scalar
    # line search to the separated correction.  Because zero scaling is
    # admissible, this step guarantees that an accepted mode cannot increase
    # the weighted residual norm.
    mode_correction = (
        final_mode.plastic_strain_rate_correction()
        - directions.H_sigma
        * final_mode.stress_correction()
    )
    correction_norm_squared = _weighted_space_time_inner_product(
        first=mode_correction,
        second=mode_correction,
        time=time_array,
        directions=directions,
        element_volumes=element_volumes,
    )

    numerical_zero = float(np.finfo(np.float64).eps)
    if correction_norm_squared <= numerical_zero:
        optimal_scale = 0.0
    else:
        optimal_scale = -(
            _weighted_space_time_inner_product(
                first=residual_array,
                second=mode_correction,
                time=time_array,
                directions=directions,
                element_volumes=element_volumes,
            )
            / correction_norm_squared
        )

    if not np.isfinite(optimal_scale):
        raise FloatingPointError(
            "The PGD enrichment line-search scale is non-finite."
        )

    final_mode.temporal_amplitude *= optimal_scale
    final_mode.temporal_rate *= optimal_scale

    final_residual = _mode_correction_residual(
        residual=residual_array,
        mode=final_mode,
        directions=directions,
    )

    residual_norm_after = _weighted_space_time_norm(
        field=final_residual,
        time=time_array,
        directions=directions,
        element_volumes=element_volumes,
    )
    relative_residual = float(
        residual_norm_after / residual_norm_before
    )

    accepted = bool(
        np.isfinite(relative_residual)
        and residual_norm_after
        < residual_norm_before
        * (1.0 - acceptance_tolerance)
    )

    return PGDEnrichmentResult(
        mode=final_mode,
        residual=final_residual,
        residual_norm_before=residual_norm_before,
        residual_norm_after=residual_norm_after,
        relative_residual=relative_residual,
        fixed_point_history=np.asarray(
            fixed_point_values,
            dtype=np.float64,
        ),
        iterations=len(fixed_point_values),
        converged=converged,
        accepted=accepted,
    )

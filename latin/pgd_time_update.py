# -*- coding: utf-8 -*-
"""
Update of PGD temporal functions on a fixed spatial basis.

At the beginning of a LATIN global stage, the spatial PGD functions generated
during previous iterations are reused and only their temporal coefficients are
updated. In the one-dimensional setting, the reduced residual is

    r(x, t) =
        P(x) lambda_dot(t)
        - H_sigma(x, t) S(x) lambda(t)
        - f(x, t),

where

    P[:, j] = eps_p_bar_j,
    S[:, j] = C(E - I) eps_p_bar_j,

and f is the known forcing of the plastic global correction.

The temporal coefficients minimise the H_sigma^{-1}-weighted mechanical
residual. A backward-Euler discretisation is used for lambda_dot, which is
the discrete counterpart of the zero-order discontinuous Galerkin time
treatment described in the reference formulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from latin.pgd_basis import PGDBasis1D
from latin.search_directions import DescentSearchDirections


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class PGDTimeUpdateResult:
    """Updated basis and diagnostics of the reduced temporal problem."""

    basis: PGDBasis1D
    residual: FloatArray
    weighted_residual_norm: float
    relative_residual: float
    condition_history: FloatArray

    @property
    def n_modes(self) -> int:
        """Number of reused spatial PGD modes."""
        return self.basis.n_modes


def _validate_inputs(
    basis: PGDBasis1D,
    time: FloatArray,
    directions: DescentSearchDirections,
    forcing: FloatArray,
    mesh: BarMesh1D,
    area: float,
    rcond: float,
) -> Tuple[FloatArray, FloatArray]:
    """Validate and normalise the reduced temporal-update inputs."""
    time_array = np.asarray(time, dtype=np.float64)
    forcing_array = np.asarray(forcing, dtype=np.float64)

    if basis.n_modes < 1:
        raise ValueError(
            "At least one PGD mode is required to update time functions."
        )
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
    if forcing_array.shape != basis.field_shape:
        raise ValueError(
            "forcing must have shape (n_time, n_elements)."
        )
    if not np.all(np.isfinite(forcing_array)):
        raise ValueError("forcing contains non-finite values.")
    if directions.field_shape != basis.field_shape:
        raise ValueError(
            "directions and basis must use the same space-time grid."
        )
    if mesh.n_elements != basis.n_elements:
        raise ValueError(
            "mesh and basis must contain the same number of elements."
        )
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if rcond <= 0.0 or not np.isfinite(rcond):
        raise ValueError("rcond must be positive and finite.")

    return time_array, forcing_array


def _weighted_least_squares(
    matrix: FloatArray,
    right_hand_side: FloatArray,
    weights: FloatArray,
    rcond: float,
) -> Tuple[FloatArray, float]:
    """Solve one weighted least-squares problem and estimate conditioning."""
    square_root_weights = np.sqrt(weights)
    weighted_matrix = square_root_weights[:, np.newaxis] * matrix
    weighted_right_hand_side = (
        square_root_weights * right_hand_side
    )

    solution, _, _, singular_values = np.linalg.lstsq(
        weighted_matrix,
        weighted_right_hand_side,
        rcond=rcond,
    )

    if singular_values.size == 0:
        condition_number = np.inf
    elif singular_values[-1] <= np.finfo(float).eps:
        condition_number = np.inf
    else:
        condition_number = float(
            singular_values[0] / singular_values[-1]
        )

    if not np.all(np.isfinite(solution)):
        raise FloatingPointError(
            "The reduced temporal least-squares solution is non-finite."
        )

    return solution, condition_number


def _space_time_weighted_norm(
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
            "The squared reduced residual norm became negative."
        )

    return float(np.sqrt(max(squared_norm, 0.0)))


def update_pgd_time_functions(
    basis: PGDBasis1D,
    time: FloatArray,
    directions: DescentSearchDirections,
    forcing: FloatArray,
    mesh: BarMesh1D,
    area: float,
    *,
    rcond: float = 1.0e-12,
) -> PGDTimeUpdateResult:
    """
    Update all temporal PGD coefficients while keeping the spatial basis fixed.

    Parameters
    ----------
    basis:
        Existing PGD basis containing at least one spatial mode.
    time:
        LATIN time grid.
    directions:
        Current descent search-direction fields.
    forcing:
        Known right-hand side f in

            Delta eps_p_dot - H_sigma Delta sigma_plastic = f.

    mesh, area:
        Spatial integration data.
    rcond:
        Relative singular-value cut-off used by ``numpy.linalg.lstsq``.

    Returns
    -------
    PGDTimeUpdateResult
        A deep-copied basis with updated temporal amplitudes and rates,
        together with the complete residual field and conditioning diagnostics.
    """
    time_array, forcing_array = _validate_inputs(
        basis=basis,
        time=time,
        directions=directions,
        forcing=forcing,
        mesh=mesh,
        area=area,
        rcond=rcond,
    )

    spatial_plastic = basis.spatial_plastic_strain_matrix()
    spatial_stress = basis.spatial_stress_matrix()

    n_time = basis.n_time
    n_modes = basis.n_modes

    amplitudes = np.zeros(
        (n_time, n_modes),
        dtype=np.float64,
    )
    rates = np.zeros(
        (n_time, n_modes),
        dtype=np.float64,
    )
    condition_history = np.zeros(
        n_time,
        dtype=np.float64,
    )

    element_volumes = area * mesh.element_lengths

    initial_weights = (
        element_volumes / directions.H_sigma[0, :]
    )
    rates[0, :], condition_history[0] = _weighted_least_squares(
        matrix=spatial_plastic,
        right_hand_side=forcing_array[0, :],
        weights=initial_weights,
        rcond=rcond,
    )

    for step in range(1, n_time):
        time_step = float(
            time_array[step] - time_array[step - 1]
        )
        H_sigma = directions.H_sigma[step, :]

        reduced_matrix = (
            spatial_plastic / time_step
            - H_sigma[:, np.newaxis] * spatial_stress
        )
        right_hand_side = (
            forcing_array[step, :]
            + spatial_plastic
            @ (amplitudes[step - 1, :] / time_step)
        )
        weights = element_volumes / H_sigma

        amplitudes[step, :], condition_history[step] = (
            _weighted_least_squares(
                matrix=reduced_matrix,
                right_hand_side=right_hand_side,
                weights=weights,
                rcond=rcond,
            )
        )
        rates[step, :] = (
            amplitudes[step, :] - amplitudes[step - 1, :]
        ) / time_step

    updated_basis = basis.copy()

    for mode_index, mode in enumerate(updated_basis.modes):
        mode.temporal_amplitude = amplitudes[:, mode_index].copy()
        mode.temporal_rate = rates[:, mode_index].copy()

    reconstructed_rate = (
        updated_basis.plastic_strain_rate_correction()
    )
    reconstructed_stress = updated_basis.stress_correction()

    residual = (
        reconstructed_rate
        - directions.H_sigma * reconstructed_stress
        - forcing_array
    )

    residual_norm = _space_time_weighted_norm(
        field=residual,
        time=time_array,
        directions=directions,
        element_volumes=element_volumes,
    )
    forcing_norm = _space_time_weighted_norm(
        field=forcing_array,
        time=time_array,
        directions=directions,
        element_volumes=element_volumes,
    )

    if forcing_norm <= np.finfo(float).eps:
        relative_residual = (
            0.0
            if residual_norm <= np.finfo(float).eps
            else np.inf
        )
    else:
        relative_residual = float(residual_norm / forcing_norm)

    return PGDTimeUpdateResult(
        basis=updated_basis,
        residual=residual,
        weighted_residual_norm=residual_norm,
        relative_residual=relative_residual,
        condition_history=condition_history,
    )

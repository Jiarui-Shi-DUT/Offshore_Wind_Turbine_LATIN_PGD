# -*- coding: utf-8 -*-
"""
Fixed-spatial-basis temporal update for the tower LATIN-PGD global stage.

For an existing tower PGD basis

    P = [p_1, ..., p_m],
    S = [s_1, ..., s_m],

the reduced plastic correction is

    Delta eps_p_dot(t) = P lambda_dot(t),
    Delta sigma'(t) = S lambda(t).

The Eq. (58)-(59) mechanical residual is

    r(t) = P lambda_dot(t) - H_sigma(t) S lambda(t) - f(t),

and the temporal coordinates minimise the H_sigma^{-1} material-point metric

    integral r(t)^T M H_sigma(t)^(-1) r(t) dt.

Backward Euler is used for lambda_dot, matching the validated one-dimensional
implementation.  No enrichment, Gram-Schmidt, state construction, relaxation,
LATIN indicator, or transaction commit occurs in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from latin.pgd_basis import PGDBasisTower
from latin.tower_equilibrium_operator import (
    EquilibriumProjectionTower,
    MaterialPointMetric,
    TowerEquilibriumOperator,
)


FloatArray = NDArray[np.float64]


def _readonly_array(
    values: FloatArray,
    name: str,
    ndim: int,
) -> FloatArray:
    """Return a detached finite read-only float64 array."""
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != ndim:
        raise ValueError(
            name
            + " must have "
            + str(ndim)
            + " dimension(s)."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " contains non-finite values.")
    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class FixedBasisPGDResult:
    """Immutable result of the tower Eq. (58)-(59) temporal update."""

    basis: PGDBasisTower
    plastic_strain_correction: FloatArray
    plastic_strain_rate_correction: FloatArray
    plastic_projection: EquilibriumProjectionTower
    mechanical_residual: FloatArray
    weighted_residual_norm: float
    relative_residual: float
    forcing_norm: float
    reduced_converged: bool
    condition_history: FloatArray
    time_functions_updated: bool

    def __post_init__(self) -> None:
        if not isinstance(self.basis, PGDBasisTower):
            raise TypeError("basis must be a PGDBasisTower.")
        if not isinstance(
            self.plastic_projection,
            EquilibriumProjectionTower,
        ):
            raise TypeError(
                "plastic_projection must be an "
                "EquilibriumProjectionTower."
            )

        shape = self.basis.field_shape
        arrays = (
            (
                "plastic_strain_correction",
                self.plastic_strain_correction,
            ),
            (
                "plastic_strain_rate_correction",
                self.plastic_strain_rate_correction,
            ),
            (
                "mechanical_residual",
                self.mechanical_residual,
            ),
        )

        for name, values in arrays:
            array = _readonly_array(
                values,
                name,
                ndim=2,
            )
            if array.shape != shape:
                raise ValueError(
                    name
                    + " must have shape "
                    + str(shape)
                    + "."
                )
            object.__setattr__(self, name, array)

        if self.plastic_projection.source_strain.shape != shape:
            raise ValueError(
                "plastic_projection must use the basis field shape."
            )

        condition_history = _readonly_array(
            self.condition_history,
            "condition_history",
            ndim=1,
        )
        if condition_history.shape != (self.basis.n_time,):
            raise ValueError(
                "condition_history must contain one value per time point."
            )

        weighted_residual_norm = float(
            self.weighted_residual_norm
        )
        forcing_norm = float(self.forcing_norm)
        relative_residual = float(self.relative_residual)

        if (
            not np.isfinite(weighted_residual_norm)
            or weighted_residual_norm < 0.0
        ):
            raise ValueError(
                "weighted_residual_norm must be finite and non-negative."
            )
        if not np.isfinite(forcing_norm) or forcing_norm < 0.0:
            raise ValueError(
                "forcing_norm must be finite and non-negative."
            )
        if (
            not np.isfinite(relative_residual)
            and not np.isinf(relative_residual)
        ):
            raise ValueError(
                "relative_residual must be finite or infinity."
            )
        if relative_residual < 0.0:
            raise ValueError(
                "relative_residual must be non-negative."
            )
        if not isinstance(
            self.reduced_converged,
            (bool, np.bool_),
        ):
            raise TypeError("reduced_converged must be Boolean.")
        if not isinstance(
            self.time_functions_updated,
            (bool, np.bool_),
        ):
            raise TypeError(
                "time_functions_updated must be Boolean."
            )

        object.__setattr__(
            self,
            "basis",
            self.basis.copy(),
        )
        object.__setattr__(
            self,
            "weighted_residual_norm",
            weighted_residual_norm,
        )
        object.__setattr__(
            self,
            "relative_residual",
            relative_residual,
        )
        object.__setattr__(
            self,
            "forcing_norm",
            forcing_norm,
        )
        object.__setattr__(
            self,
            "reduced_converged",
            bool(self.reduced_converged),
        )
        object.__setattr__(
            self,
            "condition_history",
            condition_history,
        )
        object.__setattr__(
            self,
            "time_functions_updated",
            bool(self.time_functions_updated),
        )

    @property
    def n_modes(self) -> int:
        """Return the reduced-basis dimension after temporal updating."""
        return self.basis.n_modes


def _validate_inputs(
    basis: PGDBasisTower,
    time: FloatArray,
    forcing: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
    reduced_tolerance: float,
    rcond: float,
) -> Tuple[FloatArray, FloatArray, FloatArray]:
    """Validate fixed-basis tower temporal-update inputs."""
    if not isinstance(basis, PGDBasisTower):
        raise TypeError("basis must be a PGDBasisTower.")
    if not isinstance(metric, MaterialPointMetric):
        raise TypeError(
            "metric must be a MaterialPointMetric."
        )
    if not isinstance(
        equilibrium_operator,
        TowerEquilibriumOperator,
    ):
        raise TypeError(
            "equilibrium_operator must be a "
            "TowerEquilibriumOperator."
        )

    time_array = np.asarray(time, dtype=np.float64)
    forcing_array = np.asarray(
        forcing,
        dtype=np.float64,
    )
    H_sigma_array = np.asarray(
        H_sigma,
        dtype=np.float64,
    )

    if time_array.ndim != 1:
        raise ValueError("time must be one-dimensional.")
    if time_array.shape != (basis.n_time,):
        raise ValueError(
            "time and basis must contain the same number of points."
        )
    if np.any(~np.isfinite(time_array)):
        raise ValueError("time contains non-finite values.")
    if np.any(np.diff(time_array) <= 0.0):
        raise ValueError("time must be strictly increasing.")

    if forcing_array.shape != basis.field_shape:
        raise ValueError(
            "forcing must have shape "
            + str(basis.field_shape)
            + "."
        )
    if np.any(~np.isfinite(forcing_array)):
        raise ValueError(
            "forcing contains non-finite values."
        )

    if H_sigma_array.shape != basis.field_shape:
        raise ValueError(
            "H_sigma must have shape "
            + str(basis.field_shape)
            + "."
        )
    if np.any(~np.isfinite(H_sigma_array)):
        raise ValueError(
            "H_sigma contains non-finite values."
        )
    if np.any(H_sigma_array <= 0.0):
        raise ValueError(
            "H_sigma must be strictly positive."
        )

    if (
        metric.n_material_points
        != basis.n_material_points
    ):
        raise ValueError(
            "metric and basis must use the same material-point count."
        )
    if (
        equilibrium_operator.n_material_points
        != basis.n_material_points
    ):
        raise ValueError(
            "equilibrium_operator and basis must use the same "
            "material-point count."
        )
    if not np.array_equal(
        metric.weights,
        equilibrium_operator.metric.weights,
    ):
        raise ValueError(
            "metric must match the metric bound to "
            "equilibrium_operator."
        )

    reduced_tolerance = float(reduced_tolerance)
    rcond = float(rcond)
    if (
        not np.isfinite(reduced_tolerance)
        or reduced_tolerance <= 0.0
    ):
        raise ValueError(
            "reduced_tolerance must be positive and finite."
        )
    if not np.isfinite(rcond) or rcond <= 0.0:
        raise ValueError(
            "rcond must be positive and finite."
        )

    return time_array, forcing_array, H_sigma_array


def _validate_spatial_stress_pairs(
    basis: PGDBasisTower,
    equilibrium_operator: TowerEquilibriumOperator,
) -> None:
    """
    Verify that stored stress modes are associated equilibrium projections.

    This is a tower-v1 representation invariant rather than an additional PGD
    equation: s_j must equal C0(E-I)p_j for every stored spatial pair.
    """
    for mode_index, mode in enumerate(basis.modes):
        expected = equilibrium_operator.apply_spatial(
            mode.spatial_plastic_strain
        ).stress
        difference = float(
            np.linalg.norm(
                mode.spatial_stress - expected
            )
        )
        scale = max(
            1.0,
            float(np.linalg.norm(expected)),
        )
        if difference > 1.0e-10 * scale:
            raise ValueError(
                "Stored spatial_stress of tower PGD mode "
                + str(mode_index)
                + " is inconsistent with the reference "
                "equilibrium operator."
            )


def _weighted_least_squares(
    matrix: FloatArray,
    right_hand_side: FloatArray,
    weights: FloatArray,
    rcond: float,
) -> Tuple[FloatArray, float]:
    """Solve one H_sigma^{-1} material-point weighted least squares."""
    square_root_weights = np.sqrt(weights)
    weighted_matrix = (
        square_root_weights[:, np.newaxis]
        * matrix
    )
    weighted_right_hand_side = (
        square_root_weights * right_hand_side
    )

    solution, _, _, singular_values = np.linalg.lstsq(
        weighted_matrix,
        weighted_right_hand_side,
        rcond=rcond,
    )

    if singular_values.size == 0:
        condition_number = 0.0
    elif singular_values[-1] <= np.finfo(np.float64).eps:
        condition_number = float(
            np.finfo(np.float64).max
        )
    else:
        condition_number = float(
            singular_values[0]
            / singular_values[-1]
        )

    if np.any(~np.isfinite(solution)):
        raise FloatingPointError(
            "The tower reduced temporal solution is non-finite."
        )
    if not np.isfinite(condition_number):
        condition_number = float(
            np.finfo(np.float64).max
        )

    return (
        np.asarray(solution, dtype=np.float64),
        condition_number,
    )


def _space_time_weighted_norm(
    field: FloatArray,
    time: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
) -> float:
    """Evaluate the discrete H_sigma^{-1} tower space-time norm."""
    density = field**2 / H_sigma
    space_integral = density @ metric.weights
    if hasattr(np, "trapezoid"):
        squared_norm = float(
            np.trapezoid(
                space_integral,
                x=time,
            )
        )
    else:
        squared_norm = float(
            np.trapz(
                space_integral,
                x=time,
            )
        )

    tolerance = 1.0e-12 * max(
        1.0,
        float(np.max(np.abs(space_integral))),
    )
    if squared_norm < -tolerance:
        raise FloatingPointError(
            "The squared tower reduced residual norm became negative."
        )

    return float(
        np.sqrt(
            max(0.0, squared_norm)
        )
    )


def _relative_residual(
    residual_norm: float,
    forcing_norm: float,
) -> float:
    """Return the mechanical residual norm relative to the forcing norm."""
    numerical_zero = float(
        np.finfo(np.float64).eps
    )
    if forcing_norm <= numerical_zero:
        return (
            0.0
            if residual_norm <= numerical_zero
            else np.inf
        )
    return float(residual_norm / forcing_norm)


def update_tower_pgd_time_functions(
    basis: PGDBasisTower,
    time: FloatArray,
    forcing: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
    *,
    reduced_tolerance: float = 1.0e-4,
    rcond: float = 1.0e-12,
) -> FixedBasisPGDResult:
    """
    Solve the tower Eq. (58)-(59) temporal problem on a fixed spatial basis.

    Empty-basis semantics are valid: all reduced corrections are zero and the
    mechanical residual is exactly -forcing, allowing the outer solver to
    decide whether first-mode enrichment is required.
    """
    (
        time_array,
        forcing_array,
        H_sigma_array,
    ) = _validate_inputs(
        basis=basis,
        time=time,
        forcing=forcing,
        H_sigma=H_sigma,
        metric=metric,
        equilibrium_operator=equilibrium_operator,
        reduced_tolerance=reduced_tolerance,
        rcond=rcond,
    )

    _validate_spatial_stress_pairs(
        basis=basis,
        equilibrium_operator=equilibrium_operator,
    )

    if basis.n_modes == 0:
        zero_correction = np.zeros(
            basis.field_shape,
            dtype=np.float64,
        )
        plastic_projection = (
            equilibrium_operator.apply_history(
                zero_correction
            )
        )
        mechanical_residual = -forcing_array.copy()

        residual_norm = _space_time_weighted_norm(
            field=mechanical_residual,
            time=time_array,
            H_sigma=H_sigma_array,
            metric=metric,
        )
        forcing_norm = _space_time_weighted_norm(
            field=forcing_array,
            time=time_array,
            H_sigma=H_sigma_array,
            metric=metric,
        )
        relative_residual = _relative_residual(
            residual_norm=residual_norm,
            forcing_norm=forcing_norm,
        )

        return FixedBasisPGDResult(
            basis=basis,
            plastic_strain_correction=zero_correction,
            plastic_strain_rate_correction=zero_correction,
            plastic_projection=plastic_projection,
            mechanical_residual=mechanical_residual,
            weighted_residual_norm=residual_norm,
            relative_residual=relative_residual,
            forcing_norm=forcing_norm,
            reduced_converged=bool(
                np.isfinite(relative_residual)
                and relative_residual
                <= reduced_tolerance
            ),
            condition_history=np.zeros(
                basis.n_time,
                dtype=np.float64,
            ),
            time_functions_updated=False,
        )

    spatial_plastic = (
        basis.spatial_plastic_strain_matrix()
    )
    spatial_stress = (
        basis.spatial_stress_matrix()
    )

    amplitudes = np.zeros(
        (
            basis.n_time,
            basis.n_modes,
        ),
        dtype=np.float64,
    )
    rates = np.zeros_like(amplitudes)
    condition_history = np.zeros(
        basis.n_time,
        dtype=np.float64,
    )

    initial_weights = (
        metric.weights
        / H_sigma_array[0, :]
    )
    (
        rates[0, :],
        condition_history[0],
    ) = _weighted_least_squares(
        matrix=spatial_plastic,
        right_hand_side=forcing_array[0, :],
        weights=initial_weights,
        rcond=rcond,
    )

    for step in range(1, basis.n_time):
        time_step = float(
            time_array[step]
            - time_array[step - 1]
        )
        H_step = H_sigma_array[step, :]

        reduced_matrix = (
            spatial_plastic / time_step
            - H_step[:, np.newaxis]
            * spatial_stress
        )
        right_hand_side = (
            forcing_array[step, :]
            + spatial_plastic
            @ (
                amplitudes[step - 1, :]
                / time_step
            )
        )
        weights = metric.weights / H_step

        (
            amplitudes[step, :],
            condition_history[step],
        ) = _weighted_least_squares(
            matrix=reduced_matrix,
            right_hand_side=right_hand_side,
            weights=weights,
            rcond=rcond,
        )
        rates[step, :] = (
            amplitudes[step, :]
            - amplitudes[step - 1, :]
        ) / time_step

    updated_basis = basis.with_temporal_coordinates(
        temporal_amplitudes=amplitudes,
        temporal_rates=rates,
    )

    plastic_strain_correction = (
        updated_basis.plastic_strain_correction()
    )
    plastic_strain_rate_correction = (
        updated_basis.plastic_strain_rate_correction()
    )
    stored_stress_correction = (
        updated_basis.stress_correction()
    )
    plastic_projection = (
        equilibrium_operator.apply_history(
            plastic_strain_correction
        )
    )

    stress_difference = float(
        np.linalg.norm(
            stored_stress_correction
            - plastic_projection.stress
        )
    )
    stress_scale = max(
        1.0,
        float(
            np.linalg.norm(
                plastic_projection.stress
            )
        ),
    )
    if stress_difference > 1.0e-10 * stress_scale:
        raise FloatingPointError(
            "Reconstructed tower PGD stress correction is inconsistent "
            "with the reference equilibrium projection."
        )

    mechanical_residual = (
        plastic_strain_rate_correction
        - H_sigma_array
        * plastic_projection.stress
        - forcing_array
    )

    residual_norm = _space_time_weighted_norm(
        field=mechanical_residual,
        time=time_array,
        H_sigma=H_sigma_array,
        metric=metric,
    )
    forcing_norm = _space_time_weighted_norm(
        field=forcing_array,
        time=time_array,
        H_sigma=H_sigma_array,
        metric=metric,
    )
    relative_residual = _relative_residual(
        residual_norm=residual_norm,
        forcing_norm=forcing_norm,
    )

    return FixedBasisPGDResult(
        basis=updated_basis,
        plastic_strain_correction=plastic_strain_correction,
        plastic_strain_rate_correction=(
            plastic_strain_rate_correction
        ),
        plastic_projection=plastic_projection,
        mechanical_residual=mechanical_residual,
        weighted_residual_norm=residual_norm,
        relative_residual=relative_residual,
        forcing_norm=forcing_norm,
        reduced_converged=bool(
            np.isfinite(relative_residual)
            and relative_residual
            <= reduced_tolerance
        ),
        condition_history=condition_history,
        time_functions_updated=True,
    )

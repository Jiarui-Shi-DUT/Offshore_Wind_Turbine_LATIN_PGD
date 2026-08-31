# -*- coding: utf-8 -*-
"""
Tower-compatible nonlinear LATIN local stage on canonical material points.

The constitutive algorithm is intentionally the same as the validated
one-dimensional local stage.  The only structural generalisation is the space
index:

    element  ->  canonical material point q.

For the ascent-direction choice used by the project,

    (B_plus)^(-1) = 0,
    (b_plus)^(-1) = 0,

the thermodynamic-force histories remain fixed during the local projection:

    sigma_hat = sigma_i,
    beta_hat = beta_i,
    R_bar_hat = R_bar_i,
    Y_hat = Y_i.

The four integrated internal histories [eps_p, alpha, r_bar, D] are advanced
sequentially in time.  Distinct q-points are independent.

For homogeneous material-point parameters, OPT-4 advances all q-points in one
NumPy batch at each time step while preserving the sequential time history.
For heterogeneous material-point parameters, the validated pointwise path is
retained as a reference/fallback implementation.

The returned LatinStateTower is a new immutable value; the accepted input state
is never used as scratch storage.
"""

from __future__ import annotations

from typing import Sequence, Tuple, Union

import numpy as np

from latin.local_stage import (
    _integrate_one_local_step,
    local_rates_from_forces,
    unilateral_elastic_strain,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


MaterialInput = Union[
    MaterialParameters,
    Sequence[MaterialParameters],
]

_FLOAT_EPS = float(np.finfo(float).eps)


def _material_sequence(
    materials: MaterialInput,
    n_material_points: int,
) -> Tuple[MaterialParameters, ...]:
    """Broadcast one material or validate one material object per q-point."""
    if isinstance(materials, MaterialParameters):
        material_tuple = (materials,) * n_material_points
    else:
        material_tuple = tuple(materials)
        if len(material_tuple) == 1:
            material_tuple = material_tuple * n_material_points

    if len(material_tuple) != n_material_points:
        raise ValueError(
            "materials must contain one MaterialParameters object or one "
            "object per material point."
        )
    if any(
        not isinstance(material, MaterialParameters)
        for material in material_tuple
    ):
        raise TypeError(
            "materials must contain only MaterialParameters objects."
        )

    return material_tuple


def _materials_are_homogeneous(
    materials: Tuple[MaterialParameters, ...],
) -> bool:
    """Return True when all q-points have identical material parameters."""
    first = materials[0]
    return all(material == first for material in materials[1:])


def _vectorized_local_rates_from_forces(
    stress: np.ndarray,
    beta: np.ndarray,
    transformed_force: np.ndarray,
    energy_release_rate: np.ndarray,
    damage: np.ndarray,
    material: MaterialParameters,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Evaluate the local evolution laws for all homogeneous q-points at once.

    This is the array counterpart of latin.local_stage.local_rates_from_forces.
    The constitutive equations are unchanged; only the material-point execution
    granularity changes from scalar Python calls to NumPy array operations.
    """
    if material.R_inf <= 0.0:
        raise ValueError("R_inf must be positive.")
    if material.gamma <= 0.0:
        raise ValueError("gamma must be positive.")

    damage_safe = np.clip(
        np.asarray(damage, dtype=np.float64),
        0.0,
        material.damage_upper_bound,
    )
    one_minus_damage = 1.0 - damage_safe

    relative_effective_stress = (
        np.asarray(stress, dtype=np.float64)
        / one_minus_damage
        - np.asarray(beta, dtype=np.float64)
    )

    transformed_force_array = np.asarray(
        transformed_force,
        dtype=np.float64,
    )
    q_value = (
        transformed_force_array
        * material.sqrt_gamma
        / (2.0 * material.R_inf)
    )
    isotropic_force = (
        material.R_inf
        * q_value
        * (2.0 - q_value)
    )

    beta_array = np.asarray(beta, dtype=np.float64)
    yield_function = (
        np.abs(relative_effective_stress)
        + material.a * beta_array**2 / (2.0 * material.C)
        - isotropic_force
        - material.sigma_y
    )

    plastic_multiplier = (
        material.k_viscoplastic
        * np.maximum(yield_function, 0.0) ** material.n
    )

    flow_direction = np.where(
        np.abs(relative_effective_stress) <= _FLOAT_EPS,
        0.0,
        np.sign(relative_effective_stress),
    )

    plastic_strain_rate = (
        plastic_multiplier
        * flow_direction
        / one_minus_damage
    )

    alpha_rate = plastic_multiplier * (
        flow_direction
        - material.a * beta_array / material.C
    )

    r_bar_rate = plastic_multiplier * (
        material.sqrt_gamma
        - (
            material.gamma
            * transformed_force_array
            / (2.0 * material.R_inf)
        )
    )

    energy_array = np.asarray(
        energy_release_rate,
        dtype=np.float64,
    )
    damage_rate = (
        material.k_damage
        * np.maximum(
            energy_array - material.Y0,
            0.0,
        ) ** material.n_damage
    )

    return (
        np.asarray(plastic_strain_rate, dtype=np.float64),
        np.asarray(alpha_rate, dtype=np.float64),
        np.asarray(r_bar_rate, dtype=np.float64),
        np.asarray(damage_rate, dtype=np.float64),
    )


def _vectorized_unilateral_elastic_strain(
    stress: np.ndarray,
    damage: np.ndarray,
    material: MaterialParameters,
) -> np.ndarray:
    """
    Evaluate the unilateral elastic law for all homogeneous q-points at once.
    """
    stress_array = np.asarray(stress, dtype=np.float64)
    damage_safe = np.clip(
        np.asarray(damage, dtype=np.float64),
        0.0,
        material.damage_upper_bound,
    )

    stiffness_factor = np.where(
        stress_array >= 0.0,
        1.0 - damage_safe,
        1.0 - material.h * damage_safe,
    )

    if np.any(stiffness_factor <= 0.0):
        raise FloatingPointError(
            "The damaged elastic stiffness is non-positive."
        )

    return np.asarray(
        stress_array / (material.E * stiffness_factor),
        dtype=np.float64,
    )


def _integrate_homogeneous_material_points_one_step(
    state_old: np.ndarray,
    time_step: float,
    stress_old: np.ndarray,
    stress_new: np.ndarray,
    beta_old: np.ndarray,
    beta_new: np.ndarray,
    transformed_force_old: np.ndarray,
    transformed_force_new: np.ndarray,
    energy_old: np.ndarray,
    energy_new: np.ndarray,
    material: MaterialParameters,
) -> np.ndarray:
    """
    Advance all homogeneous q-points through one time step with classical RK4.

    The state shape is (4, Nq), corresponding to
    [eps_p, alpha, r_bar, D].  Time remains strictly sequential; only the
    independent q-direction is batched.
    """
    if time_step <= 0.0:
        raise ValueError("time_step must be positive.")

    state_start = np.asarray(
        state_old,
        dtype=np.float64,
    ).copy()

    if state_start.ndim != 2 or state_start.shape[0] != 4:
        raise ValueError(
            "Batched local internal state must have shape (4, Nq)."
        )

    # Preserve the OPT-3 interpolation structure at the three distinct RK4
    # abscissae: start, midpoint, end.  k2 and k3 share the midpoint forces.
    stress_start = (
        stress_old + 0.0 * (stress_new - stress_old)
    )
    stress_mid = (
        stress_old + 0.5 * (stress_new - stress_old)
    )
    stress_end = (
        stress_old + 1.0 * (stress_new - stress_old)
    )

    beta_start = (
        beta_old + 0.0 * (beta_new - beta_old)
    )
    beta_mid = (
        beta_old + 0.5 * (beta_new - beta_old)
    )
    beta_end = (
        beta_old + 1.0 * (beta_new - beta_old)
    )

    transformed_force_start = (
        transformed_force_old
        + 0.0 * (transformed_force_new - transformed_force_old)
    )
    transformed_force_mid = (
        transformed_force_old
        + 0.5 * (transformed_force_new - transformed_force_old)
    )
    transformed_force_end = (
        transformed_force_old
        + 1.0 * (transformed_force_new - transformed_force_old)
    )

    energy_start = (
        energy_old + 0.0 * (energy_new - energy_old)
    )
    energy_mid = (
        energy_old + 0.5 * (energy_new - energy_old)
    )
    energy_end = (
        energy_old + 1.0 * (energy_new - energy_old)
    )

    k1 = np.vstack(
        _vectorized_local_rates_from_forces(
            stress=stress_start,
            beta=beta_start,
            transformed_force=transformed_force_start,
            energy_release_rate=energy_start,
            damage=state_start[3, :],
            material=material,
        )
    )

    state_k2 = state_start + 0.5 * time_step * k1
    k2 = np.vstack(
        _vectorized_local_rates_from_forces(
            stress=stress_mid,
            beta=beta_mid,
            transformed_force=transformed_force_mid,
            energy_release_rate=energy_mid,
            damage=state_k2[3, :],
            material=material,
        )
    )

    state_k3 = state_start + 0.5 * time_step * k2
    k3 = np.vstack(
        _vectorized_local_rates_from_forces(
            stress=stress_mid,
            beta=beta_mid,
            transformed_force=transformed_force_mid,
            energy_release_rate=energy_mid,
            damage=state_k3[3, :],
            material=material,
        )
    )

    state_k4 = state_start + time_step * k3
    k4 = np.vstack(
        _vectorized_local_rates_from_forces(
            stress=stress_end,
            beta=beta_end,
            transformed_force=transformed_force_end,
            energy_release_rate=energy_end,
            damage=state_k4[3, :],
            material=material,
        )
    )

    state_new = state_start + (time_step / 6.0) * (
        k1 + 2.0 * k2 + 2.0 * k3 + k4
    )
    state_new = np.asarray(state_new, dtype=np.float64)
    state_new[3, :] = np.clip(
        state_new[3, :],
        0.0,
        material.damage_upper_bound,
    )

    if not np.all(np.isfinite(state_new)):
        raise FloatingPointError(
            "Non-finite internal variable in the LATIN local stage."
        )

    return state_new


def _solve_tower_local_stage_pointwise(
    global_state: LatinStateTower,
    material_tuple: Tuple[MaterialParameters, ...],
    fields: dict,
) -> None:
    """
    Validated pointwise reference/fallback path.

    This preserves the pre-OPT-4 execution structure for heterogeneous
    materials and for direct numerical cross-checking of the batched path.
    """
    for q, material in enumerate(material_tuple):
        internal_state = np.array(
            [
                global_state.plastic_strain[0, q],
                global_state.alpha[0, q],
                global_state.r_bar[0, q],
                global_state.damage[0, q],
            ],
            dtype=np.float64,
        )

        fields["plastic_strain"][0, q] = internal_state[0]
        fields["alpha"][0, q] = internal_state[1]
        fields["r_bar"][0, q] = internal_state[2]
        fields["damage"][0, q] = internal_state[3]

        initial_rates = local_rates_from_forces(
            stress=float(fields["stress"][0, q]),
            beta=float(fields["beta"][0, q]),
            transformed_force=float(fields["R_bar"][0, q]),
            energy_release_rate=float(
                fields["energy_release_rate"][0, q]
            ),
            damage=float(fields["damage"][0, q]),
            material=material,
        )

        (
            fields["plastic_strain_rate"][0, q],
            fields["alpha_rate"][0, q],
            fields["r_bar_rate"][0, q],
            fields["damage_rate"][0, q],
        ) = initial_rates

        fields["elastic_strain"][0, q] = unilateral_elastic_strain(
            stress=float(fields["stress"][0, q]),
            damage=float(fields["damage"][0, q]),
            material=material,
        )

        for step in range(1, global_state.n_time):
            time_step = float(
                global_state.time[step]
                - global_state.time[step - 1]
            )

            internal_state = _integrate_one_local_step(
                state_old=internal_state,
                time_step=time_step,
                stress_old=float(
                    fields["stress"][step - 1, q]
                ),
                stress_new=float(
                    fields["stress"][step, q]
                ),
                beta_old=float(
                    fields["beta"][step - 1, q]
                ),
                beta_new=float(
                    fields["beta"][step, q]
                ),
                transformed_force_old=float(
                    fields["R_bar"][step - 1, q]
                ),
                transformed_force_new=float(
                    fields["R_bar"][step, q]
                ),
                energy_old=float(
                    fields["energy_release_rate"][
                        step - 1,
                        q,
                    ]
                ),
                energy_new=float(
                    fields["energy_release_rate"][step, q]
                ),
                material=material,
            )

            fields["plastic_strain"][step, q] = internal_state[0]
            fields["alpha"][step, q] = internal_state[1]
            fields["r_bar"][step, q] = internal_state[2]
            fields["damage"][step, q] = internal_state[3]

            current_rates = local_rates_from_forces(
                stress=float(fields["stress"][step, q]),
                beta=float(fields["beta"][step, q]),
                transformed_force=float(
                    fields["R_bar"][step, q]
                ),
                energy_release_rate=float(
                    fields["energy_release_rate"][step, q]
                ),
                damage=float(fields["damage"][step, q]),
                material=material,
            )

            (
                fields["plastic_strain_rate"][step, q],
                fields["alpha_rate"][step, q],
                fields["r_bar_rate"][step, q],
                fields["damage_rate"][step, q],
            ) = current_rates

            fields["elastic_strain"][step, q] = (
                unilateral_elastic_strain(
                    stress=float(fields["stress"][step, q]),
                    damage=float(fields["damage"][step, q]),
                    material=material,
                )
            )


def _solve_tower_local_stage_homogeneous(
    global_state: LatinStateTower,
    material: MaterialParameters,
    fields: dict,
) -> None:
    """
    Batched homogeneous-material path.

    The q-direction is vectorized, while the time direction remains sequential
    because the internal variables carry history.
    """
    internal_state = np.vstack(
        (
            global_state.plastic_strain[0, :],
            global_state.alpha[0, :],
            global_state.r_bar[0, :],
            global_state.damage[0, :],
        )
    ).astype(np.float64, copy=True)

    fields["plastic_strain"][0, :] = internal_state[0, :]
    fields["alpha"][0, :] = internal_state[1, :]
    fields["r_bar"][0, :] = internal_state[2, :]
    fields["damage"][0, :] = internal_state[3, :]

    initial_rates = _vectorized_local_rates_from_forces(
        stress=fields["stress"][0, :],
        beta=fields["beta"][0, :],
        transformed_force=fields["R_bar"][0, :],
        energy_release_rate=fields["energy_release_rate"][0, :],
        damage=fields["damage"][0, :],
        material=material,
    )

    fields["plastic_strain_rate"][0, :] = initial_rates[0]
    fields["alpha_rate"][0, :] = initial_rates[1]
    fields["r_bar_rate"][0, :] = initial_rates[2]
    fields["damage_rate"][0, :] = initial_rates[3]

    fields["elastic_strain"][0, :] = (
        _vectorized_unilateral_elastic_strain(
            stress=fields["stress"][0, :],
            damage=fields["damage"][0, :],
            material=material,
        )
    )

    for step in range(1, global_state.n_time):
        time_step = float(
            global_state.time[step]
            - global_state.time[step - 1]
        )

        internal_state = _integrate_homogeneous_material_points_one_step(
            state_old=internal_state,
            time_step=time_step,
            stress_old=fields["stress"][step - 1, :],
            stress_new=fields["stress"][step, :],
            beta_old=fields["beta"][step - 1, :],
            beta_new=fields["beta"][step, :],
            transformed_force_old=fields["R_bar"][step - 1, :],
            transformed_force_new=fields["R_bar"][step, :],
            energy_old=fields["energy_release_rate"][step - 1, :],
            energy_new=fields["energy_release_rate"][step, :],
            material=material,
        )

        fields["plastic_strain"][step, :] = internal_state[0, :]
        fields["alpha"][step, :] = internal_state[1, :]
        fields["r_bar"][step, :] = internal_state[2, :]
        fields["damage"][step, :] = internal_state[3, :]

        current_rates = _vectorized_local_rates_from_forces(
            stress=fields["stress"][step, :],
            beta=fields["beta"][step, :],
            transformed_force=fields["R_bar"][step, :],
            energy_release_rate=fields["energy_release_rate"][step, :],
            damage=fields["damage"][step, :],
            material=material,
        )

        fields["plastic_strain_rate"][step, :] = current_rates[0]
        fields["alpha_rate"][step, :] = current_rates[1]
        fields["r_bar_rate"][step, :] = current_rates[2]
        fields["damage_rate"][step, :] = current_rates[3]

        fields["elastic_strain"][step, :] = (
            _vectorized_unilateral_elastic_strain(
                stress=fields["stress"][step, :],
                damage=fields["damage"][step, :],
                material=material,
            )
        )


def solve_tower_local_stage(
    global_state: LatinStateTower,
    materials: MaterialInput,
) -> LatinStateTower:
    """
    Project one immutable tower global state from A onto local manifold Gamma.

    Time remains sequential because the integrated constitutive variables carry
    history.  Distinct q-points are independent.  Homogeneous q-points use the
    OPT-4 batched path; heterogeneous q-points retain the validated pointwise
    reference/fallback path.
    """
    if not isinstance(global_state, LatinStateTower):
        raise TypeError(
            "global_state must be a LatinStateTower."
        )

    material_tuple = _material_sequence(
        materials=materials,
        n_material_points=global_state.n_material_points,
    )

    fields = {
        name: np.array(
            getattr(global_state, name),
            dtype=np.float64,
            copy=True,
        )
        for name in LatinStateTower.MATERIAL_FIELD_NAMES
    }

    # The ascent search direction keeps the thermodynamic forces fixed.
    fields["stress"][:, :] = global_state.stress
    fields["beta"][:, :] = global_state.beta
    fields["R_bar"][:, :] = global_state.R_bar
    fields["energy_release_rate"][:, :] = (
        global_state.energy_release_rate
    )

    if _materials_are_homogeneous(material_tuple):
        _solve_tower_local_stage_homogeneous(
            global_state=global_state,
            material=material_tuple[0],
            fields=fields,
        )
    else:
        _solve_tower_local_stage_pointwise(
            global_state=global_state,
            material_tuple=material_tuple,
            fields=fields,
        )

    return LatinStateTower(
        time=global_state.time,
        **fields,
    )

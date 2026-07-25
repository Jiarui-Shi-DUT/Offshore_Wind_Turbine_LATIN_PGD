# -*- coding: utf-8 -*-
"""
Nonlinear local stage of the one-dimensional LATIN solver.

At a known global state s_i in the admissibility manifold A, the local stage
computes a state s_hat_(i+1/2) in the constitutive manifold Gamma.

Following the ascent-direction choice used in the reference paper,

    (B_plus)^(-1) = 0
    (b_plus)^(-1) = 0

the thermodynamic-force histories remain fixed during the local projection:

    sigma_hat = sigma_i
    beta_hat = beta_i
    R_bar_hat = R_bar_i
    Y_hat = Y_i

The viscoplastic and damage evolution laws are then integrated independently
at every material point over the complete time interval. The nonlinear
unilateral elastic law is evaluated from the fixed stress history and the
updated damage history.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


def _positive_part(value: float) -> float:
    """Return max(value, 0)."""
    return float(max(value, 0.0))


def isotropic_force_from_transformed_force(
    transformed_force: float,
    material: MaterialParameters,
) -> float:
    """
    Recover the nonlinear isotropic hardening force R from R_bar.

    The partial-normal formulation uses

        R_bar = R_inf * r_bar

    together with the transformed relation from the paper.
    """
    if material.R_inf <= 0.0:
        raise ValueError("R_inf must be positive.")
    if material.gamma <= 0.0:
        raise ValueError("gamma must be positive.")

    q = (
        transformed_force
        * np.sqrt(material.gamma)
        / (2.0 * material.R_inf)
    )
    return float(material.R_inf * q * (2.0 - q))


def local_rates_from_forces(
    stress: float,
    beta: float,
    transformed_force: float,
    energy_release_rate: float,
    damage: float,
    material: MaterialParameters,
) -> Tuple[float, float, float, float]:
    """
    Evaluate the one-dimensional local evolution laws.

    Returns
    -------
    plastic_strain_rate, alpha_rate, r_bar_rate, damage_rate
    """
    damage_safe = float(
        np.clip(
            damage,
            0.0,
            material.damage_upper_bound,
        )
    )
    one_minus_damage = 1.0 - damage_safe

    relative_effective_stress = (
        stress / one_minus_damage - beta
    )

    isotropic_force = isotropic_force_from_transformed_force(
        transformed_force=transformed_force,
        material=material,
    )

    yield_function = (
        abs(relative_effective_stress)
        + material.a * beta**2 / (2.0 * material.C)
        - isotropic_force
        - material.sigma_y
    )

    plastic_multiplier = (
        material.k_viscoplastic
        * _positive_part(yield_function) ** material.n
    )

    if abs(relative_effective_stress) <= np.finfo(float).eps:
        flow_direction = 0.0
    else:
        flow_direction = float(np.sign(relative_effective_stress))

    plastic_strain_rate = (
        plastic_multiplier
        * flow_direction
        / one_minus_damage
    )

    alpha_rate = plastic_multiplier * (
        flow_direction
        - material.a * beta / material.C
    )

    r_bar_rate = plastic_multiplier * (
        np.sqrt(material.gamma)
        - (
            material.gamma
            * transformed_force
            / (2.0 * material.R_inf)
        )
    )

    damage_rate = (
        material.k_damage
        * _positive_part(
            energy_release_rate - material.Y0
        ) ** material.n_damage
    )

    return (
        float(plastic_strain_rate),
        float(alpha_rate),
        float(r_bar_rate),
        float(damage_rate),
    )


def unilateral_elastic_strain(
    stress: float,
    damage: float,
    material: MaterialParameters,
) -> float:
    """
    Evaluate the one-dimensional unilateral elastic state law.

    Tension:
        eps_e = sigma / [E (1 - D)]

    Compression:
        eps_e = sigma / [E (1 - h D)]
    """
    damage_safe = float(
        np.clip(
            damage,
            0.0,
            material.damage_upper_bound,
        )
    )

    if stress >= 0.0:
        stiffness_factor = 1.0 - damage_safe
    else:
        stiffness_factor = 1.0 - material.h * damage_safe

    if stiffness_factor <= 0.0:
        raise FloatingPointError(
            "The damaged elastic stiffness is non-positive."
        )

    return float(stress / (material.E * stiffness_factor))


def _interpolate(
    value_old: float,
    value_new: float,
    fraction: float,
) -> float:
    """Linearly interpolate a prescribed force history inside one step."""
    return float(value_old + fraction * (value_new - value_old))


def _local_state_rate(
    state: FloatArray,
    fraction: float,
    stress_old: float,
    stress_new: float,
    beta_old: float,
    beta_new: float,
    transformed_force_old: float,
    transformed_force_new: float,
    energy_old: float,
    energy_new: float,
    material: MaterialParameters,
) -> FloatArray:
    """Rate of [eps_p, alpha, r_bar, D] at one RK4 stage."""
    stress = _interpolate(stress_old, stress_new, fraction)
    beta = _interpolate(beta_old, beta_new, fraction)
    transformed_force = _interpolate(
        transformed_force_old,
        transformed_force_new,
        fraction,
    )
    energy = _interpolate(energy_old, energy_new, fraction)

    rates = local_rates_from_forces(
        stress=stress,
        beta=beta,
        transformed_force=transformed_force,
        energy_release_rate=energy,
        damage=float(state[3]),
        material=material,
    )
    return np.asarray(rates, dtype=np.float64)


def _integrate_one_local_step(
    state_old: FloatArray,
    time_step: float,
    stress_old: float,
    stress_new: float,
    beta_old: float,
    beta_new: float,
    transformed_force_old: float,
    transformed_force_new: float,
    energy_old: float,
    energy_new: float,
    material: MaterialParameters,
) -> FloatArray:
    """Advance one local material history step with classical RK4."""
    if time_step <= 0.0:
        raise ValueError("time_step must be positive.")

    state_start = np.asarray(
        state_old,
        dtype=np.float64,
    ).copy()

    if state_start.shape != (4,):
        raise ValueError("Local internal state must have shape (4,).")

    def rate(state: FloatArray, fraction: float) -> FloatArray:
        return _local_state_rate(
            state=state,
            fraction=fraction,
            stress_old=stress_old,
            stress_new=stress_new,
            beta_old=beta_old,
            beta_new=beta_new,
            transformed_force_old=transformed_force_old,
            transformed_force_new=transformed_force_new,
            energy_old=energy_old,
            energy_new=energy_new,
            material=material,
        )

    k1 = rate(state_start, 0.0)
    k2 = rate(
        state_start + 0.5 * time_step * k1,
        0.5,
    )
    k3 = rate(
        state_start + 0.5 * time_step * k2,
        0.5,
    )
    k4 = rate(
        state_start + time_step * k3,
        1.0,
    )

    state_new = state_start + (time_step / 6.0) * (
        k1 + 2.0 * k2 + 2.0 * k3 + k4
    )
    state_new = np.asarray(state_new, dtype=np.float64)
    state_new[3] = np.clip(
        state_new[3],
        0.0,
        material.damage_upper_bound,
    )

    if not np.all(np.isfinite(state_new)):
        raise FloatingPointError(
            "Non-finite internal variable in the LATIN local stage."
        )

    return state_new


def solve_local_stage(
    global_state: LatinState,
    materials: Sequence[MaterialParameters],
) -> LatinState:
    """
    Project one global LATIN state from A onto the local manifold Gamma.

    The operation is local in space: every element/material point is treated
    independently. Time remains sequential inside each material point because
    the integrated internal variables carry history.
    """
    if len(materials) != global_state.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )

    local_state = global_state.copy()

    # Infinite ascent directions keep all thermodynamic forces unchanged.
    local_state.stress[:, :] = global_state.stress
    local_state.beta[:, :] = global_state.beta
    local_state.R_bar[:, :] = global_state.R_bar
    local_state.energy_release_rate[:, :] = (
        global_state.energy_release_rate
    )

    for element, material in enumerate(materials):
        internal_state = np.array(
            [
                global_state.plastic_strain[0, element],
                global_state.alpha[0, element],
                global_state.r_bar[0, element],
                global_state.damage[0, element],
            ],
            dtype=np.float64,
        )

        local_state.plastic_strain[0, element] = internal_state[0]
        local_state.alpha[0, element] = internal_state[1]
        local_state.r_bar[0, element] = internal_state[2]
        local_state.damage[0, element] = internal_state[3]

        initial_rates = local_rates_from_forces(
            stress=float(local_state.stress[0, element]),
            beta=float(local_state.beta[0, element]),
            transformed_force=float(local_state.R_bar[0, element]),
            energy_release_rate=float(
                local_state.energy_release_rate[0, element]
            ),
            damage=float(local_state.damage[0, element]),
            material=material,
        )

        (
            local_state.plastic_strain_rate[0, element],
            local_state.alpha_rate[0, element],
            local_state.r_bar_rate[0, element],
            local_state.damage_rate[0, element],
        ) = initial_rates

        local_state.elastic_strain[0, element] = (
            unilateral_elastic_strain(
                stress=float(local_state.stress[0, element]),
                damage=float(local_state.damage[0, element]),
                material=material,
            )
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
                    local_state.stress[step - 1, element]
                ),
                stress_new=float(
                    local_state.stress[step, element]
                ),
                beta_old=float(
                    local_state.beta[step - 1, element]
                ),
                beta_new=float(
                    local_state.beta[step, element]
                ),
                transformed_force_old=float(
                    local_state.R_bar[step - 1, element]
                ),
                transformed_force_new=float(
                    local_state.R_bar[step, element]
                ),
                energy_old=float(
                    local_state.energy_release_rate[
                        step - 1,
                        element,
                    ]
                ),
                energy_new=float(
                    local_state.energy_release_rate[step, element]
                ),
                material=material,
            )

            local_state.plastic_strain[step, element] = (
                internal_state[0]
            )
            local_state.alpha[step, element] = internal_state[1]
            local_state.r_bar[step, element] = internal_state[2]
            local_state.damage[step, element] = internal_state[3]

            current_rates = local_rates_from_forces(
                stress=float(local_state.stress[step, element]),
                beta=float(local_state.beta[step, element]),
                transformed_force=float(
                    local_state.R_bar[step, element]
                ),
                energy_release_rate=float(
                    local_state.energy_release_rate[step, element]
                ),
                damage=float(local_state.damage[step, element]),
                material=material,
            )

            (
                local_state.plastic_strain_rate[step, element],
                local_state.alpha_rate[step, element],
                local_state.r_bar_rate[step, element],
                local_state.damage_rate[step, element],
            ) = current_rates

            local_state.elastic_strain[step, element] = (
                unilateral_elastic_strain(
                    stress=float(local_state.stress[step, element]),
                    damage=float(local_state.damage[step, element]),
                    material=material,
                )
            )

    return local_state

# -*- coding:utf-8 -*-
"""
作者：Shi Jiarui
日期：2026年07月24日
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from typing import Callable, Optional, Tuple, cast

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
StrainFunction = Callable[[float], float]

_FLOAT_EPS = float(np.finfo(float).eps)


@dataclass(frozen=True)
class MaterialParameters:
    """Parameters of the one-dimensional cyclic viscoplastic-damage model."""

    # Elasticity
    E: float = 134_000.0
    nu: float = 0.3

    # Isotropic hardening
    R_inf: float = 30.0
    gamma: float = 2.0

    # Kinematic hardening
    C: float = 5_500.0
    a: float = 250.0

    # Norton viscoplasticity
    K: float = 1_220.0
    n: float = 2.5

    # Damage
    k_damage: float = 2.778
    n_damage: float = 2.0
    h: float = 0.2

    # Yield stress
    sigma_y: float = 80.0

    # Numerical protection
    damage_upper_bound: float = 0.999

    @cached_property
    def k_viscoplastic(self) -> float:
        """Norton coefficient k = K^(-n), evaluated once per material object."""
        return float(self.K ** (-self.n))

    @cached_property
    def Y0(self) -> float:
        """Damage threshold Y0 = sigma_y^2 / (2E), evaluated once."""
        return float(self.sigma_y**2 / (2.0 * self.E))

    @cached_property
    def sqrt_gamma(self) -> float:
        """Square root of gamma, evaluated once per material object."""
        return float(np.sqrt(self.gamma))


@dataclass
class MaterialState:
    """Internal state variables of one material point."""

    plastic_strain: float = 0.0
    alpha: float = 0.0
    r_bar: float = 0.0
    damage: float = 0.0

    def to_array(self) -> FloatArray:
        return np.array(
            [
                self.plastic_strain,
                self.alpha,
                self.r_bar,
                self.damage,
            ],
            dtype=np.float64,
        )

    @classmethod
    def from_array(cls, values: FloatArray) -> "MaterialState":
        if values.shape != (4,):
            raise ValueError("Material state array must have shape (4,).")
        return cls(
            plastic_strain=float(values[0]),
            alpha=float(values[1]),
            r_bar=float(values[2]),
            damage=float(values[3]),
        )


@dataclass
class MaterialResponse:
    """Complete time history returned by a material-point simulation."""

    time: FloatArray
    strain: FloatArray
    elastic_strain: FloatArray
    stress: FloatArray
    plastic_strain: FloatArray
    plastic_strain_rate: FloatArray
    alpha: FloatArray
    alpha_rate: FloatArray
    beta: FloatArray
    r_bar: FloatArray
    r_bar_rate: FloatArray
    R_bar: FloatArray
    R: FloatArray
    damage: FloatArray
    damage_rate: FloatArray
    energy_release_rate: FloatArray
    effective_relative_stress: FloatArray
    yield_function: FloatArray
    plastic_multiplier: FloatArray


def positive_part(value: float) -> float:
    """Return the Macaulay positive part <value>+."""
    return max(float(value), 0.0)


def safe_damage(
    damage: float,
    material: MaterialParameters,
) -> float:
    """Clamp damage to the numerically admissible interval."""
    value = float(damage)
    upper_bound = float(material.damage_upper_bound)
    if value < 0.0:
        return 0.0
    if value > upper_bound:
        return upper_bound
    return value


def prescribed_strain(
    time: float,
    amplitude: float = 1.2e-3,
    period: float = 10.0,
) -> float:
    """Sinusoidal strain used by the reference material-point example."""
    return float(
        amplitude * np.sin(2.0 * np.pi * time / period)
    )


def stress_from_state(
    total_strain: float,
    plastic_strain: float,
    damage: float,
    material: MaterialParameters,
) -> float:
    """Compute one-dimensional stress with unilateral damaged elasticity."""
    damage_safe = safe_damage(damage, material)
    elastic_strain = float(total_strain - plastic_strain)

    if elastic_strain >= 0.0:
        effective_modulus = material.E * (1.0 - damage_safe)
    else:
        effective_modulus = material.E * (
            1.0 - material.h * damage_safe
        )

    return float(effective_modulus * elastic_strain)


def isotropic_hardening_force(
    r_bar: float,
    material: MaterialParameters,
) -> Tuple[float, float]:
    """Return transformed force R_bar and physical hardening force R."""
    R_bar = float(material.R_inf * r_bar)
    eta = float(0.5 * material.sqrt_gamma * r_bar)
    R = float(material.R_inf * eta * (2.0 - eta))
    return R_bar, R


def damage_energy_release_rate(
    stress: float,
    damage: float,
    material: MaterialParameters,
) -> float:
    """Compute the unilateral damage energy-release rate."""
    damage_safe = safe_damage(damage, material)

    if stress >= 0.0:
        denominator = (
            2.0
            * material.E
            * (1.0 - damage_safe) ** 2
        )
        return float(stress**2 / denominator)

    denominator = (
        2.0
        * material.E
        * (1.0 - material.h * damage_safe) ** 2
    )
    return float(material.h * stress**2 / denominator)


def evaluate_state(
    total_strain: float,
    state: FloatArray,
    material: MaterialParameters,
) -> Tuple[float, float, float, float, float, float, float]:
    """Evaluate stress, hardening forces, yield function and damage force."""
    if state.shape != (4,):
        raise ValueError("Material state array must have shape (4,).")

    plastic_strain = float(state[0])
    alpha = float(state[1])
    r_bar = float(state[2])
    damage = safe_damage(float(state[3]), material)

    stress = stress_from_state(
        total_strain,
        plastic_strain,
        damage,
        material,
    )

    beta = float(material.C * alpha)
    R_bar, R = isotropic_hardening_force(r_bar, material)

    effective_relative_stress = float(
        stress / (1.0 - damage) - beta
    )

    yield_function = float(
        abs(effective_relative_stress)
        + material.a * beta**2 / (2.0 * material.C)
        - R
        - material.sigma_y
    )

    energy_release_rate = damage_energy_release_rate(
        stress,
        damage,
        material,
    )

    return (
        stress,
        beta,
        R_bar,
        R,
        effective_relative_stress,
        yield_function,
        energy_release_rate,
    )


def _scalar_state_rate_components(
    total_strain: float,
    plastic_strain: float,
    alpha: float,
    r_bar: float,
    damage: float,
    material: MaterialParameters,
) -> Tuple[float, float, float, float]:
    # Scalar form of the same constitutive rate equations.
    damage = safe_damage(damage, material)

    elastic_strain = float(total_strain - plastic_strain)
    if elastic_strain >= 0.0:
        effective_modulus = material.E * (1.0 - damage)
    else:
        effective_modulus = material.E * (
            1.0 - material.h * damage
        )
    stress = float(effective_modulus * elastic_strain)

    beta = float(material.C * alpha)
    eta = float(0.5 * material.sqrt_gamma * r_bar)
    hardening_force = float(
        material.R_inf * eta * (2.0 - eta)
    )

    effective_relative_stress = float(
        stress / (1.0 - damage) - beta
    )
    yield_function = float(
        abs(effective_relative_stress)
        + material.a * beta**2 / (2.0 * material.C)
        - hardening_force
        - material.sigma_y
    )

    if stress >= 0.0:
        denominator = (
            2.0 * material.E * (1.0 - damage) ** 2
        )
        energy_release_rate = float(
            stress**2 / denominator
        )
    else:
        denominator = (
            2.0
            * material.E
            * (1.0 - material.h * damage) ** 2
        )
        energy_release_rate = float(
            material.h * stress**2 / denominator
        )

    positive_yield = (
        yield_function if yield_function > 0.0 else 0.0
    )
    plastic_multiplier = float(
        material.k_viscoplastic
        * positive_yield ** material.n
    )

    if abs(effective_relative_stress) <= _FLOAT_EPS:
        flow_direction = 0.0
    elif effective_relative_stress > 0.0:
        flow_direction = 1.0
    else:
        flow_direction = -1.0

    plastic_strain_rate = float(
        plastic_multiplier
        * flow_direction
        / (1.0 - damage)
    )

    alpha_rate = float(
        plastic_multiplier
        * (
            flow_direction
            - material.a * beta / material.C
        )
    )

    r_bar_rate = float(
        plastic_multiplier
        * (
            material.sqrt_gamma
            - 0.5 * material.gamma * r_bar
        )
    )

    damage_drive = energy_release_rate - material.Y0
    positive_damage_drive = (
        damage_drive if damage_drive > 0.0 else 0.0
    )
    damage_rate = float(
        material.k_damage
        * positive_damage_drive ** material.n_damage
    )

    return (
        plastic_strain_rate,
        alpha_rate,
        r_bar_rate,
        damage_rate,
    )


def state_rate(
    total_strain: float,
    state: FloatArray,
    material: MaterialParameters,
) -> FloatArray:
    """Compute rates [dot_eps_p, dot_alpha, dot_r_bar, dot_D]."""
    damage = safe_damage(float(state[3]), material)

    (
        _stress,
        beta,
        _R_bar,
        _R,
        effective_relative_stress,
        yield_function,
        energy_release_rate,
    ) = evaluate_state(total_strain, state, material)

    plastic_multiplier = float(
        material.k_viscoplastic
        * positive_part(yield_function) ** material.n
    )

    if abs(effective_relative_stress) <= np.finfo(float).eps:
        flow_direction = 0.0
    else:
        flow_direction = float(np.sign(effective_relative_stress))

    plastic_strain_rate = float(
        plastic_multiplier
        * flow_direction
        / (1.0 - damage)
    )

    alpha_rate = float(
        plastic_multiplier
        * (
            flow_direction
            - material.a * beta / material.C
        )
    )

    r_bar_rate = float(
        plastic_multiplier
        * (
            material.sqrt_gamma
            - 0.5 * material.gamma * float(state[2])
        )
    )

    damage_rate = float(
        material.k_damage
        * positive_part(
            energy_release_rate - material.Y0
        ) ** material.n_damage
    )

    return np.array(
        [
            plastic_strain_rate,
            alpha_rate,
            r_bar_rate,
            damage_rate,
        ],
        dtype=np.float64,
    )


def rk4_step(
    time: float,
    state: FloatArray,
    time_step: float,
    strain_function: StrainFunction,
    material: MaterialParameters,
) -> FloatArray:
    # Classical RK4 with scalar state components and identical equations.
    if time_step <= 0.0:
        raise ValueError("time_step must be positive.")

    dt = float(time_step)
    half_dt = 0.5 * dt

    eps_p0 = float(state[0])
    alpha0 = float(state[1])
    r_bar0 = float(state[2])
    damage0 = float(state[3])

    strain1 = strain_function(time)
    k1_eps_p, k1_alpha, k1_r_bar, k1_damage = (
        _scalar_state_rate_components(
            strain1,
            eps_p0,
            alpha0,
            r_bar0,
            damage0,
            material,
        )
    )

    strain2 = strain_function(time + half_dt)
    k2_eps_p, k2_alpha, k2_r_bar, k2_damage = (
        _scalar_state_rate_components(
            strain2,
            eps_p0 + half_dt * k1_eps_p,
            alpha0 + half_dt * k1_alpha,
            r_bar0 + half_dt * k1_r_bar,
            damage0 + half_dt * k1_damage,
            material,
        )
    )

    k3_eps_p, k3_alpha, k3_r_bar, k3_damage = (
        _scalar_state_rate_components(
            strain2,
            eps_p0 + half_dt * k2_eps_p,
            alpha0 + half_dt * k2_alpha,
            r_bar0 + half_dt * k2_r_bar,
            damage0 + half_dt * k2_damage,
            material,
        )
    )

    strain4 = strain_function(time + dt)
    k4_eps_p, k4_alpha, k4_r_bar, k4_damage = (
        _scalar_state_rate_components(
            strain4,
            eps_p0 + dt * k3_eps_p,
            alpha0 + dt * k3_alpha,
            r_bar0 + dt * k3_r_bar,
            damage0 + dt * k3_damage,
            material,
        )
    )

    sixth_dt = dt / 6.0
    new_state = np.array(
        [
            eps_p0 + sixth_dt * (
                k1_eps_p
                + 2.0 * k2_eps_p
                + 2.0 * k3_eps_p
                + k4_eps_p
            ),
            alpha0 + sixth_dt * (
                k1_alpha
                + 2.0 * k2_alpha
                + 2.0 * k3_alpha
                + k4_alpha
            ),
            r_bar0 + sixth_dt * (
                k1_r_bar
                + 2.0 * k2_r_bar
                + 2.0 * k3_r_bar
                + k4_r_bar
            ),
            damage0 + sixth_dt * (
                k1_damage
                + 2.0 * k2_damage
                + 2.0 * k3_damage
                + k4_damage
            ),
        ],
        dtype=np.float64,
    )

    new_state[3] = safe_damage(
        float(new_state[3]),
        material,
    )

    if not np.all(np.isfinite(new_state)):
        raise FloatingPointError(
            "Non-finite state encountered during RK4 integration."
        )

    return new_state


def run_material_point(
    material: MaterialParameters,
    total_time: float = 200.0,
    time_step: float = 0.1,
    strain_function: StrainFunction = prescribed_strain,
    initial_state: Optional[MaterialState] = None,
) -> MaterialResponse:
    """Integrate the complete history of one material point."""
    if total_time <= 0.0:
        raise ValueError("total_time must be positive.")
    if time_step <= 0.0:
        raise ValueError("time_step must be positive.")

    number_of_steps = int(round(total_time / time_step)) + 1
    time = np.linspace(
        0.0,
        total_time,
        number_of_steps,
        dtype=np.float64,
    )

    states = np.zeros((number_of_steps, 4), dtype=np.float64)
    if initial_state is None:
        initial_state = MaterialState()
    states[0] = initial_state.to_array()

    strain = np.zeros(number_of_steps, dtype=np.float64)
    elastic_strain = np.zeros(number_of_steps, dtype=np.float64)
    stress = np.zeros(number_of_steps, dtype=np.float64)
    plastic_strain_rate = np.zeros(number_of_steps, dtype=np.float64)
    alpha_rate = np.zeros(number_of_steps, dtype=np.float64)
    beta = np.zeros(number_of_steps, dtype=np.float64)
    r_bar_rate = np.zeros(number_of_steps, dtype=np.float64)
    R_bar = np.zeros(number_of_steps, dtype=np.float64)
    R = np.zeros(number_of_steps, dtype=np.float64)
    damage_rate = np.zeros(number_of_steps, dtype=np.float64)
    energy_release_rate = np.zeros(number_of_steps, dtype=np.float64)
    effective_relative_stress = np.zeros(number_of_steps, dtype=np.float64)
    yield_function = np.zeros(number_of_steps, dtype=np.float64)
    plastic_multiplier = np.zeros(number_of_steps, dtype=np.float64)

    for index in range(number_of_steps):
        current_time = float(time[index])
        current_state = states[index]

        strain_value = cast(float, strain_function(current_time))
        strain[index] = strain_value
        elastic_strain[index] = float(strain_value - current_state[0])

        (
            stress[index],
            beta[index],
            R_bar[index],
            R[index],
            effective_relative_stress[index],
            yield_function[index],
            energy_release_rate[index],
        ) = evaluate_state(strain_value, current_state, material)

        rates = state_rate(strain_value, current_state, material)
        plastic_strain_rate[index] = float(rates[0])
        alpha_rate[index] = float(rates[1])
        r_bar_rate[index] = float(rates[2])
        damage_rate[index] = float(rates[3])

        plastic_multiplier[index] = float(
            material.k_viscoplastic
            * positive_part(yield_function[index]) ** material.n
        )

        if index < number_of_steps - 1:
            states[index + 1] = rk4_step(
                current_time,
                current_state,
                time_step,
                strain_function,
                material,
            )

    return MaterialResponse(
        time=time,
        strain=strain,
        elastic_strain=elastic_strain,
        stress=stress,
        plastic_strain=states[:, 0].copy(),
        plastic_strain_rate=plastic_strain_rate,
        alpha=states[:, 1].copy(),
        alpha_rate=alpha_rate,
        beta=beta,
        r_bar=states[:, 2].copy(),
        r_bar_rate=r_bar_rate,
        R_bar=R_bar,
        R=R,
        damage=states[:, 3].copy(),
        damage_rate=damage_rate,
        energy_release_rate=energy_release_rate,
        effective_relative_stress=effective_relative_stress,
        yield_function=yield_function,
        plastic_multiplier=plastic_multiplier,
    )

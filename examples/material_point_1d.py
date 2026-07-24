
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, List, Tuple, cast

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
StrainFunction = Callable[[float], float]
PlotItem = Tuple[str, FloatArray, str]


@dataclass(frozen=True)
class MaterialParameters:
    """
    One-dimensional material parameters for the cyclic
    viscoplastic-damage model used in Section 5.1.
    """

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

    @property
    def k_viscoplastic(self) -> float:
        """Norton coefficient k = K^(-n)."""
        return float(self.K ** (-self.n))

    @property
    def Y0(self) -> float:
        """Damage threshold Y0 = sigma_y^2 / (2E)."""
        return float(self.sigma_y**2 / (2.0 * self.E))


@dataclass
class MaterialResponse:
    """History variables returned by the material-point simulation."""

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
    """Macaulay bracket <x>_+."""
    return max(float(value), 0.0)


def safe_damage(
    damage: float,
    material: MaterialParameters,
) -> float:
    """Clamp damage only for numerical protection."""
    return float(
        np.clip(
            damage,
            0.0,
            material.damage_upper_bound,
        )
    )


def prescribed_strain(
    time: float,
    amplitude: float = 1.2e-3,
    period: float = 10.0,
) -> float:
    """
    Prescribed sinusoidal strain.

    For the three-material bar:
        displacement amplitude = 1.2e-3 * L
    so the corresponding nominal strain amplitude is 1.2e-3.
    """
    return float(
        amplitude * np.sin(2.0 * np.pi * time / period)
    )


def stress_from_state(
    total_strain: float,
    plastic_strain: float,
    damage: float,
    material: MaterialParameters,
) -> float:
    """
    One-dimensional unilateral damaged elasticity.

    Tension:
        sigma = E * (1 - D) * eps_e

    Compression:
        sigma = E * (1 - hD) * eps_e
    """
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
    """
    Return transformed force R_bar and physical isotropic
    hardening force R.

    R_bar = R_inf * r_bar

    eta = sqrt(gamma) * r_bar / 2
    R = R_inf * eta * (2 - eta)
    """
    R_bar = float(material.R_inf * r_bar)

    eta = float(
        0.5 * np.sqrt(material.gamma) * r_bar
    )
    R = float(
        material.R_inf * eta * (2.0 - eta)
    )

    return R_bar, R


def damage_energy_release_rate(
    stress: float,
    damage: float,
    material: MaterialParameters,
) -> float:
    """
    One-dimensional unilateral damage energy release rate.

    Tension:
        Y = sigma^2 / [2E(1-D)^2]

    Compression:
        Y = h*sigma^2 / [2E(1-hD)^2]
    """
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
    return float(
        material.h * stress**2 / denominator
    )


def evaluate_state(
    total_strain: float,
    state: FloatArray,
    material: MaterialParameters,
) -> Tuple[
    float,
    float,
    float,
    float,
    float,
    float,
    float,
]:
    """
    Evaluate the current constitutive quantities.

    State vector:
        state[0] = eps_p
        state[1] = alpha
        state[2] = r_bar
        state[3] = D

    Returns:
        stress,
        beta,
        R_bar,
        R,
        q,
        f_p,
        Y
    """
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
    R_bar, R = isotropic_hardening_force(
        r_bar,
        material,
    )

    q = float(
        stress / (1.0 - damage) - beta
    )

    yield_function = float(
        abs(q)
        + material.a * beta**2
        / (2.0 * material.C)
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
        q,
        yield_function,
        energy_release_rate,
    )


def state_rate(
    total_strain: float,
    state: FloatArray,
    material: MaterialParameters,
) -> FloatArray:
    """
    Compute rates of the internal variables.

    State:
        [eps_p, alpha, r_bar, D]

    Rate:
        [dot_eps_p, dot_alpha, dot_r_bar, dot_D]
    """
    damage = safe_damage(float(state[3]), material)

    (
        _stress,
        beta,
        _R_bar,
        _R,
        q,
        yield_function,
        energy_release_rate,
    ) = evaluate_state(
        total_strain,
        state,
        material,
    )

    plastic_multiplier = float(
        material.k_viscoplastic
        * positive_part(yield_function) ** material.n
    )

    if abs(q) <= np.finfo(float).eps:
        flow_direction = 0.0
    else:
        flow_direction = float(np.sign(q))

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
            np.sqrt(material.gamma)
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
    """Advance one fixed time step with classical RK4."""

    k1 = state_rate(
        strain_function(time),
        state,
        material,
    )

    k2 = state_rate(
        strain_function(time + 0.5 * time_step),
        state + 0.5 * time_step * k1,
        material,
    )

    k3 = state_rate(
        strain_function(time + 0.5 * time_step),
        state + 0.5 * time_step * k2,
        material,
    )

    k4 = state_rate(
        strain_function(time + time_step),
        state + time_step * k3,
        material,
    )

    new_state = state + (
        time_step / 6.0
    ) * (
        k1
        + 2.0 * k2
        + 2.0 * k3
        + k4
    )

    new_state = np.asarray(
        new_state,
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
) -> MaterialResponse:
    """Run the full material-point history."""

    if total_time <= 0.0:
        raise ValueError("total_time must be positive.")

    if time_step <= 0.0:
        raise ValueError("time_step must be positive.")

    number_of_steps = (
        int(round(total_time / time_step)) + 1
    )

    time = np.linspace(
        0.0,
        total_time,
        number_of_steps,
        dtype=np.float64,
    )

    states = np.zeros(
        (number_of_steps, 4),
        dtype=np.float64,
    )

    strain = np.zeros(number_of_steps, dtype=np.float64)
    elastic_strain = np.zeros(number_of_steps, dtype=np.float64)
    stress = np.zeros(number_of_steps, dtype=np.float64)

    plastic_strain_rate = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    alpha_rate = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    beta = np.zeros(number_of_steps, dtype=np.float64)

    r_bar_rate = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    R_bar = np.zeros(number_of_steps, dtype=np.float64)
    R = np.zeros(number_of_steps, dtype=np.float64)

    damage_rate = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    energy_release_rate = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    effective_relative_stress = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    yield_function = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    plastic_multiplier = np.zeros(
        number_of_steps,
        dtype=np.float64,
    )

    for index in range(number_of_steps):
        current_time: float = float(time[index])
        current_state = states[index]

        strain_value: float = cast(
            float,
            strain_function(current_time),
        )
        strain[index] = strain_value

        elastic_strain[index] = float(
            strain[index] - current_state[0]
        )

        (
            stress[index],
            beta[index],
            R_bar[index],
            R[index],
            effective_relative_stress[index],
            yield_function[index],
            energy_release_rate[index],
        ) = evaluate_state(
            strain[index],
            current_state,
            material,
        )

        rates = state_rate(
            strain[index],
            current_state,
            material,
        )

        plastic_strain_rate[index] = float(rates[0])
        alpha_rate[index] = float(rates[1])
        r_bar_rate[index] = float(rates[2])
        damage_rate[index] = float(rates[3])

        plastic_multiplier[index] = float(
            material.k_viscoplastic
            * positive_part(
                yield_function[index]
            ) ** material.n
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


def plot_single_history(
    time: FloatArray,
    values: FloatArray,
    title: str,
    ylabel: str,
) -> None:
    """Create one independent history plot."""
    plt.figure(figsize=(8.0, 4.5))
    plt.plot(time, values)
    plt.xlabel(r"$t$ (s)")
    plt.ylabel(ylabel)
    plt.title(title)
    plt.grid(True)
    plt.tight_layout()


def plot_response(
    response: MaterialResponse,
) -> None:
    """Plot the principal response histories."""

    figures: List[PlotItem] = [
        (
            "Total strain",
            response.strain,
            r"$\varepsilon$",
        ),
        (
            "Stress",
            response.stress,
            r"$\sigma$ (MPa)",
        ),
        (
            "Plastic strain rate",
            response.plastic_strain_rate,
            r"$\dot{\varepsilon}^{p}$ (s$^{-1}$)",
        ),
        (
            "Backstress",
            response.beta,
            r"$\beta$ (MPa)",
        ),
        (
            "Transformed isotropic hardening force",
            response.R_bar,
            r"$\bar{R}$ (MPa)",
        ),
        (
            "Damage",
            response.damage,
            r"$D$",
        ),
        (
            "Damage rate",
            response.damage_rate,
            r"$\dot{D}$ (s$^{-1}$)",
        ),
        (
            "Energy release rate",
            response.energy_release_rate,
            r"$Y$ (MPa)",
        ),
        (
            "Yield function",
            response.yield_function,
            r"$f^{p}$ (MPa)",
        ),
    ]

    for title, values, ylabel in figures:
        plot_single_history(
            response.time,
            values,
            title,
            ylabel,
        )

    plt.show()


def print_summary(
    response: MaterialResponse,
    material: MaterialParameters,
) -> None:
    """Print important scalar results."""

    print("=" * 60)
    print("One-dimensional material-point simulation")
    print("=" * 60)
    print(f"Yield stress: {material.sigma_y:.3f} MPa")
    print(
        f"Norton coefficient k: "
        f"{material.k_viscoplastic:.6e}"
    )
    print(
        f"Damage threshold Y0: "
        f"{material.Y0:.6e} MPa"
    )
    print(
        f"Final plastic strain: "
        f"{response.plastic_strain[-1]:.6e}"
    )
    print(
        f"Final backstress beta: "
        f"{response.beta[-1]:.6f} MPa"
    )
    print(
        f"Final transformed hardening R_bar: "
        f"{response.R_bar[-1]:.6f} MPa"
    )
    print(
        f"Final physical hardening R: "
        f"{response.R[-1]:.6f} MPa"
    )
    print(
        f"Final damage: "
        f"{response.damage[-1]:.6f}"
    )
    print(
        f"Maximum damage: "
        f"{np.max(response.damage):.6f}"
    )
    print(
        f"Maximum tensile stress: "
        f"{np.max(response.stress):.6f} MPa"
    )
    print(
        f"Minimum compressive stress: "
        f"{np.min(response.stress):.6f} MPa"
    )
    print(
        f"Maximum energy release rate: "
        f"{np.max(response.energy_release_rate):.6e} MPa"
    )
    print("=" * 60)


def main() -> None:
    material = MaterialParameters(
        sigma_y=80.0,
    )

    response = run_material_point(
        material=material,
        total_time=200.0,
        time_step=0.1,
    )

    print_summary(
        response,
        material,
    )

    plot_response(response)


if __name__ == "__main__":
    main()

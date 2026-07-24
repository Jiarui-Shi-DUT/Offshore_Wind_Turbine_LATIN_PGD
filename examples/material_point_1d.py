from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np

from material.viscoplastic_damage_1d import (
    MaterialParameters,
    MaterialResponse,
    run_material_point,
)


def plot_response(response: MaterialResponse) -> None:
    """Plot the main material-point histories."""

    histories = [
        ("Stress", response.stress, r"$\sigma$ (MPa)"),
        (
            "Plastic strain rate",
            response.plastic_strain_rate,
            r"$\dot{\varepsilon}^{p}$ (s$^{-1}$)",
        ),
        ("Backstress", response.beta, r"$\beta$ (MPa)"),
        ("Isotropic hardening", response.R_bar, r"$\bar{R}$ (MPa)"),
        ("Damage", response.damage, r"$D$"),
        (
            "Energy release rate",
            response.energy_release_rate,
            r"$Y$ (MPa)",
        ),
    ]

    for title, values, ylabel in histories:
        plt.figure()
        plt.plot(response.time, values)
        plt.xlabel(r"$t$ (s)")
        plt.ylabel(ylabel)
        plt.title(title)
        plt.grid(True)
        plt.tight_layout()

    plt.show()


def print_summary(
    material: MaterialParameters,
    response: MaterialResponse,
) -> None:
    """Print the main scalar results."""

    print("=" * 60)
    print("One-dimensional material-point simulation")
    print("=" * 60)

    print(f"Yield stress: {material.sigma_y:.3f} MPa")
    print(
        "Norton coefficient k: "
        f"{material.k_viscoplastic:.6e}"
    )
    print(f"Damage threshold Y0: {material.Y0:.6e} MPa")
    print(
        "Final plastic strain: "
        f"{response.plastic_strain[-1]:.6e}"
    )
    print(
        "Final backstress beta: "
        f"{response.beta[-1]:.6f} MPa"
    )
    print(
        "Final transformed hardening R_bar: "
        f"{response.R_bar[-1]:.6f} MPa"
    )
    print(
        "Final physical hardening R: "
        f"{response.R[-1]:.6f} MPa"
    )
    print(f"Final damage: {response.damage[-1]:.6f}")
    print(f"Maximum damage: {np.max(response.damage):.6f}")
    print(
        "Maximum tensile stress: "
        f"{np.max(response.stress):.6f} MPa"
    )
    print(
        "Minimum compressive stress: "
        f"{np.min(response.stress):.6f} MPa"
    )
    print(
        "Maximum energy release rate: "
        f"{np.max(response.energy_release_rate):.6e} MPa"
    )
    print("=" * 60)


def main() -> None:
    material = MaterialParameters(sigma_y=80.0)

    response = run_material_point(
        material=material,
        total_time=200.0,
        time_step=0.1,
    )

    print_summary(material, response)
    plot_response(response)


if __name__ == "__main__":
    main()

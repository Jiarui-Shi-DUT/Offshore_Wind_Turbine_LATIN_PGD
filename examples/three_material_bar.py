# -*- coding: utf-8 -*-
"""
Reference three-material bar benchmark from Section 5.1 of the LATIN-PGD paper.

Full benchmark:
    length = 1000 mm
    area = 100 mm^2
    90 linear bar elements
    three equal material regions
    yield stresses = 80, 82.5, 85 MPa
    sinusoidal end displacement amplitude = 1.2 mm
    period = 10 s
    20 cycles
    time step = 0.1 s

Use:
    python -m examples.three_material_bar --quick
    python -m examples.three_material_bar
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import List, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np

from fem.bar_1d import (
    BarMesh1D,
    BarResponse,
    create_uniform_bar_mesh,
    solve_displacement_controlled_bar,
)
from material.viscoplastic_damage_1d import MaterialParameters


@dataclass(frozen=True)
class BenchmarkConfiguration:
    """Geometry, loading and discretisation of the bar benchmark."""

    length: float = 1000.0
    area: float = 100.0
    displacement_amplitude: float = 1.2
    period: float = 10.0
    cycles: int = 20
    time_step: float = 0.1
    n_elements: int = 90

    @property
    def total_time(self) -> float:
        return float(self.cycles * self.period)


def prescribed_displacement(
    time: np.ndarray,
    amplitude: float,
    period: float,
) -> np.ndarray:
    """Return the sinusoidal displacement prescribed at x = L."""
    return amplitude * np.sin(2.0 * np.pi * time / period)


def create_three_material_distribution(
    n_elements: int,
) -> List[MaterialParameters]:
    """Assign one material to each third of the bar."""
    if n_elements < 3 or n_elements % 3 != 0:
        raise ValueError(
            "n_elements must be a positive multiple of three."
        )

    elements_per_material = n_elements // 3

    material_1 = MaterialParameters(sigma_y=80.0)
    material_2 = MaterialParameters(sigma_y=82.5)
    material_3 = MaterialParameters(sigma_y=85.0)

    return (
        [material_1] * elements_per_material
        + [material_2] * elements_per_material
        + [material_3] * elements_per_material
    )


def create_time_grid(
    total_time: float,
    time_step: float,
) -> np.ndarray:
    """Create a time grid containing both t = 0 and t = total_time."""
    number_of_steps = int(round(total_time / time_step))

    if not np.isclose(
        number_of_steps * time_step,
        total_time,
        rtol=0.0,
        atol=1.0e-12,
    ):
        raise ValueError(
            "total_time must be an integer multiple of time_step."
        )

    return np.linspace(
        0.0,
        total_time,
        number_of_steps + 1,
        dtype=np.float64,
    )


def run_benchmark(
    configuration: BenchmarkConfiguration,
) -> Tuple[BarMesh1D, BarResponse]:
    """Build and solve the displacement-controlled three-material bar."""
    mesh = create_uniform_bar_mesh(
        length=configuration.length,
        n_elements=configuration.n_elements,
    )
    materials = create_three_material_distribution(
        configuration.n_elements
    )
    time = create_time_grid(
        total_time=configuration.total_time,
        time_step=configuration.time_step,
    )
    right_displacement = prescribed_displacement(
        time=time,
        amplitude=configuration.displacement_amplitude,
        period=configuration.period,
    )

    response = solve_displacement_controlled_bar(
        mesh=mesh,
        area=configuration.area,
        materials=materials,
        time=time,
        right_displacement=right_displacement,
    )
    return mesh, response


def material_region_slices(
    n_elements: int,
) -> Sequence[slice]:
    """Return the element slices occupied by the three materials."""
    elements_per_material = n_elements // 3
    return (
        slice(0, elements_per_material),
        slice(elements_per_material, 2 * elements_per_material),
        slice(2 * elements_per_material, 3 * elements_per_material),
    )


def region_average(
    values: np.ndarray,
    region: slice,
) -> np.ndarray:
    """Average an element history over one material region."""
    return np.mean(values[:, region], axis=1)


def print_summary(
    configuration: BenchmarkConfiguration,
    response: BarResponse,
) -> None:
    """Print the main numerical checks and final damage values."""
    regions = material_region_slices(configuration.n_elements)
    damage = response.state[:, :, 3]

    final_damage = [
        float(np.mean(damage[-1, region]))
        for region in regions
    ]

    maximum_stress_spread = float(
        np.max(np.ptp(response.stress, axis=1))
    )
    maximum_force_imbalance = float(
        np.max(
            np.abs(
                response.reaction_left
                + response.reaction_right
            )
        )
    )

    print("=" * 68)
    print("Three-material bar benchmark")
    print("=" * 68)
    print(f"Elements: {configuration.n_elements}")
    print(f"Cycles: {configuration.cycles}")
    print(f"Time step: {configuration.time_step:.6g} s")
    print(f"Time increments: {response.time.size - 1}")
    print(
        "Final regional damage: "
        f"mat. 1 = {final_damage[0]:.6f}, "
        f"mat. 2 = {final_damage[1]:.6f}, "
        f"mat. 3 = {final_damage[2]:.6f}"
    )
    print(
        "Maximum spatial stress spread: "
        f"{maximum_stress_spread:.6e} MPa"
    )
    print(
        "Maximum end-force imbalance: "
        f"{maximum_force_imbalance:.6e} N"
    )
    print(
        "Maximum Newton iterations: "
        f"{int(np.max(response.newton_iterations))}"
    )
    print("=" * 68)


def plot_results(
    mesh: BarMesh1D,
    configuration: BenchmarkConfiguration,
    response: BarResponse,
) -> None:
    """Plot the prescribed loading and principal benchmark responses."""
    regions = material_region_slices(configuration.n_elements)
    damage = response.state[:, :, 3]
    element_centres = 0.5 * (
        mesh.coordinates[:-1] + mesh.coordinates[1:]
    )

    plt.figure()
    plt.plot(response.time, response.displacement[:, -1])
    plt.xlabel(r"$t$ (s)")
    plt.ylabel(r"$u_d$ (mm)")
    plt.title("Prescribed end displacement")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(response.time, np.mean(response.stress, axis=1))
    plt.xlabel(r"$t$ (s)")
    plt.ylabel(r"$\sigma$ (MPa)")
    plt.title("Axial stress")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    for material_number, region in enumerate(regions, start=1):
        plt.plot(
            response.time,
            region_average(damage, region),
            label=f"mat. {material_number}",
        )
    plt.xlabel(r"$t$ (s)")
    plt.ylabel(r"$D$")
    plt.title("Damage evolution")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(element_centres, damage[-1])
    plt.xlabel(r"$x$ (mm)")
    plt.ylabel(r"$D$")
    plt.title("Final damage distribution")
    plt.grid(True)
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description="Run the three-material cyclic bar benchmark."
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help=(
            "Run a one-cycle, nine-element smoke test instead "
            "of the full paper benchmark."
        ),
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Run the analysis without opening figures.",
    )
    return parser.parse_args()


def main() -> None:
    arguments = parse_arguments()

    if arguments.quick:
        configuration = BenchmarkConfiguration(
            cycles=1,
            n_elements=9,
        )
    else:
        configuration = BenchmarkConfiguration()

    mesh, response = run_benchmark(configuration)
    print_summary(configuration, response)

    if not arguments.no_plot:
        plot_results(mesh, configuration, response)


if __name__ == "__main__":
    main()

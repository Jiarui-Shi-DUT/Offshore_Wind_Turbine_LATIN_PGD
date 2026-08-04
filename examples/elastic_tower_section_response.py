# -*- coding: utf-8 -*-
"""
Elastic section and fiber response of the tapered offshore-wind-turbine tower.

This example extends the global elastic tower solution by recovering, at
every beam Gauss point:

    - axial strain;
    - curvature;
    - axial force;
    - bending moment;
    - minimum and maximum fiber strain;
    - minimum and maximum fiber stress.

Use:
    python -m examples.elastic_tower_section_response --no-plot
    python -m examples.elastic_tower_section_response
"""

from __future__ import annotations

import argparse
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from examples.elastic_tapered_tower import (
    TowerConfiguration,
    run_analysis,
)
from fem.tower_response_2d import (
    ElasticTowerResponse,
    recover_elastic_tower_response,
)


def recover_response(
    configuration: TowerConfiguration,
) -> ElasticTowerResponse:
    """Solve the tower and recover all Gauss-point and fiber responses."""
    _, tower_geometry, assembly, solution = run_analysis(
        configuration
    )

    return recover_elastic_tower_response(
        assembly=assembly,
        solution=solution,
        tower_geometry=tower_geometry,
        elastic_modulus=configuration.elastic_modulus,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )


def critical_response_indices(
    response: ElasticTowerResponse,
) -> Tuple[int, int]:
    """Return the element and Gauss indices of maximum absolute fiber stress."""
    maximum_absolute_stress = np.maximum(
        np.abs(response.minimum_fiber_stresses),
        np.abs(response.maximum_fiber_stresses),
    )
    flat_index = int(np.argmax(maximum_absolute_stress))
    element_index, gauss_index = np.unravel_index(
        flat_index,
        maximum_absolute_stress.shape,
    )
    return int(element_index), int(gauss_index)


def print_response_summary(
    configuration: TowerConfiguration,
    response: ElasticTowerResponse,
) -> None:
    """Print the critical elastic section and fiber response."""
    element_index, gauss_index = critical_response_indices(
        response
    )

    height = float(
        response.gauss_heights[element_index, gauss_index]
    )
    axial_strain = float(
        response.axial_strains[element_index, gauss_index]
    )
    curvature = float(
        response.curvatures[element_index, gauss_index]
    )
    axial_force = float(
        response.axial_forces[element_index, gauss_index]
    )
    bending_moment = float(
        response.bending_moments[element_index, gauss_index]
    )
    minimum_strain = float(
        response.minimum_fiber_strains[
            element_index,
            gauss_index,
        ]
    )
    maximum_strain = float(
        response.maximum_fiber_strains[
            element_index,
            gauss_index,
        ]
    )
    minimum_stress = float(
        response.minimum_fiber_stresses[
            element_index,
            gauss_index,
        ]
    )
    maximum_stress = float(
        response.maximum_fiber_stresses[
            element_index,
            gauss_index,
        ]
    )

    expected_bending_moment = (
        -configuration.horizontal_force
        * (configuration.height - height)
    )
    moment_error = bending_moment - expected_bending_moment
    moment_relative_error = (
        abs(moment_error)
        / max(1.0, abs(expected_bending_moment))
    )

    maximum_absolute_stress = max(
        abs(minimum_stress),
        abs(maximum_stress),
    )

    print("=" * 76)
    print("Elastic tapered-tower section and fiber response")
    print("=" * 76)
    print(f"Tower height: {configuration.height:.6g} m")
    print(f"Beam elements: {configuration.n_elements}")
    print(f"Gauss points per element: {configuration.n_gauss}")
    print(
        "Fibers per section: "
        f"{configuration.n_circumferential} x "
        f"{configuration.n_radial} = "
        f"{response.n_fibers}"
    )
    print(
        "Top horizontal force: "
        f"{configuration.horizontal_force:.6e} N"
    )
    print("-" * 76)
    print("Critical discrete integration point")
    print(f"Element index: {element_index}")
    print(f"Gauss-point index: {gauss_index}")
    print(f"Height above tower base: {height:.6e} m")
    print("-" * 76)
    print(f"Section axial strain: {axial_strain:.6e}")
    print(f"Section curvature: {curvature:.6e} 1/m")
    print(f"Section axial force: {axial_force:.6e} N")
    print(
        "Recovered bending moment: "
        f"{bending_moment:.6e} N m"
    )
    print(
        "Analytical bending moment: "
        f"{expected_bending_moment:.6e} N m"
    )
    print(
        "Bending-moment relative error: "
        f"{moment_relative_error:.6e}"
    )
    print("-" * 76)
    print(f"Minimum fiber strain: {minimum_strain:.6e}")
    print(f"Maximum fiber strain: {maximum_strain:.6e}")
    print(
        "Minimum fiber stress: "
        f"{minimum_stress / 1.0e6:.6e} MPa"
    )
    print(
        "Maximum fiber stress: "
        f"{maximum_stress / 1.0e6:.6e} MPa"
    )
    print(
        "Maximum absolute fiber stress: "
        f"{maximum_absolute_stress / 1.0e6:.6e} MPa"
    )
    print("=" * 76)


def plot_response(
    configuration: TowerConfiguration,
    response: ElasticTowerResponse,
) -> None:
    """Plot the recovered elastic response along the tower height."""
    heights = response.flattened_gauss_heights
    bending_moments = (
        response.bending_moments.reshape(-1) / 1.0e6
    )
    curvatures = response.curvatures.reshape(-1)
    minimum_stresses = (
        response.minimum_fiber_stresses.reshape(-1)
        / 1.0e6
    )
    maximum_stresses = (
        response.maximum_fiber_stresses.reshape(-1)
        / 1.0e6
    )
    analytical_moments = (
        -configuration.horizontal_force
        * (configuration.height - heights)
        / 1.0e6
    )

    plt.figure()
    plt.plot(bending_moments, heights, label="Recovered")
    plt.plot(
        analytical_moments,
        heights,
        linestyle="--",
        label="Analytical",
    )
    plt.xlabel("Bending moment (MN m)")
    plt.ylabel("Height above tower base (m)")
    plt.title("Tower bending-moment distribution")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(curvatures, heights)
    plt.xlabel("Curvature (1/m)")
    plt.ylabel("Height above tower base (m)")
    plt.title("Tower curvature distribution")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        minimum_stresses,
        heights,
        label="Minimum fiber stress",
    )
    plt.plot(
        maximum_stresses,
        heights,
        label="Maximum fiber stress",
    )
    plt.xlabel("Fiber stress (MPa)")
    plt.ylabel("Height above tower base (m)")
    plt.title("Extreme fiber-stress distribution")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Recover elastic section and fiber response of the "
            "tapered offshore-wind-turbine tower."
        )
    )
    parser.add_argument(
        "--elements",
        type=int,
        default=40,
        help="Number of uniform Euler-Bernoulli beam elements.",
    )
    parser.add_argument(
        "--force",
        type=float,
        default=1.0e6,
        help="Tower-top horizontal force in newtons.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Run the analysis without opening figures.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the elastic section-response example."""
    arguments = parse_arguments()
    configuration = TowerConfiguration(
        horizontal_force=arguments.force,
        n_elements=arguments.elements,
    )
    response = recover_response(configuration)

    print_response_summary(
        configuration=configuration,
        response=response,
    )

    if not arguments.no_plot:
        plot_response(
            configuration=configuration,
            response=response,
        )


if __name__ == "__main__":
    main()

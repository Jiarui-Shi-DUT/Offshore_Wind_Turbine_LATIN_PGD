# -*- coding: utf-8 -*-
"""
Linear-elastic tapered offshore-wind-turbine tower example.

The model uses:
    - NREL 5 MW nominal tower geometry;
    - 2D Euler-Bernoulli beam elements;
    - annular fiber sections;
    - four Gauss points per beam element;
    - fixed tower base;
    - horizontal concentrated force at the tower top.

Use:
    python -m examples.elastic_tapered_tower --no-plot
    python -m examples.elastic_tapered_tower
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from fem.beam_column_2d import (
    BeamMesh2D,
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.tower_system_2d import (
    ElasticTowerAssembly,
    LinearStaticSolution,
    solve_elastic_tower_top_force,
)


@dataclass(frozen=True)
class TowerConfiguration:
    """Geometry, material, loading, and discretisation parameters."""

    height: float = 87.6
    base_outer_diameter: float = 6.0
    top_outer_diameter: float = 3.87
    base_thickness: float = 0.027
    top_thickness: float = 0.019
    elastic_modulus: float = 210.0e9
    horizontal_force: float = 1.0e6
    n_elements: int = 40
    n_gauss: int = 4
    n_circumferential: int = 32
    n_radial: int = 2


def create_tower_geometry(
    configuration: TowerConfiguration,
) -> LinearTaperedTowerGeometry:
    """Create the linearly tapered annular tower geometry."""
    return LinearTaperedTowerGeometry(
        height=configuration.height,
        base_outer_diameter=configuration.base_outer_diameter,
        top_outer_diameter=configuration.top_outer_diameter,
        base_thickness=configuration.base_thickness,
        top_thickness=configuration.top_thickness,
    )


def run_analysis(
    configuration: TowerConfiguration,
) -> Tuple[
    BeamMesh2D,
    LinearTaperedTowerGeometry,
    ElasticTowerAssembly,
    LinearStaticSolution,
]:
    """Build and solve the elastic tower under a top horizontal force."""
    mesh = create_uniform_vertical_tower_mesh(
        height=configuration.height,
        n_elements=configuration.n_elements,
    )
    tower_geometry = create_tower_geometry(configuration)

    assembly, solution = solve_elastic_tower_top_force(
        mesh=mesh,
        tower_geometry=tower_geometry,
        elastic_modulus=configuration.elastic_modulus,
        horizontal_force=configuration.horizontal_force,
        n_gauss=configuration.n_gauss,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )

    return mesh, tower_geometry, assembly, solution


def print_summary(
    configuration: TowerConfiguration,
    mesh: BeamMesh2D,
    solution: LinearStaticSolution,
) -> None:
    """Print the principal displacement and equilibrium results."""
    top_horizontal_displacement = float(
        solution.displacements[-3]
    )
    top_vertical_displacement = float(
        solution.displacements[-2]
    )
    top_rotation = float(solution.displacements[-1])

    base_horizontal_reaction = float(solution.reactions[0])
    base_vertical_reaction = float(solution.reactions[1])
    base_moment_reaction = float(solution.reactions[2])

    expected_base_horizontal_reaction = (
        -configuration.horizontal_force
    )
    expected_base_moment_reaction = (
        configuration.horizontal_force * configuration.height
    )

    free_dof_residual = float(
        np.max(
            np.abs(
                solution.reactions[solution.free_dofs]
            )
        )
    )

    print("=" * 72)
    print("Elastic tapered offshore-wind-turbine tower")
    print("=" * 72)
    print(f"Tower height: {configuration.height:.6g} m")
    print(f"Beam elements: {configuration.n_elements}")
    print(f"Nodes: {mesh.n_nodes}")
    print(f"Total DOFs: {mesh.n_dof}")
    print(f"Gauss points per element: {configuration.n_gauss}")
    print(
        "Fibers per section: "
        f"{configuration.n_circumferential} x "
        f"{configuration.n_radial} = "
        f"{configuration.n_circumferential * configuration.n_radial}"
    )
    print(
        "Top horizontal force: "
        f"{configuration.horizontal_force:.6e} N"
    )
    print("-" * 72)
    print(
        "Top horizontal displacement: "
        f"{top_horizontal_displacement:.6e} m"
    )
    print(
        "Top vertical displacement: "
        f"{top_vertical_displacement:.6e} m"
    )
    print(f"Top rotation: {top_rotation:.6e} rad")
    print("-" * 72)
    print(
        "Base horizontal reaction: "
        f"{base_horizontal_reaction:.6e} N"
    )
    print(
        "Expected horizontal reaction: "
        f"{expected_base_horizontal_reaction:.6e} N"
    )
    print(
        "Base vertical reaction: "
        f"{base_vertical_reaction:.6e} N"
    )
    print(
        "Base moment reaction: "
        f"{base_moment_reaction:.6e} N m"
    )
    print(
        "Expected base moment: "
        f"{expected_base_moment_reaction:.6e} N m"
    )
    print(
        "Maximum free-DOF equilibrium residual: "
        f"{free_dof_residual:.6e} N"
    )
    print("=" * 72)


def plot_results(
    mesh: BeamMesh2D,
    solution: LinearStaticSolution,
) -> None:
    """Plot horizontal displacement and rotation along the tower."""
    height = mesh.tower_axis_coordinates
    horizontal_displacement = solution.displacements[0::3]
    rotation = solution.displacements[2::3]

    plt.figure()
    plt.plot(horizontal_displacement, height)
    plt.xlabel(r"$u_x$ (m)")
    plt.ylabel(r"$z$ (m)")
    plt.title("Tower horizontal displacement")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(rotation, height)
    plt.xlabel(r"$\theta_z$ (rad)")
    plt.ylabel(r"$z$ (m)")
    plt.title("Tower section rotation")
    plt.grid(True)
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Run the linear-elastic tapered offshore-wind-turbine "
            "tower example."
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
    arguments = parse_arguments()

    configuration = TowerConfiguration(
        horizontal_force=arguments.force,
        n_elements=arguments.elements,
    )

    mesh, _, _, solution = run_analysis(configuration)
    print_summary(
        configuration=configuration,
        mesh=mesh,
        solution=solution,
    )

    if not arguments.no_plot:
        plot_results(mesh=mesh, solution=solution)


if __name__ == "__main__":
    main()

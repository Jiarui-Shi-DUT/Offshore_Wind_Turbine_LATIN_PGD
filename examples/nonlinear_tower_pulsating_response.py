# -*- coding: utf-8 -*-
"""
One positive-mean nonlinear pulsating cycle of the NREL 5 MW tower.

The example uses the default viscoplastic-damage material model and
records both global tower response and every fiber state. The resulting
data provide:

    - tower-top force-displacement hysteresis;
    - maximum plastic-strain history;
    - maximum damage history;
    - stress-strain response of one fixed critical fiber.

The default discretisation is intentionally moderate for verification.
A finer model can be selected with command-line options.

Use:
    python -m examples.nonlinear_tower_pulsating_response --no-plot
    python -m examples.nonlinear_tower_pulsating_response
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from examples.elastic_tapered_tower import (
    TowerConfiguration,
    create_tower_geometry,
)
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import (
    PulsatingTopForceHistory,
    create_pulsating_top_force_history,
)
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
    ViscoplasticTowerResponse,
    solve_nonlinear_tower_load_step,
)
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
FiberLocation = Tuple[int, int, int]


@dataclass(frozen=True)
class NonlinearPulsatingResponse:
    """Global and fiber histories from one sequential tower analysis."""

    loading: PulsatingTopForceHistory
    analysis_times: FloatArray
    top_displacements: FloatArray
    top_rotations: FloatArray
    base_horizontal_reactions: FloatArray
    base_moment_reactions: FloatArray
    iterations: IntArray
    residual_norms: FloatArray
    fiber_strains: FloatArray
    fiber_stresses: FloatArray
    fiber_plastic_strains: FloatArray
    fiber_damages: FloatArray
    critical_location: FiberLocation
    critical_height: float
    critical_y_coordinate: float

    def __post_init__(self) -> None:
        if not isinstance(
            self.loading,
            PulsatingTopForceHistory,
        ):
            raise TypeError(
                "loading must be a PulsatingTopForceHistory."
            )

        n_points = self.loading.n_time_points
        analysis_times = np.asarray(
            self.analysis_times,
            dtype=np.float64,
        )
        top_displacements = np.asarray(
            self.top_displacements,
            dtype=np.float64,
        )
        top_rotations = np.asarray(
            self.top_rotations,
            dtype=np.float64,
        )
        base_horizontal_reactions = np.asarray(
            self.base_horizontal_reactions,
            dtype=np.float64,
        )
        base_moment_reactions = np.asarray(
            self.base_moment_reactions,
            dtype=np.float64,
        )
        iterations = np.asarray(
            self.iterations,
            dtype=np.int64,
        )
        residual_norms = np.asarray(
            self.residual_norms,
            dtype=np.float64,
        )

        one_dimensional_arrays = (
            analysis_times,
            top_displacements,
            top_rotations,
            base_horizontal_reactions,
            base_moment_reactions,
            iterations,
            residual_norms,
        )
        for array in one_dimensional_arrays:
            if array.shape != (n_points,):
                raise ValueError(
                    "Every global history must match loading.times."
                )

        fiber_strains = np.asarray(
            self.fiber_strains,
            dtype=np.float64,
        )
        fiber_stresses = np.asarray(
            self.fiber_stresses,
            dtype=np.float64,
        )
        fiber_plastic_strains = np.asarray(
            self.fiber_plastic_strains,
            dtype=np.float64,
        )
        fiber_damages = np.asarray(
            self.fiber_damages,
            dtype=np.float64,
        )

        if fiber_strains.ndim != 4:
            raise ValueError(
                "Fiber histories must have shape "
                "(n_points, n_elements, n_gauss, n_fibers)."
            )
        if fiber_strains.shape[0] != n_points:
            raise ValueError(
                "Fiber histories must match loading.times."
            )

        fiber_arrays = (
            fiber_strains,
            fiber_stresses,
            fiber_plastic_strains,
            fiber_damages,
        )
        for array in fiber_arrays:
            if array.shape != fiber_strains.shape:
                raise ValueError(
                    "All fiber histories must have equal shape."
                )
            if np.any(~np.isfinite(array)):
                raise ValueError(
                    "All fiber histories must be finite."
                )

        if np.any(~np.isfinite(analysis_times)):
            raise ValueError("analysis_times must be finite.")
        if np.any(np.diff(analysis_times) <= 0.0):
            raise ValueError(
                "analysis_times must be strictly increasing."
            )
        if np.any(~np.isfinite(top_displacements)):
            raise ValueError(
                "top_displacements must be finite."
            )
        if np.any(~np.isfinite(top_rotations)):
            raise ValueError("top_rotations must be finite.")
        if np.any(~np.isfinite(base_horizontal_reactions)):
            raise ValueError(
                "base_horizontal_reactions must be finite."
            )
        if np.any(~np.isfinite(base_moment_reactions)):
            raise ValueError(
                "base_moment_reactions must be finite."
            )
        if np.any(iterations < 1):
            raise ValueError(
                "Every Newton iteration count must be positive."
            )
        if np.any(~np.isfinite(residual_norms)):
            raise ValueError("residual_norms must be finite.")
        if np.any(residual_norms < 0.0):
            raise ValueError(
                "residual_norms must be non-negative."
            )
        if np.any(fiber_damages < 0.0):
            raise ValueError(
                "Fiber damage must be non-negative."
            )

        location = tuple(int(value) for value in self.critical_location)
        if len(location) != 3:
            raise ValueError(
                "critical_location must contain three indices."
            )
        element_index, gauss_index, fiber_index = location
        if (
            element_index < 0
            or element_index >= fiber_strains.shape[1]
            or gauss_index < 0
            or gauss_index >= fiber_strains.shape[2]
            or fiber_index < 0
            or fiber_index >= fiber_strains.shape[3]
        ):
            raise ValueError(
                "critical_location contains an invalid index."
            )

        critical_height = float(self.critical_height)
        critical_y_coordinate = float(
            self.critical_y_coordinate
        )
        if not np.isfinite(critical_height):
            raise ValueError("critical_height must be finite.")
        if not np.isfinite(critical_y_coordinate):
            raise ValueError(
                "critical_y_coordinate must be finite."
            )

        object.__setattr__(
            self,
            "analysis_times",
            analysis_times.copy(),
        )
        object.__setattr__(
            self,
            "top_displacements",
            top_displacements.copy(),
        )
        object.__setattr__(
            self,
            "top_rotations",
            top_rotations.copy(),
        )
        object.__setattr__(
            self,
            "base_horizontal_reactions",
            base_horizontal_reactions.copy(),
        )
        object.__setattr__(
            self,
            "base_moment_reactions",
            base_moment_reactions.copy(),
        )
        object.__setattr__(
            self,
            "iterations",
            iterations.copy(),
        )
        object.__setattr__(
            self,
            "residual_norms",
            residual_norms.copy(),
        )
        object.__setattr__(
            self,
            "fiber_strains",
            fiber_strains.copy(),
        )
        object.__setattr__(
            self,
            "fiber_stresses",
            fiber_stresses.copy(),
        )
        object.__setattr__(
            self,
            "fiber_plastic_strains",
            fiber_plastic_strains.copy(),
        )
        object.__setattr__(
            self,
            "fiber_damages",
            fiber_damages.copy(),
        )
        object.__setattr__(
            self,
            "critical_location",
            location,
        )
        object.__setattr__(
            self,
            "critical_height",
            critical_height,
        )
        object.__setattr__(
            self,
            "critical_y_coordinate",
            critical_y_coordinate,
        )

    @property
    def maximum_absolute_stresses(self) -> FloatArray:
        """Return maximum absolute fiber stress at every time point."""
        return np.max(
            np.abs(self.fiber_stresses),
            axis=(1, 2, 3),
        )

    @property
    def maximum_absolute_plastic_strains(self) -> FloatArray:
        """Return maximum absolute plastic strain at every time point."""
        return np.max(
            np.abs(self.fiber_plastic_strains),
            axis=(1, 2, 3),
        )

    @property
    def maximum_damages(self) -> FloatArray:
        """Return maximum fiber damage at every time point."""
        return np.max(
            self.fiber_damages,
            axis=(1, 2, 3),
        )

    @property
    def critical_fiber_strains(self) -> FloatArray:
        """Return total strain history at the fixed critical fiber."""
        element_index, gauss_index, fiber_index = (
            self.critical_location
        )
        return self.fiber_strains[
            :,
            element_index,
            gauss_index,
            fiber_index,
        ].copy()

    @property
    def critical_fiber_stresses(self) -> FloatArray:
        """Return stress history at the fixed critical fiber."""
        element_index, gauss_index, fiber_index = (
            self.critical_location
        )
        return self.fiber_stresses[
            :,
            element_index,
            gauss_index,
            fiber_index,
        ].copy()

    @property
    def critical_fiber_plastic_strains(self) -> FloatArray:
        """Return plastic-strain history at the fixed critical fiber."""
        element_index, gauss_index, fiber_index = (
            self.critical_location
        )
        return self.fiber_plastic_strains[
            :,
            element_index,
            gauss_index,
            fiber_index,
        ].copy()

    @property
    def critical_fiber_damages(self) -> FloatArray:
        """Return damage history at the fixed critical fiber."""
        element_index, gauss_index, fiber_index = (
            self.critical_location
        )
        return self.fiber_damages[
            :,
            element_index,
            gauss_index,
            fiber_index,
        ].copy()


def snapshot_fiber_fields(
    response: ViscoplasticTowerResponse,
) -> Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    """Stack every element, Gauss-point, and fiber response."""
    strains = []
    stresses = []
    plastic_strains = []
    damages = []

    for element_response in response.element_responses:
        element_strains = []
        element_stresses = []
        element_plastic_strains = []
        element_damages = []

        for section_response in element_response.section_responses:
            element_strains.append(
                section_response.fiber_strains
            )
            element_stresses.append(
                section_response.fiber_stresses
            )
            element_plastic_strains.append(
                section_response.plastic_strains
            )
            element_damages.append(
                section_response.damages
            )

        strains.append(np.stack(element_strains, axis=0))
        stresses.append(np.stack(element_stresses, axis=0))
        plastic_strains.append(
            np.stack(element_plastic_strains, axis=0)
        )
        damages.append(np.stack(element_damages, axis=0))

    return (
        np.stack(strains, axis=0),
        np.stack(stresses, axis=0),
        np.stack(plastic_strains, axis=0),
        np.stack(damages, axis=0),
    )


def locate_final_critical_fiber(
    fiber_damages: FloatArray,
    fiber_plastic_strains: FloatArray,
) -> FiberLocation:
    """
    Select one fixed critical fiber for path-dependent plots.

    Final damage is the primary criterion. If all final damage values
    are zero, final absolute plastic strain is used instead.
    """
    final_damage = fiber_damages[-1]
    if float(np.max(final_damage)) > 0.0:
        flat_index = int(np.argmax(final_damage))
        location = np.unravel_index(
            flat_index,
            final_damage.shape,
        )
    else:
        final_plastic = np.abs(
            fiber_plastic_strains[-1]
        )
        flat_index = int(np.argmax(final_plastic))
        location = np.unravel_index(
            flat_index,
            final_plastic.shape,
        )

    return (
        int(location[0]),
        int(location[1]),
        int(location[2]),
    )


def run_nonlinear_pulsating_analysis(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    loading: PulsatingTopForceHistory,
    max_iterations: int = 40,
) -> NonlinearPulsatingResponse:
    """Build the tower and advance one complete load history."""
    mesh = create_uniform_vertical_tower_mesh(
        height=configuration.height,
        n_elements=configuration.n_elements,
    )
    geometry = create_tower_geometry(configuration)
    system = ViscoplasticDamageTowerSystem2D(
        mesh=mesh,
        tower_geometry=geometry,
        material=material,
        n_gauss=configuration.n_gauss,
        n_circumferential=(
            configuration.n_circumferential
        ),
        n_radial=configuration.n_radial,
    )

    preload_duration = loading.time_increment
    analysis_times = (
        preload_duration + loading.times
    )
    n_points = loading.n_time_points

    top_displacements = np.empty(
        n_points,
        dtype=np.float64,
    )
    top_rotations = np.empty(
        n_points,
        dtype=np.float64,
    )
    base_horizontal_reactions = np.empty(
        n_points,
        dtype=np.float64,
    )
    base_moment_reactions = np.empty(
        n_points,
        dtype=np.float64,
    )
    iterations = np.empty(
        n_points,
        dtype=np.int64,
    )
    residual_norms = np.empty(
        n_points,
        dtype=np.float64,
    )

    strain_snapshots = []
    stress_snapshots = []
    plastic_snapshots = []
    damage_snapshots = []

    for time_index in range(n_points):
        load_vector = top_horizontal_load_vector(
            mesh=mesh,
            horizontal_force=float(
                loading.forces[time_index]
            ),
        )
        solution = solve_nonlinear_tower_load_step(
            system=system,
            time=float(analysis_times[time_index]),
            load_vector=load_vector,
            max_iterations=max_iterations,
            relative_residual_tolerance=1.0e-8,
            absolute_residual_tolerance=1.0e-4,
        )

        top_displacements[time_index] = (
            solution.displacements[-3]
        )
        top_rotations[time_index] = (
            solution.displacements[-1]
        )
        base_horizontal_reactions[time_index] = (
            solution.reactions[0]
        )
        base_moment_reactions[time_index] = (
            solution.reactions[2]
        )
        iterations[time_index] = solution.iterations
        residual_norms[time_index] = (
            solution.residual_norm
        )

        (
            fiber_strains,
            fiber_stresses,
            fiber_plastic_strains,
            fiber_damages,
        ) = snapshot_fiber_fields(solution.response)

        strain_snapshots.append(fiber_strains)
        stress_snapshots.append(fiber_stresses)
        plastic_snapshots.append(fiber_plastic_strains)
        damage_snapshots.append(fiber_damages)

    fiber_strain_history = np.stack(
        strain_snapshots,
        axis=0,
    )
    fiber_stress_history = np.stack(
        stress_snapshots,
        axis=0,
    )
    fiber_plastic_history = np.stack(
        plastic_snapshots,
        axis=0,
    )
    fiber_damage_history = np.stack(
        damage_snapshots,
        axis=0,
    )

    critical_location = locate_final_critical_fiber(
        fiber_damages=fiber_damage_history,
        fiber_plastic_strains=fiber_plastic_history,
    )
    element_index, gauss_index, fiber_index = (
        critical_location
    )
    critical_element = system.elements[element_index]
    critical_section = critical_element.sections[gauss_index]

    return NonlinearPulsatingResponse(
        loading=loading,
        analysis_times=analysis_times,
        top_displacements=top_displacements,
        top_rotations=top_rotations,
        base_horizontal_reactions=(
            base_horizontal_reactions
        ),
        base_moment_reactions=base_moment_reactions,
        iterations=iterations,
        residual_norms=residual_norms,
        fiber_strains=fiber_strain_history,
        fiber_stresses=fiber_stress_history,
        fiber_plastic_strains=fiber_plastic_history,
        fiber_damages=fiber_damage_history,
        critical_location=critical_location,
        critical_height=float(
            critical_element.gauss_heights[gauss_index]
        ),
        critical_y_coordinate=float(
            critical_section.section.y_coordinates[
                fiber_index
            ]
        ),
    )


def print_summary(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    response: NonlinearPulsatingResponse,
) -> None:
    """Print loading, convergence, and nonlinear-response results."""
    loading = response.loading
    quarter = loading.increments_per_cycle // 4
    indices = (
        0,
        quarter,
        2 * quarter,
        3 * quarter,
        4 * quarter,
    )

    print("=" * 108)
    print("Nonlinear NREL 5 MW tower under one pulsating cycle")
    print("=" * 108)
    print(
        f"F_max={loading.maximum_force / 1.0e6:.6g} MN, "
        f"F_min={loading.minimum_force / 1.0e6:.6g} MN, "
        f"R_F={loading.force_ratio:.6g}"
    )
    print(
        f"Period={loading.period:.6g}, "
        f"increments/cycle={loading.increments_per_cycle}"
    )
    print(
        f"Material: E={material.E:.6g} MPa, "
        f"sigma_y={material.sigma_y:.6g} MPa, "
        f"k_damage={material.k_damage:.6g}"
    )
    print(
        "Discretisation: "
        f"{configuration.n_elements} elements, "
        f"{configuration.n_gauss} Gauss points/element, "
        f"{configuration.n_circumferential} x "
        f"{configuration.n_radial} fibers/section"
    )
    print("-" * 108)
    print(
        " t/T      force(MN)    top ux(m)    "
        "max|stress|(MPa)   max|eps_p|      max D      Newton"
    )

    maximum_stresses = (
        response.maximum_absolute_stresses
    )
    maximum_plastic_strains = (
        response.maximum_absolute_plastic_strains
    )
    maximum_damages = response.maximum_damages

    for index in indices:
        print(
            f"{loading.times[index] / loading.period:6.2f}   "
            f"{loading.forces[index] / 1.0e6:10.4f}   "
            f"{response.top_displacements[index]:10.4e}   "
            f"{maximum_stresses[index]:17.6e}   "
            f"{maximum_plastic_strains[index]:11.4e}   "
            f"{maximum_damages[index]:9.4e}   "
            f"{response.iterations[index]:6d}"
        )

    element_index, gauss_index, fiber_index = (
        response.critical_location
    )
    print("-" * 108)
    print(
        "Fixed critical fiber: "
        f"element {element_index + 1}, "
        f"Gauss point {gauss_index + 1}, "
        f"fiber {fiber_index + 1}"
    )
    print(
        f"Critical height={response.critical_height:.6e} m, "
        f"y={response.critical_y_coordinate:.6e} m"
    )
    print(
        "Final maximum plastic strain: "
        f"{maximum_plastic_strains[-1]:.6e}"
    )
    print(
        "Final maximum damage: "
        f"{maximum_damages[-1]:.6e}"
    )
    print(
        "Maximum Newton iterations in one step: "
        f"{int(np.max(response.iterations))}"
    )
    print(
        "Maximum free-DOF residual norm: "
        f"{float(np.max(response.residual_norms)):.6e} N"
    )
    print("=" * 108)


def plot_results(
    response: NonlinearPulsatingResponse,
) -> None:
    """Plot tower and fixed-fiber cyclic responses."""
    normalized_time = (
        response.loading.times / response.loading.period
    )

    plt.figure()
    plt.plot(
        response.top_displacements,
        response.loading.forces / 1.0e6,
        marker="o",
        markersize=3,
    )
    plt.xlabel("Tower-top horizontal displacement (m)")
    plt.ylabel("Tower-top horizontal force (MN)")
    plt.title("Tower force-displacement response")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        response.maximum_absolute_plastic_strains,
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Maximum absolute plastic strain")
    plt.title("Tower maximum plastic-strain history")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        response.maximum_damages,
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Maximum damage")
    plt.title("Tower maximum damage history")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        response.critical_fiber_strains,
        response.critical_fiber_stresses,
        marker="o",
        markersize=3,
    )
    plt.xlabel("Critical-fiber total strain")
    plt.ylabel("Critical-fiber stress (MPa)")
    plt.title("Fixed critical-fiber stress-strain response")
    plt.grid(True)
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Run one nonlinear positive-mean pulsating "
            "cycle of the NREL 5 MW tower."
        )
    )
    parser.add_argument(
        "--maximum-force",
        type=float,
        default=1.0e6,
        help="Maximum tower-top horizontal force in newtons.",
    )
    parser.add_argument(
        "--force-ratio",
        type=float,
        default=0.1,
        help="Positive force ratio F_min/F_max.",
    )
    parser.add_argument(
        "--period",
        type=float,
        default=10.0,
        help="Duration of one pulsating cycle.",
    )
    parser.add_argument(
        "--increments",
        type=int,
        default=20,
        help=(
            "Equal increments per cycle; must be divisible by four."
        ),
    )
    parser.add_argument(
        "--elements",
        type=int,
        default=10,
        help="Number of uniform tower beam elements.",
    )
    parser.add_argument(
        "--gauss",
        type=int,
        default=2,
        help="Gauss points per beam element.",
    )
    parser.add_argument(
        "--circumferential",
        type=int,
        default=16,
        help="Circumferential fibers per annular section.",
    )
    parser.add_argument(
        "--radial",
        type=int,
        default=1,
        help="Radial fiber layers through the tower wall.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Run without opening figures.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the nonlinear pulsating tower example."""
    arguments = parse_arguments()

    configuration = TowerConfiguration(
        horizontal_force=arguments.maximum_force,
        n_elements=arguments.elements,
        n_gauss=arguments.gauss,
        n_circumferential=arguments.circumferential,
        n_radial=arguments.radial,
    )
    material = MaterialParameters()
    loading = create_pulsating_top_force_history(
        maximum_force=arguments.maximum_force,
        force_ratio=arguments.force_ratio,
        period=arguments.period,
        n_cycles=1,
        increments_per_cycle=arguments.increments,
    )

    response = run_nonlinear_pulsating_analysis(
        configuration=configuration,
        material=material,
        loading=loading,
    )
    print_summary(
        configuration=configuration,
        material=material,
        response=response,
    )

    if not arguments.no_plot:
        plot_results(response)


if __name__ == "__main__":
    main()

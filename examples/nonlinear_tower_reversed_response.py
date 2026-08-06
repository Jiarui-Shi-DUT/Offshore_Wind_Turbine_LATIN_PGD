# -*- coding: utf-8 -*-
"""
Nonlinear fully reversed response of the NREL 5 MW tower.

The tower-top force is

    F(t) = F_a * sin(2*pi*t/T),

so one complete cycle follows

    0 -> +F_a -> 0 -> -F_a -> 0.

The example records the global tower response and the complete internal
state of every fiber. It is intended to verify:

    - force and stress reversal;
    - alternating tensile and compressive response at one fixed fiber;
    - reversal of the plastic-flow direction;
    - kinematic hardening and backstress evolution;
    - unilateral damage response;
    - residual tower displacement after one complete reversed cycle.

Use:
    python -m examples.nonlinear_tower_reversed_response --no-plot
    python -m examples.nonlinear_tower_reversed_response
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
    ReversedTopForceHistory,
    create_reversed_top_force_history,
)
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
    ViscoplasticTowerResponse,
    solve_nonlinear_tower_load_step,
)
from material.viscoplastic_damage_1d import (
    MaterialParameters,
    evaluate_state,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
FiberLocation = Tuple[int, int, int]


@dataclass(frozen=True)
class NonlinearReversedResponse:
    """Global and fiber histories from a fully reversed tower analysis."""

    loading: ReversedTopForceHistory
    material: MaterialParameters
    analysis_times: FloatArray
    top_displacements: FloatArray
    top_rotations: FloatArray
    base_horizontal_reactions: FloatArray
    base_moment_reactions: FloatArray
    iterations: IntArray
    residual_norms: FloatArray
    fiber_strains: FloatArray
    fiber_stresses: FloatArray
    fiber_states: FloatArray
    critical_location: FiberLocation
    critical_height: float
    critical_y_coordinate: float

    def __post_init__(self) -> None:
        if not isinstance(
            self.loading,
            ReversedTopForceHistory,
        ):
            raise TypeError(
                "loading must be a ReversedTopForceHistory."
            )
        if not isinstance(
            self.material,
            MaterialParameters,
        ):
            raise TypeError(
                "material must be a MaterialParameters object."
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

        global_arrays = (
            analysis_times,
            top_displacements,
            top_rotations,
            base_horizontal_reactions,
            base_moment_reactions,
            iterations,
            residual_norms,
        )
        for array in global_arrays:
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
        fiber_states = np.asarray(
            self.fiber_states,
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
        if fiber_stresses.shape != fiber_strains.shape:
            raise ValueError(
                "fiber_stresses must match fiber_strains."
            )
        if fiber_states.shape != fiber_strains.shape + (4,):
            raise ValueError(
                "fiber_states must append four material variables."
            )

        for array in (
            analysis_times,
            top_displacements,
            top_rotations,
            base_horizontal_reactions,
            base_moment_reactions,
            residual_norms,
            fiber_strains,
            fiber_stresses,
            fiber_states,
        ):
            if np.any(~np.isfinite(array)):
                raise ValueError(
                    "All response histories must be finite."
                )

        if np.any(np.diff(analysis_times) <= 0.0):
            raise ValueError(
                "analysis_times must be strictly increasing."
            )
        if np.any(iterations < 1):
            raise ValueError(
                "Every Newton iteration count must be positive."
            )
        if np.any(residual_norms < 0.0):
            raise ValueError(
                "residual_norms must be non-negative."
            )
        if np.any(fiber_states[..., 3] < 0.0):
            raise ValueError(
                "Fiber damage must be non-negative."
            )
        if np.any(fiber_states[..., 3] >= 1.0):
            raise ValueError(
                "Fiber damage must be smaller than one."
            )

        location = tuple(
            int(value) for value in self.critical_location
        )
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
            raise ValueError(
                "critical_height must be finite."
            )
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
            "fiber_states",
            fiber_states.copy(),
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
    def fiber_plastic_strains(self) -> FloatArray:
        """Return all plastic-strain histories."""
        return self.fiber_states[..., 0].copy()

    @property
    def fiber_alphas(self) -> FloatArray:
        """Return all kinematic internal-variable histories."""
        return self.fiber_states[..., 1].copy()

    @property
    def fiber_r_bars(self) -> FloatArray:
        """Return all isotropic internal-variable histories."""
        return self.fiber_states[..., 2].copy()

    @property
    def fiber_damages(self) -> FloatArray:
        """Return all damage histories."""
        return self.fiber_states[..., 3].copy()

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
            np.abs(self.fiber_states[..., 0]),
            axis=(1, 2, 3),
        )

    @property
    def maximum_damages(self) -> FloatArray:
        """Return maximum fiber damage at every time point."""
        return np.max(
            self.fiber_states[..., 3],
            axis=(1, 2, 3),
        )

    def _critical_history(
        self,
        values: FloatArray,
    ) -> FloatArray:
        """Return one history at the fixed critical location."""
        element_index, gauss_index, fiber_index = (
            self.critical_location
        )
        return values[
            :,
            element_index,
            gauss_index,
            fiber_index,
        ].copy()

    @property
    def critical_fiber_strains(self) -> FloatArray:
        """Return total-strain history at the fixed critical fiber."""
        return self._critical_history(self.fiber_strains)

    @property
    def critical_fiber_stresses(self) -> FloatArray:
        """Return stress history at the fixed critical fiber."""
        return self._critical_history(self.fiber_stresses)

    @property
    def critical_fiber_states(self) -> FloatArray:
        """Return all four state histories at the critical fiber."""
        element_index, gauss_index, fiber_index = (
            self.critical_location
        )
        return self.fiber_states[
            :,
            element_index,
            gauss_index,
            fiber_index,
            :,
        ].copy()

    @property
    def critical_fiber_plastic_strains(self) -> FloatArray:
        """Return plastic-strain history at the critical fiber."""
        return self.critical_fiber_states[:, 0]

    @property
    def critical_fiber_alphas(self) -> FloatArray:
        """Return alpha history at the critical fiber."""
        return self.critical_fiber_states[:, 1]

    @property
    def critical_fiber_r_bars(self) -> FloatArray:
        """Return r_bar history at the critical fiber."""
        return self.critical_fiber_states[:, 2]

    @property
    def critical_fiber_damages(self) -> FloatArray:
        """Return damage history at the critical fiber."""
        return self.critical_fiber_states[:, 3]

    @property
    def critical_fiber_backstresses(self) -> FloatArray:
        """Return beta = C * alpha at the critical fiber."""
        return (
            self.material.C
            * self.critical_fiber_alphas
        )

    def _critical_evaluated_quantity(
        self,
        result_index: int,
    ) -> FloatArray:
        """Evaluate one material response quantity at every time point."""
        values = np.empty(
            self.loading.n_time_points,
            dtype=np.float64,
        )
        strains = self.critical_fiber_strains
        states = self.critical_fiber_states

        for time_index in range(self.loading.n_time_points):
            response = evaluate_state(
                total_strain=float(strains[time_index]),
                state=states[time_index],
                material=self.material,
            )
            values[time_index] = float(
                response[result_index]
            )
        return values

    @property
    def critical_fiber_effective_relative_stresses(
        self,
    ) -> FloatArray:
        """Return sigma/(1-D)-beta at the critical fiber."""
        return self._critical_evaluated_quantity(4)

    @property
    def critical_fiber_yield_functions(self) -> FloatArray:
        """Return the yield-function history at the critical fiber."""
        return self._critical_evaluated_quantity(5)

    @property
    def critical_fiber_energy_release_rates(
        self,
    ) -> FloatArray:
        """Return the damage energy-release-rate history."""
        return self._critical_evaluated_quantity(6)


def snapshot_fiber_fields(
    response: ViscoplasticTowerResponse,
) -> Tuple[FloatArray, FloatArray, FloatArray]:
    """Stack fiber strains, stresses, and complete states."""
    strains = []
    stresses = []
    states = []

    for element_response in response.element_responses:
        element_strains = []
        element_stresses = []
        element_states = []

        for section_response in element_response.section_responses:
            element_strains.append(
                section_response.fiber_strains
            )
            element_stresses.append(
                section_response.fiber_stresses
            )
            element_states.append(
                section_response.fiber_states
            )

        strains.append(
            np.stack(element_strains, axis=0)
        )
        stresses.append(
            np.stack(element_stresses, axis=0)
        )
        states.append(
            np.stack(element_states, axis=0)
        )

    return (
        np.stack(strains, axis=0),
        np.stack(stresses, axis=0),
        np.stack(states, axis=0),
    )


def locate_critical_fiber(
    fiber_states: FloatArray,
) -> FiberLocation:
    """
    Select one fixed critical fiber for reversed-path diagnostics.

    Peak damage over the complete history is the primary criterion.
    If no damage occurs, peak absolute plastic strain is used.
    """
    damages = fiber_states[..., 3]
    peak_damage = np.max(damages, axis=0)

    if float(np.max(peak_damage)) > 0.0:
        flat_index = int(np.argmax(peak_damage))
        location = np.unravel_index(
            flat_index,
            peak_damage.shape,
        )
    else:
        plastic_strains = fiber_states[..., 0]
        peak_plastic = np.max(
            np.abs(plastic_strains),
            axis=0,
        )
        flat_index = int(np.argmax(peak_plastic))
        location = np.unravel_index(
            flat_index,
            peak_plastic.shape,
        )

    return (
        int(location[0]),
        int(location[1]),
        int(location[2]),
    )


def run_nonlinear_reversed_analysis(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    loading: ReversedTopForceHistory,
    max_iterations: int = 40,
) -> NonlinearReversedResponse:
    """Build the tower and advance the complete reversed history."""
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
    state_snapshots = []

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
            fiber_states,
        ) = snapshot_fiber_fields(solution.response)

        strain_snapshots.append(fiber_strains)
        stress_snapshots.append(fiber_stresses)
        state_snapshots.append(fiber_states)

    fiber_strain_history = np.stack(
        strain_snapshots,
        axis=0,
    )
    fiber_stress_history = np.stack(
        stress_snapshots,
        axis=0,
    )
    fiber_state_history = np.stack(
        state_snapshots,
        axis=0,
    )

    critical_location = locate_critical_fiber(
        fiber_states=fiber_state_history,
    )
    element_index, gauss_index, fiber_index = (
        critical_location
    )
    critical_element = system.elements[element_index]
    critical_section = critical_element.sections[gauss_index]

    return NonlinearReversedResponse(
        loading=loading,
        material=material,
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
        fiber_states=fiber_state_history,
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
    response: NonlinearReversedResponse,
) -> None:
    """Print loading, reversal, and nonlinear-response results."""
    loading = response.loading
    quarter = loading.increments_per_cycle // 4
    indices = (
        0,
        quarter,
        2 * quarter,
        3 * quarter,
        4 * quarter,
    )

    maximum_stresses = (
        response.maximum_absolute_stresses
    )
    maximum_plastic_strains = (
        response.maximum_absolute_plastic_strains
    )
    maximum_damages = response.maximum_damages
    critical_stresses = (
        response.critical_fiber_stresses
    )
    critical_plastic_strains = (
        response.critical_fiber_plastic_strains
    )
    critical_backstresses = (
        response.critical_fiber_backstresses
    )
    critical_yield_functions = (
        response.critical_fiber_yield_functions
    )

    print("=" * 132)
    print(
        "Nonlinear NREL 5 MW tower under fully reversed loading"
    )
    print("=" * 132)
    print(
        f"F_a={loading.force_amplitude / 1.0e6:.6g} MN, "
        f"F_max={loading.maximum_force / 1.0e6:.6g} MN, "
        f"F_min={loading.minimum_force / 1.0e6:.6g} MN"
    )
    print(
        f"Period={loading.period:.6g}, "
        f"cycles={loading.n_cycles}, "
        f"increments/cycle={loading.increments_per_cycle}"
    )
    print(
        f"Material: E={response.material.E:.6g} MPa, "
        f"sigma_y={response.material.sigma_y:.6g} MPa, "
        f"C={response.material.C:.6g} MPa"
    )
    print(
        "Discretisation: "
        f"{configuration.n_elements} elements, "
        f"{configuration.n_gauss} Gauss points/element, "
        f"{configuration.n_circumferential} x "
        f"{configuration.n_radial} fibers/section"
    )
    print("-" * 132)
    print(
        " t/T   force(MN)   top ux(m)   max|stress|   "
        "critical stress   critical beta   critical eps_p   "
        "critical f      max D   Newton"
    )

    for index in indices:
        print(
            f"{loading.times[index] / loading.period:5.2f}  "
            f"{loading.forces[index] / 1.0e6:10.4f}  "
            f"{response.top_displacements[index]:10.4e}  "
            f"{maximum_stresses[index]:11.4e}  "
            f"{critical_stresses[index]:15.6e}  "
            f"{critical_backstresses[index]:13.6e}  "
            f"{critical_plastic_strains[index]:14.6e}  "
            f"{critical_yield_functions[index]:11.4e}  "
            f"{maximum_damages[index]:8.3e}  "
            f"{response.iterations[index]:6d}"
        )

    positive_peak = quarter
    negative_peak = 3 * quarter
    stress_reversal = bool(
        critical_stresses[positive_peak]
        * critical_stresses[negative_peak]
        < 0.0
    )
    positive_plastic_increment = float(
        critical_plastic_strains[positive_peak]
        - critical_plastic_strains[0]
    )
    negative_plastic_increment = float(
        critical_plastic_strains[negative_peak]
        - critical_plastic_strains[2 * quarter]
    )
    plastic_direction_reversal = bool(
        positive_plastic_increment
        * negative_plastic_increment
        < 0.0
    )

    element_index, gauss_index, fiber_index = (
        response.critical_location
    )
    print("-" * 132)
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
        "Stress reversal between positive and negative peaks: "
        f"{stress_reversal}"
    )
    print(
        "Plastic-flow direction reversal over the two half-cycles: "
        f"{plastic_direction_reversal}"
    )
    print(
        "Positive-half plastic increment: "
        f"{positive_plastic_increment:.6e}"
    )
    print(
        "Negative-half plastic increment: "
        f"{negative_plastic_increment:.6e}"
    )
    print(
        "Final tower-top residual displacement: "
        f"{response.top_displacements[-1]:.6e} m"
    )
    print(
        "Final maximum absolute plastic strain: "
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
    print("=" * 132)


def plot_results(
    response: NonlinearReversedResponse,
) -> None:
    """Plot global and fixed-fiber reversed responses."""
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
    plt.title("Tower fully reversed force-displacement response")
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
    plt.title("Fixed critical-fiber reversed stress-strain response")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        response.critical_fiber_backstresses,
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Critical-fiber backstress, beta (MPa)")
    plt.title("Kinematic-hardening backstress history")
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

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Run fully reversed nonlinear cycles "
            "of the NREL 5 MW tower."
        )
    )
    parser.add_argument(
        "--force-amplitude",
        type=float,
        default=1.0e6,
        help="Tower-top force amplitude in newtons.",
    )
    parser.add_argument(
        "--period",
        type=float,
        default=10.0,
        help="Duration of one reversed cycle.",
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=1,
        help="Number of complete reversed cycles.",
    )
    parser.add_argument(
        "--increments",
        type=int,
        default=40,
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
    """Run the nonlinear fully reversed tower example."""
    arguments = parse_arguments()

    configuration = TowerConfiguration(
        horizontal_force=arguments.force_amplitude,
        n_elements=arguments.elements,
        n_gauss=arguments.gauss,
        n_circumferential=arguments.circumferential,
        n_radial=arguments.radial,
    )
    material = MaterialParameters()
    loading = create_reversed_top_force_history(
        force_amplitude=arguments.force_amplitude,
        period=arguments.period,
        n_cycles=arguments.cycles,
        increments_per_cycle=arguments.increments,
    )

    response = run_nonlinear_reversed_analysis(
        configuration=configuration,
        material=material,
        loading=loading,
    )
    print_summary(
        configuration=configuration,
        response=response,
    )

    if not arguments.no_plot:
        plot_results(response)


if __name__ == "__main__":
    main()

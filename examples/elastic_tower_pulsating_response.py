# -*- coding: utf-8 -*-
"""
Linear-elastic response of the tapered tower under a pulsating top force.

The efficient time history is obtained by scaling one reference solution at
F_max. Direct finite-element solutions are also performed at five checkpoint
times to verify that displacement, curvature, and fiber stress scale linearly
with the applied force.

Use:
    python -m examples.elastic_tower_pulsating_response --no-plot
    python -m examples.elastic_tower_pulsating_response
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
    run_analysis,
)
from fem.tower_loading import (
    PulsatingTopForceHistory,
    create_pulsating_top_force_history,
)
from fem.tower_response_2d import (
    ElasticTowerResponse,
    recover_elastic_tower_response,
)


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ElasticPulsatingResponse:
    """Selected global and local elastic response histories."""

    loading: PulsatingTopForceHistory
    top_displacements: FloatArray
    top_rotations: FloatArray
    critical_curvatures: FloatArray
    critical_minimum_stresses: FloatArray
    critical_maximum_stresses: FloatArray
    critical_element_index: int
    critical_gauss_index: int
    critical_height: float

    def __post_init__(self) -> None:
        expected_shape = self.loading.times.shape

        arrays = (
            np.asarray(self.top_displacements, dtype=np.float64),
            np.asarray(self.top_rotations, dtype=np.float64),
            np.asarray(self.critical_curvatures, dtype=np.float64),
            np.asarray(
                self.critical_minimum_stresses,
                dtype=np.float64,
            ),
            np.asarray(
                self.critical_maximum_stresses,
                dtype=np.float64,
            ),
        )

        for array in arrays:
            if array.shape != expected_shape:
                raise ValueError(
                    "Every response history must match "
                    "the loading time vector."
                )
            if np.any(~np.isfinite(array)):
                raise ValueError(
                    "Every response history must be finite."
                )

        if self.critical_element_index < 0:
            raise ValueError(
                "critical_element_index must be non-negative."
            )
        if self.critical_gauss_index < 0:
            raise ValueError(
                "critical_gauss_index must be non-negative."
            )
        if not np.isfinite(self.critical_height):
            raise ValueError("critical_height must be finite.")

        object.__setattr__(self, "top_displacements", arrays[0])
        object.__setattr__(self, "top_rotations", arrays[1])
        object.__setattr__(self, "critical_curvatures", arrays[2])
        object.__setattr__(
            self,
            "critical_minimum_stresses",
            arrays[3],
        )
        object.__setattr__(
            self,
            "critical_maximum_stresses",
            arrays[4],
        )


@dataclass(frozen=True)
class LinearityCheck:
    """Maximum relative errors from direct checkpoint solutions."""

    top_displacement_error: float
    top_rotation_error: float
    curvature_error: float
    minimum_stress_error: float
    maximum_stress_error: float


def configuration_with_force(
    configuration: TowerConfiguration,
    horizontal_force: float,
) -> TowerConfiguration:
    """Copy a tower configuration while replacing the top force."""
    return TowerConfiguration(
        height=configuration.height,
        base_outer_diameter=configuration.base_outer_diameter,
        top_outer_diameter=configuration.top_outer_diameter,
        base_thickness=configuration.base_thickness,
        top_thickness=configuration.top_thickness,
        elastic_modulus=configuration.elastic_modulus,
        horizontal_force=float(horizontal_force),
        n_elements=configuration.n_elements,
        n_gauss=configuration.n_gauss,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )


def critical_response_indices(
    response: ElasticTowerResponse,
) -> Tuple[int, int]:
    """Locate the Gauss point with maximum absolute fiber stress."""
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


def build_scaled_history(
    configuration: TowerConfiguration,
    loading: PulsatingTopForceHistory,
) -> ElasticPulsatingResponse:
    """Build the exact linear-elastic response by scaling the F_max solution."""
    reference_configuration = configuration_with_force(
        configuration=configuration,
        horizontal_force=loading.maximum_force,
    )
    _, tower_geometry, assembly, solution = run_analysis(
        reference_configuration
    )
    response = recover_elastic_tower_response(
        assembly=assembly,
        solution=solution,
        tower_geometry=tower_geometry,
        elastic_modulus=configuration.elastic_modulus,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )

    element_index, gauss_index = critical_response_indices(response)
    scale = loading.forces / loading.maximum_force

    reference_top_displacement = float(solution.displacements[-3])
    reference_top_rotation = float(solution.displacements[-1])
    reference_curvature = float(
        response.curvatures[element_index, gauss_index]
    )
    reference_minimum_stress = float(
        response.minimum_fiber_stresses[element_index, gauss_index]
    )
    reference_maximum_stress = float(
        response.maximum_fiber_stresses[element_index, gauss_index]
    )
    critical_height = float(
        response.gauss_heights[element_index, gauss_index]
    )

    return ElasticPulsatingResponse(
        loading=loading,
        top_displacements=scale * reference_top_displacement,
        top_rotations=scale * reference_top_rotation,
        critical_curvatures=scale * reference_curvature,
        critical_minimum_stresses=scale * reference_minimum_stress,
        critical_maximum_stresses=scale * reference_maximum_stress,
        critical_element_index=element_index,
        critical_gauss_index=gauss_index,
        critical_height=critical_height,
    )


def relative_error(actual: float, expected: float) -> float:
    """Return an absolute relative error with a safe denominator."""
    return float(abs(actual - expected) / max(1.0, abs(expected)))


def verify_direct_checkpoint_solutions(
    configuration: TowerConfiguration,
    history: ElasticPulsatingResponse,
) -> LinearityCheck:
    """Compare scaled histories with independent FE solves at five times."""
    increments_per_cycle = history.loading.increments_per_cycle
    checkpoint_indices = np.array(
        [
            0,
            increments_per_cycle // 4,
            increments_per_cycle // 2,
            3 * increments_per_cycle // 4,
            increments_per_cycle,
        ],
        dtype=np.int64,
    )

    displacement_errors = []
    rotation_errors = []
    curvature_errors = []
    minimum_stress_errors = []
    maximum_stress_errors = []

    element_index = history.critical_element_index
    gauss_index = history.critical_gauss_index

    for time_index in checkpoint_indices:
        force = float(history.loading.forces[time_index])
        direct_configuration = configuration_with_force(
            configuration=configuration,
            horizontal_force=force,
        )
        _, tower_geometry, assembly, solution = run_analysis(
            direct_configuration
        )
        response = recover_elastic_tower_response(
            assembly=assembly,
            solution=solution,
            tower_geometry=tower_geometry,
            elastic_modulus=configuration.elastic_modulus,
            n_circumferential=configuration.n_circumferential,
            n_radial=configuration.n_radial,
        )

        displacement_errors.append(
            relative_error(
                actual=float(solution.displacements[-3]),
                expected=float(history.top_displacements[time_index]),
            )
        )
        rotation_errors.append(
            relative_error(
                actual=float(solution.displacements[-1]),
                expected=float(history.top_rotations[time_index]),
            )
        )
        curvature_errors.append(
            relative_error(
                actual=float(
                    response.curvatures[element_index, gauss_index]
                ),
                expected=float(
                    history.critical_curvatures[time_index]
                ),
            )
        )
        minimum_stress_errors.append(
            relative_error(
                actual=float(
                    response.minimum_fiber_stresses[
                        element_index,
                        gauss_index,
                    ]
                ),
                expected=float(
                    history.critical_minimum_stresses[time_index]
                ),
            )
        )
        maximum_stress_errors.append(
            relative_error(
                actual=float(
                    response.maximum_fiber_stresses[
                        element_index,
                        gauss_index,
                    ]
                ),
                expected=float(
                    history.critical_maximum_stresses[time_index]
                ),
            )
        )

    return LinearityCheck(
        top_displacement_error=max(displacement_errors),
        top_rotation_error=max(rotation_errors),
        curvature_error=max(curvature_errors),
        minimum_stress_error=max(minimum_stress_errors),
        maximum_stress_error=max(maximum_stress_errors),
    )


def print_summary(
    configuration: TowerConfiguration,
    history: ElasticPulsatingResponse,
    check: LinearityCheck,
) -> None:
    """Print loading, critical-point, and linearity-check results."""
    loading = history.loading
    quarter = loading.increments_per_cycle // 4
    checkpoint_indices = (
        0,
        quarter,
        2 * quarter,
        3 * quarter,
        4 * quarter,
    )

    print("=" * 78)
    print("Elastic tapered tower under pulsating horizontal force")
    print("=" * 78)
    print(f"Tower height: {configuration.height:.6g} m")
    print(f"Beam elements: {configuration.n_elements}")
    print(f"Force ratio F_min/F_max: {loading.force_ratio:.6g}")
    print(f"Maximum force: {loading.maximum_force:.6e} N")
    print(f"Minimum force: {loading.minimum_force:.6e} N")
    print(f"Mean force: {loading.mean_force:.6e} N")
    print(f"Force amplitude: {loading.force_amplitude:.6e} N")
    print(f"Period: {loading.period:.6e} s")
    print(
        "Time increments per cycle: "
        f"{loading.increments_per_cycle}"
    )
    print("-" * 78)
    print(
        "Critical Gauss point: element "
        f"{history.critical_element_index + 1}, point "
        f"{history.critical_gauss_index + 1}"
    )
    print(
        "Critical height above base: "
        f"{history.critical_height:.6e} m"
    )
    print("-" * 78)
    print(
        " time/T       force (N)       top u_x (m)   "
        "curvature (1/m)   max stress (MPa)"
    )
    for index in checkpoint_indices:
        normalized_time = loading.times[index] / loading.period
        print(
            f" {normalized_time:6.2f}   "
            f"{loading.forces[index]:13.6e}   "
            f"{history.top_displacements[index]:13.6e}   "
            f"{history.critical_curvatures[index]:13.6e}   "
            f"{history.critical_maximum_stresses[index] / 1.0e6:13.6e}"
        )
    print("-" * 78)
    print("Maximum relative errors against independent FE solves")
    print(
        "Top displacement: "
        f"{check.top_displacement_error:.6e}"
    )
    print(
        "Top rotation: "
        f"{check.top_rotation_error:.6e}"
    )
    print(f"Critical curvature: {check.curvature_error:.6e}")
    print(
        "Critical minimum stress: "
        f"{check.minimum_stress_error:.6e}"
    )
    print(
        "Critical maximum stress: "
        f"{check.maximum_stress_error:.6e}"
    )
    print("=" * 78)


def plot_history(history: ElasticPulsatingResponse) -> None:
    """Plot force, top displacement, and critical fiber stresses."""
    normalized_time = history.loading.times / history.loading.period

    plt.figure()
    plt.plot(
        normalized_time,
        history.loading.forces / 1.0e6,
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Top horizontal force (MN)")
    plt.title("Pulsating tower-top force")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        history.top_displacements,
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Tower-top horizontal displacement (m)")
    plt.title("Elastic tower-top displacement history")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        history.critical_minimum_stresses / 1.0e6,
        label="Minimum fiber stress",
    )
    plt.plot(
        normalized_time,
        history.critical_maximum_stresses / 1.0e6,
        label="Maximum fiber stress",
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Fiber stress (MPa)")
    plt.title("Critical-section elastic fiber-stress history")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Run the elastic tapered tower under one "
            "positive-mean pulsating top-force cycle."
        )
    )
    parser.add_argument(
        "--elements",
        type=int,
        default=40,
        help="Number of uniform beam elements.",
    )
    parser.add_argument(
        "--maximum-force",
        type=float,
        default=1.0e6,
        help="Maximum top horizontal force in newtons.",
    )
    parser.add_argument(
        "--force-ratio",
        type=float,
        default=0.1,
        help="Force ratio F_min/F_max.",
    )
    parser.add_argument(
        "--period",
        type=float,
        default=1.0,
        help="Force-cycle period.",
    )
    parser.add_argument(
        "--increments",
        type=int,
        default=100,
        help="Equal time increments per cycle.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Run without opening figures.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the pulsating linear-elastic tower example."""
    arguments = parse_arguments()

    configuration = TowerConfiguration(
        horizontal_force=arguments.maximum_force,
        n_elements=arguments.elements,
    )
    loading = create_pulsating_top_force_history(
        maximum_force=arguments.maximum_force,
        force_ratio=arguments.force_ratio,
        period=arguments.period,
        n_cycles=1,
        increments_per_cycle=arguments.increments,
    )
    history = build_scaled_history(
        configuration=configuration,
        loading=loading,
    )
    check = verify_direct_checkpoint_solutions(
        configuration=configuration,
        history=history,
    )

    print_summary(
        configuration=configuration,
        history=history,
        check=check,
    )

    if not arguments.no_plot:
        plot_history(history)


if __name__ == "__main__":
    main()

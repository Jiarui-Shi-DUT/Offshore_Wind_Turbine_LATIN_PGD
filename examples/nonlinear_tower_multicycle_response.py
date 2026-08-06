# -*- coding: utf-8 -*-
"""
Multi-cycle nonlinear response of the NREL 5 MW tower.

This example reuses the verified fully reversed nonlinear solver and the
per-cycle diagnostics module. It is intended to investigate:

    - cycle-by-cycle peak and residual tower response;
    - cyclic hardening or softening;
    - plastic-strain and backstress evolution;
    - damage accumulation;
    - force-displacement work evolution;
    - approach to a cycle-similar response.

The loading is

    F(t) = F_a * sin(2*pi*t/T),

and every cycle follows

    0 -> +F_a -> 0 -> -F_a -> 0.

The cycle-similarity indicator is a numerical diagnostic, not a proof of a
true stabilized material limit cycle. It compares two consecutive cycles in
terms of displacement range, fixed-fiber stress range, external-work
magnitude, normalized residual-displacement drift, and normalized
plastic-strain drift.

Use:
    python -m examples.nonlinear_tower_multicycle_response --no-plot
    python -m examples.nonlinear_tower_multicycle_response
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_multicycle_diagnostics import (
    CycleDiagnostics,
    MulticycleDiagnostics,
    extract_multicycle_diagnostics,
)
from examples.nonlinear_tower_reversed_response import (
    NonlinearReversedResponse,
    run_nonlinear_reversed_analysis,
)
from fem.tower_loading import (
    create_reversed_top_force_history,
)
from material.viscoplastic_damage_1d import MaterialParameters


@dataclass(frozen=True)
class ConsecutiveCycleSimilarity:
    """Dimensionless comparison between two consecutive cycles."""

    previous_cycle: int
    current_cycle: int
    displacement_range_change: float
    stress_range_change: float
    external_work_change: float
    normalized_residual_drift: float
    normalized_plastic_strain_drift: float
    maximum_indicator: float
    tolerance: float
    within_tolerance: bool

    def __post_init__(self) -> None:
        integer_names = (
            "previous_cycle",
            "current_cycle",
        )
        for name in integer_names:
            value = getattr(self, name)
            if isinstance(value, bool):
                raise TypeError(name + " must be an integer.")
            converted = int(value)
            if converted != value:
                raise TypeError(name + " must be an integer.")
            if converted < 1:
                raise ValueError(
                    name + " must be at least one."
                )
            object.__setattr__(self, name, converted)

        if self.current_cycle != self.previous_cycle + 1:
            raise ValueError(
                "current_cycle must immediately follow previous_cycle."
            )

        metric_names = (
            "displacement_range_change",
            "stress_range_change",
            "external_work_change",
            "normalized_residual_drift",
            "normalized_plastic_strain_drift",
            "maximum_indicator",
            "tolerance",
        )
        for name in metric_names:
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(name + " must be finite.")
            if value < 0.0:
                raise ValueError(
                    name + " must be non-negative."
                )
            object.__setattr__(self, name, value)

        expected_maximum = max(
            self.displacement_range_change,
            self.stress_range_change,
            self.external_work_change,
            self.normalized_residual_drift,
            self.normalized_plastic_strain_drift,
        )
        if not np.isclose(
            self.maximum_indicator,
            expected_maximum,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "maximum_indicator is inconsistent "
                "with component indicators."
            )
        if self.tolerance <= 0.0:
            raise ValueError("tolerance must be positive.")

        expected_flag = bool(
            self.maximum_indicator <= self.tolerance
        )
        if bool(self.within_tolerance) != expected_flag:
            raise ValueError(
                "within_tolerance is inconsistent "
                "with maximum_indicator and tolerance."
            )
        object.__setattr__(
            self,
            "within_tolerance",
            expected_flag,
        )


@dataclass(frozen=True)
class MulticycleTowerResult:
    """Full nonlinear history and cycle-by-cycle diagnostics."""

    response: NonlinearReversedResponse
    diagnostics: MulticycleDiagnostics
    similarities: Tuple[ConsecutiveCycleSimilarity, ...]
    similarity_tolerance: float

    def __post_init__(self) -> None:
        if not isinstance(
            self.response,
            NonlinearReversedResponse,
        ):
            raise TypeError(
                "response must be a NonlinearReversedResponse."
            )
        if not isinstance(
            self.diagnostics,
            MulticycleDiagnostics,
        ):
            raise TypeError(
                "diagnostics must be a MulticycleDiagnostics object."
            )

        if (
            self.diagnostics.n_cycles
            != self.response.loading.n_cycles
        ):
            raise ValueError(
                "diagnostics must cover every loading cycle."
            )

        similarities = tuple(self.similarities)
        expected_count = max(
            0,
            self.diagnostics.n_cycles - 1,
        )
        if len(similarities) != expected_count:
            raise ValueError(
                "similarities must compare every consecutive "
                "cycle pair."
            )
        for index, similarity in enumerate(
            similarities,
            start=2,
        ):
            if not isinstance(
                similarity,
                ConsecutiveCycleSimilarity,
            ):
                raise TypeError(
                    "Every similarity must be a "
                    "ConsecutiveCycleSimilarity object."
                )
            if similarity.previous_cycle != index - 1:
                raise ValueError(
                    "Similarity comparisons are not ordered."
                )
            if similarity.current_cycle != index:
                raise ValueError(
                    "Similarity comparisons are not ordered."
                )

        similarity_tolerance = float(
            self.similarity_tolerance
        )
        if not np.isfinite(similarity_tolerance):
            raise ValueError(
                "similarity_tolerance must be finite."
            )
        if similarity_tolerance <= 0.0:
            raise ValueError(
                "similarity_tolerance must be positive."
            )
        for similarity in similarities:
            if not np.isclose(
                similarity.tolerance,
                similarity_tolerance,
                rtol=0.0,
                atol=0.0,
            ):
                raise ValueError(
                    "Every similarity must use "
                    "similarity_tolerance."
                )

        object.__setattr__(
            self,
            "similarities",
            similarities,
        )
        object.__setattr__(
            self,
            "similarity_tolerance",
            similarity_tolerance,
        )

    @property
    def final_cycle_is_within_tolerance(self) -> bool:
        """Return the final consecutive-cycle comparison flag."""
        if len(self.similarities) == 0:
            return False
        return bool(
            self.similarities[-1].within_tolerance
        )


def relative_change(
    current: float,
    previous: float,
    scale_floor: float = 1.0e-14,
) -> float:
    """Return an absolute relative change with a protected denominator."""
    current = float(current)
    previous = float(previous)
    scale_floor = float(scale_floor)

    if not np.isfinite(current):
        raise ValueError("current must be finite.")
    if not np.isfinite(previous):
        raise ValueError("previous must be finite.")
    if not np.isfinite(scale_floor):
        raise ValueError("scale_floor must be finite.")
    if scale_floor <= 0.0:
        raise ValueError("scale_floor must be positive.")

    denominator = max(
        abs(current),
        abs(previous),
        scale_floor,
    )
    return float(
        abs(current - previous) / denominator
    )


def normalized_drift(
    current: float,
    previous: float,
    cycle_scale: float,
    scale_floor: float = 1.0e-14,
) -> float:
    """Return an absolute endpoint drift normalized by a cycle range."""
    current = float(current)
    previous = float(previous)
    cycle_scale = float(cycle_scale)
    scale_floor = float(scale_floor)

    for value, name in (
        (current, "current"),
        (previous, "previous"),
        (cycle_scale, "cycle_scale"),
        (scale_floor, "scale_floor"),
    ):
        if not np.isfinite(value):
            raise ValueError(name + " must be finite.")

    if cycle_scale < 0.0:
        raise ValueError(
            "cycle_scale must be non-negative."
        )
    if scale_floor <= 0.0:
        raise ValueError(
            "scale_floor must be positive."
        )

    denominator = max(
        abs(cycle_scale),
        scale_floor,
    )
    return float(
        abs(current - previous) / denominator
    )


def compare_consecutive_cycles(
    previous: CycleDiagnostics,
    current: CycleDiagnostics,
    tolerance: float = 1.0e-3,
) -> ConsecutiveCycleSimilarity:
    """Compare two adjacent cycles using five dimensionless indicators."""
    if not isinstance(previous, CycleDiagnostics):
        raise TypeError(
            "previous must be a CycleDiagnostics object."
        )
    if not isinstance(current, CycleDiagnostics):
        raise TypeError(
            "current must be a CycleDiagnostics object."
        )
    if current.cycle_number != previous.cycle_number + 1:
        raise ValueError(
            "current must immediately follow previous."
        )

    tolerance = float(tolerance)
    if not np.isfinite(tolerance):
        raise ValueError("tolerance must be finite.")
    if tolerance <= 0.0:
        raise ValueError("tolerance must be positive.")

    displacement_range_change = relative_change(
        current=current.displacement_range,
        previous=previous.displacement_range,
    )
    stress_range_change = relative_change(
        current=current.critical_stress_range,
        previous=previous.critical_stress_range,
    )
    external_work_change = relative_change(
        current=current.external_work_magnitude,
        previous=previous.external_work_magnitude,
    )
    normalized_residual_drift = normalized_drift(
        current=current.residual_displacement,
        previous=previous.residual_displacement,
        cycle_scale=current.displacement_range,
    )
    normalized_plastic_strain_drift = normalized_drift(
        current=current.critical_plastic_strain_at_end,
        previous=previous.critical_plastic_strain_at_end,
        cycle_scale=current.critical_plastic_strain_range,
    )

    maximum_indicator = max(
        displacement_range_change,
        stress_range_change,
        external_work_change,
        normalized_residual_drift,
        normalized_plastic_strain_drift,
    )

    return ConsecutiveCycleSimilarity(
        previous_cycle=previous.cycle_number,
        current_cycle=current.cycle_number,
        displacement_range_change=(
            displacement_range_change
        ),
        stress_range_change=stress_range_change,
        external_work_change=external_work_change,
        normalized_residual_drift=(
            normalized_residual_drift
        ),
        normalized_plastic_strain_drift=(
            normalized_plastic_strain_drift
        ),
        maximum_indicator=maximum_indicator,
        tolerance=tolerance,
        within_tolerance=bool(
            maximum_indicator <= tolerance
        ),
    )


def evaluate_cycle_similarities(
    diagnostics: MulticycleDiagnostics,
    tolerance: float = 1.0e-3,
) -> Tuple[ConsecutiveCycleSimilarity, ...]:
    """Compare every adjacent cycle pair."""
    if not isinstance(
        diagnostics,
        MulticycleDiagnostics,
    ):
        raise TypeError(
            "diagnostics must be a MulticycleDiagnostics object."
        )

    return tuple(
        compare_consecutive_cycles(
            previous=diagnostics.cycles[index - 1],
            current=diagnostics.cycles[index],
            tolerance=tolerance,
        )
        for index in range(1, diagnostics.n_cycles)
    )


def run_multicycle_tower_analysis(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    force_amplitude: float = 1.0e6,
    period: float = 10.0,
    n_cycles: int = 5,
    increments_per_cycle: int = 40,
    similarity_tolerance: float = 1.0e-3,
    max_iterations: int = 40,
) -> MulticycleTowerResult:
    """Run the fully reversed tower analysis and extract cycle metrics."""
    loading = create_reversed_top_force_history(
        force_amplitude=force_amplitude,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
    )
    response = run_nonlinear_reversed_analysis(
        configuration=configuration,
        material=material,
        loading=loading,
        max_iterations=max_iterations,
    )
    diagnostics = extract_multicycle_diagnostics(
        response=response,
    )
    similarities = evaluate_cycle_similarities(
        diagnostics=diagnostics,
        tolerance=similarity_tolerance,
    )

    return MulticycleTowerResult(
        response=response,
        diagnostics=diagnostics,
        similarities=similarities,
        similarity_tolerance=similarity_tolerance,
    )


def print_cycle_table(
    result: MulticycleTowerResult,
) -> None:
    """Print primary response and internal-variable metrics by cycle."""
    response = result.response
    diagnostics = result.diagnostics
    loading = response.loading
    element_index, gauss_index, fiber_index = (
        response.critical_location
    )

    print("=" * 174)
    print(
        "NREL 5 MW tower: multi-cycle fully reversed "
        "nonlinear response"
    )
    print("=" * 174)
    print(
        f"F_a={loading.force_amplitude / 1.0e6:.6g} MN, "
        f"period={loading.period:.6g}, "
        f"cycles={loading.n_cycles}, "
        f"increments/cycle={loading.increments_per_cycle}"
    )
    print(
        "Fixed critical fiber: "
        f"element {element_index + 1}, "
        f"Gauss point {gauss_index + 1}, "
        f"fiber {fiber_index + 1}, "
        f"z={response.critical_height:.6e} m, "
        f"y={response.critical_y_coordinate:.6e} m"
    )
    print("-" * 174)
    print(
        "cycle   ux(+Fa)      ux(-Fa)      ux,res       "
        "Delta ux      Delta sigma    eps_p,end    beta,end    "
        "r_bar,end    D_max,end    Delta D_max    |int Fdu|(J)  "
        "Newton   residual(N)"
    )

    for cycle in diagnostics.cycles:
        print(
            f"{cycle.cycle_number:5d}  "
            f"{cycle.displacement_at_positive_peak:11.4e}  "
            f"{cycle.displacement_at_negative_peak:11.4e}  "
            f"{cycle.residual_displacement:11.4e}  "
            f"{cycle.displacement_range:11.4e}  "
            f"{cycle.critical_stress_range:11.4e}  "
            f"{cycle.critical_plastic_strain_at_end:11.4e}  "
            f"{cycle.critical_backstress_at_end:10.4e}  "
            f"{cycle.critical_r_bar_at_end:10.4e}  "
            f"{cycle.maximum_damage_at_end:10.4e}  "
            f"{cycle.maximum_damage_increment:11.4e}  "
            f"{cycle.external_work_magnitude:12.4e}  "
            f"{cycle.maximum_newton_iterations:6d}  "
            f"{cycle.maximum_residual_norm:11.4e}"
        )

    print("-" * 174)
    print(
        "Cycle-to-cycle numerical similarity indicators"
    )
    print(
        "Each row compares cycle n with cycle n-1. "
        "All indicators are dimensionless."
    )
    print(
        "cycle pair   Delta(Delta ux)   Delta(Delta sigma)   "
        "Delta(|W|)     residual drift   eps_p drift    "
        "maximum       within tolerance"
    )

    if len(result.similarities) == 0:
        print(
            "Only one cycle is available; no consecutive-cycle "
            "comparison can be made."
        )
    else:
        for similarity in result.similarities:
            print(
                f"{similarity.previous_cycle:3d}->{similarity.current_cycle:<3d}  "
                f"{similarity.displacement_range_change:14.6e}  "
                f"{similarity.stress_range_change:18.6e}  "
                f"{similarity.external_work_change:12.6e}  "
                f"{similarity.normalized_residual_drift:14.6e}  "
                f"{similarity.normalized_plastic_strain_drift:12.6e}  "
                f"{similarity.maximum_indicator:12.6e}  "
                f"{str(similarity.within_tolerance):>16s}"
            )

    print("-" * 174)
    print(
        "Similarity tolerance: "
        f"{result.similarity_tolerance:.6e}"
    )
    print(
        "Final cycle pair within tolerance: "
        f"{result.final_cycle_is_within_tolerance}"
    )
    print(
        "Interpretation: this flag only indicates numerical "
        "cycle-to-cycle similarity under the selected metrics and "
        "tolerance; it is not by itself proof of a stabilized "
        "material limit cycle."
    )
    print("=" * 174)


def plot_results(
    result: MulticycleTowerResult,
) -> None:
    """Plot cycle loops and cycle-by-cycle evolution."""
    response = result.response
    diagnostics = result.diagnostics
    loading = response.loading
    increments = loading.increments_per_cycle

    plt.figure()
    for cycle_number in range(
        1,
        loading.n_cycles + 1,
    ):
        start = (
            cycle_number - 1
        ) * increments
        end = cycle_number * increments + 1
        plt.plot(
            response.top_displacements[start:end],
            loading.forces[start:end] / 1.0e6,
            label=f"Cycle {cycle_number}",
        )
    plt.xlabel("Tower-top horizontal displacement (m)")
    plt.ylabel("Tower-top horizontal force (MN)")
    plt.title("Multi-cycle tower force-displacement response")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    for cycle_number in range(
        1,
        loading.n_cycles + 1,
    ):
        start = (
            cycle_number - 1
        ) * increments
        end = cycle_number * increments + 1
        plt.plot(
            response.critical_fiber_strains[start:end],
            response.critical_fiber_stresses[start:end],
            label=f"Cycle {cycle_number}",
        )
    plt.xlabel("Fixed critical-fiber total strain")
    plt.ylabel("Fixed critical-fiber stress (MPa)")
    plt.title("Multi-cycle critical-fiber stress-strain response")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    cycle_numbers = diagnostics.cycle_numbers

    plt.figure()
    plt.plot(
        cycle_numbers,
        diagnostics.residual_displacements,
        marker="o",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Residual tower-top displacement (m)")
    plt.title("Residual displacement evolution")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        cycle_numbers,
        diagnostics.maximum_damage_ends,
        marker="o",
        label="Maximum damage at cycle end",
    )
    plt.plot(
        cycle_numbers,
        diagnostics.maximum_damage_increments,
        marker="o",
        label="Maximum-damage increment",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Damage")
    plt.title("Cycle-by-cycle damage evolution")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(
        cycle_numbers,
        diagnostics.external_work_magnitudes,
        marker="o",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("External-work magnitude, |int F du| (J)")
    plt.title("Cycle-by-cycle force-displacement work")
    plt.grid(True)
    plt.tight_layout()

    if len(result.similarities) > 0:
        comparison_cycles = np.asarray(
            [
                similarity.current_cycle
                for similarity in result.similarities
            ],
            dtype=np.int64,
        )
        maximum_indicators = np.asarray(
            [
                similarity.maximum_indicator
                for similarity in result.similarities
            ],
            dtype=np.float64,
        )

        plt.figure()
        plt.semilogy(
            comparison_cycles,
            maximum_indicators,
            marker="o",
            label="Maximum similarity indicator",
        )
        plt.axhline(
            result.similarity_tolerance,
            linestyle="--",
            label="Selected tolerance",
        )
        plt.xlabel("Current cycle number")
        plt.ylabel("Dimensionless indicator")
        plt.title("Consecutive-cycle similarity")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Run and diagnose multiple fully reversed nonlinear "
            "cycles of the NREL 5 MW tower."
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
        help="Duration of one complete reversed cycle.",
    )
    parser.add_argument(
        "--cycles",
        type=int,
        default=5,
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
        "--similarity-tolerance",
        type=float,
        default=1.0e-3,
        help=(
            "Dimensionless threshold for the maximum "
            "consecutive-cycle similarity indicator."
        ),
    )
    parser.add_argument(
        "--max-iterations",
        type=int,
        default=40,
        help="Maximum Newton iterations per load step.",
    )
    parser.add_argument(
        "--no-plot",
        action="store_true",
        help="Run without opening figures.",
    )
    return parser.parse_args()


def main() -> None:
    """Run the multi-cycle nonlinear tower example."""
    arguments = parse_arguments()

    configuration = TowerConfiguration(
        horizontal_force=arguments.force_amplitude,
        n_elements=arguments.elements,
        n_gauss=arguments.gauss,
        n_circumferential=arguments.circumferential,
        n_radial=arguments.radial,
    )
    material = MaterialParameters()

    result = run_multicycle_tower_analysis(
        configuration=configuration,
        material=material,
        force_amplitude=arguments.force_amplitude,
        period=arguments.period,
        n_cycles=arguments.cycles,
        increments_per_cycle=arguments.increments,
        similarity_tolerance=(
            arguments.similarity_tolerance
        ),
        max_iterations=arguments.max_iterations,
    )
    print_cycle_table(result)

    if not arguments.no_plot:
        plot_results(result)


if __name__ == "__main__":
    main()

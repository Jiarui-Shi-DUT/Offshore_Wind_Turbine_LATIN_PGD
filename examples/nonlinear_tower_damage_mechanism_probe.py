# -*- coding: utf-8 -*-
"""
Mechanism-separation probe for nonlinear multi-cycle tower response.

Two analyses are run under exactly the same tower discretisation, loading,
plasticity, hardening, and time integration:

    1. coupled viscoplastic-damage response;
    2. viscoplastic response with damage evolution disabled by k_damage = 0.

The comparison isolates the incremental influence of damage coupling on the
global and local cyclic response. The same physical fiber location is used in
both analyses: the critical location selected from the coupled response is
imposed on the damage-disabled response before cycle diagnostics are
extracted.

This is a numerical mechanism-separation study. It does not represent a new
material calibration, because all non-damage parameters are intentionally kept
unchanged.

Use:
    python -m examples.nonlinear_tower_damage_mechanism_probe --no-plot
    python -m examples.nonlinear_tower_damage_mechanism_probe
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass, replace
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_multicycle_diagnostics import (
    MulticycleDiagnostics,
    extract_multicycle_diagnostics,
)
from examples.nonlinear_tower_multicycle_response import (
    MulticycleTowerResult,
    evaluate_cycle_similarities,
    run_multicycle_tower_analysis,
)
from examples.nonlinear_tower_reversed_response import (
    NonlinearReversedResponse,
)
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


def _finite_scalar(value: float, name: str) -> float:
    """Return a validated finite scalar."""
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def relative_difference(
    coupled: float,
    reference: float,
    scale_floor: float = 1.0e-14,
) -> float:
    """
    Return (coupled - reference) / |reference| with protection near zero.

    Positive values mean that damage coupling increases the quantity relative
    to the damage-disabled reference. Negative values mean that it reduces the
    quantity.
    """
    coupled = _finite_scalar(coupled, "coupled")
    reference = _finite_scalar(reference, "reference")
    scale_floor = _finite_scalar(
        scale_floor,
        "scale_floor",
    )
    if scale_floor <= 0.0:
        raise ValueError("scale_floor must be positive.")

    denominator = max(abs(reference), scale_floor)
    return float(
        (coupled - reference) / denominator
    )


@dataclass(frozen=True)
class DamageMechanismCycleComparison:
    """Cycle-level coupled versus damage-disabled diagnostics."""

    cycle_number: int

    coupled_displacement_range: float
    undamaged_displacement_range: float
    displacement_range_relative_difference: float

    coupled_stress_range: float
    undamaged_stress_range: float
    stress_range_relative_difference: float

    coupled_residual_displacement: float
    undamaged_residual_displacement: float
    residual_displacement_difference: float

    coupled_plastic_strain_end: float
    undamaged_plastic_strain_end: float
    plastic_strain_end_difference: float

    coupled_external_work: float
    undamaged_external_work: float
    external_work_relative_difference: float

    coupled_maximum_damage_end: float
    coupled_maximum_damage_increment: float
    undamaged_maximum_damage_end: float

    def __post_init__(self) -> None:
        if isinstance(self.cycle_number, bool):
            raise TypeError(
                "cycle_number must be an integer."
            )
        cycle_number = int(self.cycle_number)
        if cycle_number != self.cycle_number:
            raise TypeError(
                "cycle_number must be an integer."
            )
        if cycle_number < 1:
            raise ValueError(
                "cycle_number must be at least one."
            )

        scalar_names = (
            "coupled_displacement_range",
            "undamaged_displacement_range",
            "displacement_range_relative_difference",
            "coupled_stress_range",
            "undamaged_stress_range",
            "stress_range_relative_difference",
            "coupled_residual_displacement",
            "undamaged_residual_displacement",
            "residual_displacement_difference",
            "coupled_plastic_strain_end",
            "undamaged_plastic_strain_end",
            "plastic_strain_end_difference",
            "coupled_external_work",
            "undamaged_external_work",
            "external_work_relative_difference",
            "coupled_maximum_damage_end",
            "coupled_maximum_damage_increment",
            "undamaged_maximum_damage_end",
        )
        for name in scalar_names:
            value = _finite_scalar(
                getattr(self, name),
                name,
            )
            object.__setattr__(self, name, value)

        nonnegative_names = (
            "coupled_displacement_range",
            "undamaged_displacement_range",
            "coupled_stress_range",
            "undamaged_stress_range",
            "coupled_external_work",
            "undamaged_external_work",
            "coupled_maximum_damage_end",
            "coupled_maximum_damage_increment",
            "undamaged_maximum_damage_end",
        )
        for name in nonnegative_names:
            if getattr(self, name) < 0.0:
                raise ValueError(
                    name + " must be non-negative."
                )

        expected_displacement_difference = (
            self.coupled_residual_displacement
            - self.undamaged_residual_displacement
        )
        if not np.isclose(
            self.residual_displacement_difference,
            expected_displacement_difference,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "residual_displacement_difference is inconsistent."
            )

        expected_plastic_difference = (
            self.coupled_plastic_strain_end
            - self.undamaged_plastic_strain_end
        )
        if not np.isclose(
            self.plastic_strain_end_difference,
            expected_plastic_difference,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "plastic_strain_end_difference is inconsistent."
            )

        expected_displacement_relative = (
            relative_difference(
                self.coupled_displacement_range,
                self.undamaged_displacement_range,
            )
        )
        if not np.isclose(
            self.displacement_range_relative_difference,
            expected_displacement_relative,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "displacement_range_relative_difference "
                "is inconsistent."
            )

        expected_stress_relative = relative_difference(
            self.coupled_stress_range,
            self.undamaged_stress_range,
        )
        if not np.isclose(
            self.stress_range_relative_difference,
            expected_stress_relative,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "stress_range_relative_difference "
                "is inconsistent."
            )

        expected_work_relative = relative_difference(
            self.coupled_external_work,
            self.undamaged_external_work,
        )
        if not np.isclose(
            self.external_work_relative_difference,
            expected_work_relative,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "external_work_relative_difference "
                "is inconsistent."
            )

        object.__setattr__(
            self,
            "cycle_number",
            cycle_number,
        )


@dataclass(frozen=True)
class DamageMechanismComparison:
    """Complete paired analyses and cycle-level comparisons."""

    coupled: MulticycleTowerResult
    damage_disabled: MulticycleTowerResult
    cycles: Tuple[DamageMechanismCycleComparison, ...]

    def __post_init__(self) -> None:
        if not isinstance(
            self.coupled,
            MulticycleTowerResult,
        ):
            raise TypeError(
                "coupled must be a MulticycleTowerResult."
            )
        if not isinstance(
            self.damage_disabled,
            MulticycleTowerResult,
        ):
            raise TypeError(
                "damage_disabled must be a "
                "MulticycleTowerResult."
            )

        coupled_loading = self.coupled.response.loading
        disabled_loading = (
            self.damage_disabled.response.loading
        )
        loading_scalars = (
            (
                coupled_loading.force_amplitude,
                disabled_loading.force_amplitude,
            ),
            (
                coupled_loading.period,
                disabled_loading.period,
            ),
            (
                coupled_loading.n_cycles,
                disabled_loading.n_cycles,
            ),
            (
                coupled_loading.increments_per_cycle,
                disabled_loading.increments_per_cycle,
            ),
        )
        for coupled_value, disabled_value in loading_scalars:
            if coupled_value != disabled_value:
                raise ValueError(
                    "Paired analyses must use identical loading."
                )

        if (
            self.coupled.response.critical_location
            != self.damage_disabled.response.critical_location
        ):
            raise ValueError(
                "Paired diagnostics must use the same "
                "physical fiber location."
            )

        if self.damage_disabled.response.material.k_damage != 0.0:
            raise ValueError(
                "damage_disabled must use k_damage = 0."
            )

        cycles = tuple(self.cycles)
        if len(cycles) != coupled_loading.n_cycles:
            raise ValueError(
                "cycles must compare every loading cycle."
            )
        for expected_number, cycle in enumerate(
            cycles,
            start=1,
        ):
            if not isinstance(
                cycle,
                DamageMechanismCycleComparison,
            ):
                raise TypeError(
                    "Every item must be a "
                    "DamageMechanismCycleComparison object."
                )
            if cycle.cycle_number != expected_number:
                raise ValueError(
                    "Cycle comparisons must be ordered."
                )

        object.__setattr__(self, "cycles", cycles)

    @property
    def n_cycles(self) -> int:
        """Return the number of paired cycles."""
        return len(self.cycles)

    def _history(self, attribute: str) -> FloatArray:
        """Return one comparison quantity across all cycles."""
        return np.asarray(
            [
                float(getattr(cycle, attribute))
                for cycle in self.cycles
            ],
            dtype=np.float64,
        )

    @property
    def cycle_numbers(self) -> NDArray[np.int64]:
        """Return one-based cycle numbers."""
        return np.asarray(
            [
                cycle.cycle_number
                for cycle in self.cycles
            ],
            dtype=np.int64,
        )

    @property
    def displacement_range_relative_differences(
        self,
    ) -> FloatArray:
        """Return damage-induced displacement-range changes."""
        return self._history(
            "displacement_range_relative_difference"
        )

    @property
    def stress_range_relative_differences(
        self,
    ) -> FloatArray:
        """Return damage-induced stress-range changes."""
        return self._history(
            "stress_range_relative_difference"
        )

    @property
    def external_work_relative_differences(
        self,
    ) -> FloatArray:
        """Return damage-induced external-work changes."""
        return self._history(
            "external_work_relative_difference"
        )

    @property
    def coupled_damage_ends(self) -> FloatArray:
        """Return coupled maximum damage at cycle ends."""
        return self._history(
            "coupled_maximum_damage_end"
        )

    @property
    def disabled_damage_ends(self) -> FloatArray:
        """Return disabled-case maximum damage at cycle ends."""
        return self._history(
            "undamaged_maximum_damage_end"
        )


def response_with_common_critical_location(
    response: NonlinearReversedResponse,
    reference: NonlinearReversedResponse,
) -> NonlinearReversedResponse:
    """Return response diagnostics anchored to reference's fiber location."""
    if not isinstance(
        response,
        NonlinearReversedResponse,
    ):
        raise TypeError(
            "response must be a NonlinearReversedResponse."
        )
    if not isinstance(
        reference,
        NonlinearReversedResponse,
    ):
        raise TypeError(
            "reference must be a NonlinearReversedResponse."
        )

    if response.fiber_states.shape[1:] != reference.fiber_states.shape[1:]:
        raise ValueError(
            "Responses must use matching fiber discretisations."
        )

    return replace(
        response,
        critical_location=reference.critical_location,
        critical_height=reference.critical_height,
        critical_y_coordinate=(
            reference.critical_y_coordinate
        ),
    )


def build_multicycle_result(
    response: NonlinearReversedResponse,
    similarity_tolerance: float,
) -> MulticycleTowerResult:
    """Build diagnostics and similarity comparisons for one response."""
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


def compare_cycle_diagnostics(
    coupled: MulticycleDiagnostics,
    damage_disabled: MulticycleDiagnostics,
) -> Tuple[DamageMechanismCycleComparison, ...]:
    """Build paired cycle-level mechanism comparisons."""
    if not isinstance(
        coupled,
        MulticycleDiagnostics,
    ):
        raise TypeError(
            "coupled must be a MulticycleDiagnostics object."
        )
    if not isinstance(
        damage_disabled,
        MulticycleDiagnostics,
    ):
        raise TypeError(
            "damage_disabled must be a "
            "MulticycleDiagnostics object."
        )
    if coupled.n_cycles != damage_disabled.n_cycles:
        raise ValueError(
            "Paired diagnostics must contain the same "
            "number of cycles."
        )

    comparisons = []
    for coupled_cycle, disabled_cycle in zip(
        coupled.cycles,
        damage_disabled.cycles,
    ):
        if (
            coupled_cycle.cycle_number
            != disabled_cycle.cycle_number
        ):
            raise ValueError(
                "Paired cycles must have matching numbers."
            )

        comparisons.append(
            DamageMechanismCycleComparison(
                cycle_number=coupled_cycle.cycle_number,
                coupled_displacement_range=(
                    coupled_cycle.displacement_range
                ),
                undamaged_displacement_range=(
                    disabled_cycle.displacement_range
                ),
                displacement_range_relative_difference=(
                    relative_difference(
                        coupled_cycle.displacement_range,
                        disabled_cycle.displacement_range,
                    )
                ),
                coupled_stress_range=(
                    coupled_cycle.critical_stress_range
                ),
                undamaged_stress_range=(
                    disabled_cycle.critical_stress_range
                ),
                stress_range_relative_difference=(
                    relative_difference(
                        coupled_cycle.critical_stress_range,
                        disabled_cycle.critical_stress_range,
                    )
                ),
                coupled_residual_displacement=(
                    coupled_cycle.residual_displacement
                ),
                undamaged_residual_displacement=(
                    disabled_cycle.residual_displacement
                ),
                residual_displacement_difference=(
                    coupled_cycle.residual_displacement
                    - disabled_cycle.residual_displacement
                ),
                coupled_plastic_strain_end=(
                    coupled_cycle
                    .critical_plastic_strain_at_end
                ),
                undamaged_plastic_strain_end=(
                    disabled_cycle
                    .critical_plastic_strain_at_end
                ),
                plastic_strain_end_difference=(
                    coupled_cycle
                    .critical_plastic_strain_at_end
                    - disabled_cycle
                    .critical_plastic_strain_at_end
                ),
                coupled_external_work=(
                    coupled_cycle.external_work_magnitude
                ),
                undamaged_external_work=(
                    disabled_cycle.external_work_magnitude
                ),
                external_work_relative_difference=(
                    relative_difference(
                        coupled_cycle
                        .external_work_magnitude,
                        disabled_cycle
                        .external_work_magnitude,
                    )
                ),
                coupled_maximum_damage_end=(
                    coupled_cycle.maximum_damage_at_end
                ),
                coupled_maximum_damage_increment=(
                    coupled_cycle.maximum_damage_increment
                ),
                undamaged_maximum_damage_end=(
                    disabled_cycle.maximum_damage_at_end
                ),
            )
        )

    return tuple(comparisons)


def run_damage_mechanism_comparison(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    force_amplitude: float = 1.0e6,
    period: float = 10.0,
    n_cycles: int = 5,
    increments_per_cycle: int = 40,
    similarity_tolerance: float = 1.0e-3,
    max_iterations: int = 40,
) -> DamageMechanismComparison:
    """Run paired coupled and damage-disabled multi-cycle analyses."""
    if not isinstance(
        configuration,
        TowerConfiguration,
    ):
        raise TypeError(
            "configuration must be a TowerConfiguration."
        )
    if not isinstance(
        material,
        MaterialParameters,
    ):
        raise TypeError(
            "material must be a MaterialParameters object."
        )

    coupled = run_multicycle_tower_analysis(
        configuration=configuration,
        material=material,
        force_amplitude=force_amplitude,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
        similarity_tolerance=similarity_tolerance,
        max_iterations=max_iterations,
    )

    disabled_material = replace(
        material,
        k_damage=0.0,
    )
    raw_disabled = run_multicycle_tower_analysis(
        configuration=configuration,
        material=disabled_material,
        force_amplitude=force_amplitude,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
        similarity_tolerance=similarity_tolerance,
        max_iterations=max_iterations,
    )

    common_disabled_response = (
        response_with_common_critical_location(
            response=raw_disabled.response,
            reference=coupled.response,
        )
    )
    damage_disabled = build_multicycle_result(
        response=common_disabled_response,
        similarity_tolerance=similarity_tolerance,
    )

    comparisons = compare_cycle_diagnostics(
        coupled=coupled.diagnostics,
        damage_disabled=damage_disabled.diagnostics,
    )

    return DamageMechanismComparison(
        coupled=coupled,
        damage_disabled=damage_disabled,
        cycles=comparisons,
    )


def print_comparison(
    comparison: DamageMechanismComparison,
) -> None:
    """Print cycle-by-cycle damage-mechanism separation results."""
    coupled_response = comparison.coupled.response
    loading = coupled_response.loading
    element_index, gauss_index, fiber_index = (
        coupled_response.critical_location
    )

    print("=" * 188)
    print(
        "NREL 5 MW tower: damage-coupling mechanism separation"
    )
    print("=" * 188)
    print(
        f"F_a={loading.force_amplitude / 1.0e6:.6g} MN, "
        f"period={loading.period:.6g}, "
        f"cycles={loading.n_cycles}, "
        f"increments/cycle={loading.increments_per_cycle}"
    )
    print(
        "Common diagnostic fiber: "
        f"element {element_index + 1}, "
        f"Gauss point {gauss_index + 1}, "
        f"fiber {fiber_index + 1}, "
        f"z={coupled_response.critical_height:.6e} m, "
        f"y={coupled_response.critical_y_coordinate:.6e} m"
    )
    print(
        "Cases: coupled k_damage="
        f"{coupled_response.material.k_damage:.6g}; "
        "damage-disabled k_damage=0"
    )
    print("-" * 188)
    print(
        "cycle   Delta ux,c   Delta ux,0   rel.diff ux    "
        "Delta sig,c   Delta sig,0   rel.diff sig   "
        "|W|,c(J)     |W|,0(J)     rel.diff W    "
        "D_c,end     Delta D_c    D_0,end"
    )

    for cycle in comparison.cycles:
        print(
            f"{cycle.cycle_number:5d}  "
            f"{cycle.coupled_displacement_range:11.4e}  "
            f"{cycle.undamaged_displacement_range:11.4e}  "
            f"{cycle.displacement_range_relative_difference:11.4e}  "
            f"{cycle.coupled_stress_range:11.4e}  "
            f"{cycle.undamaged_stress_range:11.4e}  "
            f"{cycle.stress_range_relative_difference:11.4e}  "
            f"{cycle.coupled_external_work:11.4e}  "
            f"{cycle.undamaged_external_work:11.4e}  "
            f"{cycle.external_work_relative_difference:11.4e}  "
            f"{cycle.coupled_maximum_damage_end:10.4e}  "
            f"{cycle.coupled_maximum_damage_increment:10.4e}  "
            f"{cycle.undamaged_maximum_damage_end:10.4e}"
        )

    print("-" * 188)
    print(
        "cycle   ux,res,c      ux,res,0      difference     "
        "eps_p,end,c   eps_p,end,0   difference"
    )
    for cycle in comparison.cycles:
        print(
            f"{cycle.cycle_number:5d}  "
            f"{cycle.coupled_residual_displacement:12.5e}  "
            f"{cycle.undamaged_residual_displacement:12.5e}  "
            f"{cycle.residual_displacement_difference:12.5e}  "
            f"{cycle.coupled_plastic_strain_end:12.5e}  "
            f"{cycle.undamaged_plastic_strain_end:12.5e}  "
            f"{cycle.plastic_strain_end_difference:12.5e}"
        )

    final_cycle = comparison.cycles[-1]
    print("-" * 188)
    print(
        "Final-cycle damage-induced displacement-range change: "
        f"{100.0 * final_cycle.displacement_range_relative_difference:.6f}%"
    )
    print(
        "Final-cycle damage-induced stress-range change: "
        f"{100.0 * final_cycle.stress_range_relative_difference:.6f}%"
    )
    print(
        "Final-cycle damage-induced external-work change: "
        f"{100.0 * final_cycle.external_work_relative_difference:.6f}%"
    )
    print(
        "Maximum damage in disabled case: "
        f"{float(np.max(comparison.disabled_damage_ends)):.6e}"
    )
    print(
        "Interpretation: differences between the paired cases isolate "
        "the numerical influence of damage evolution while plasticity, "
        "hardening, loading, discretisation, and time integration remain "
        "unchanged."
    )
    print("=" * 188)


def plot_comparison(
    comparison: DamageMechanismComparison,
) -> None:
    """Plot coupled and damage-disabled histories."""
    cycle_numbers = comparison.cycle_numbers
    coupled = comparison.coupled.diagnostics
    disabled = comparison.damage_disabled.diagnostics

    plt.figure()
    plt.plot(
        cycle_numbers,
        coupled.displacement_ranges,
        marker="o",
        label="Coupled damage",
    )
    plt.plot(
        cycle_numbers,
        disabled.displacement_ranges,
        marker="o",
        label="Damage disabled",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Tower-top displacement range (m)")
    plt.title("Damage influence on displacement range")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(
        cycle_numbers,
        coupled.critical_stress_ranges,
        marker="o",
        label="Coupled damage",
    )
    plt.plot(
        cycle_numbers,
        disabled.critical_stress_ranges,
        marker="o",
        label="Damage disabled",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Common-fiber stress range (MPa)")
    plt.title("Damage influence on fixed-fiber stress range")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(
        cycle_numbers,
        coupled.external_work_magnitudes,
        marker="o",
        label="Coupled damage",
    )
    plt.plot(
        cycle_numbers,
        disabled.external_work_magnitudes,
        marker="o",
        label="Damage disabled",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("External-work magnitude, |int F du| (J)")
    plt.title("Damage influence on force-displacement work")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(
        cycle_numbers,
        comparison.displacement_range_relative_differences,
        marker="o",
        label="Displacement range",
    )
    plt.plot(
        cycle_numbers,
        comparison.stress_range_relative_differences,
        marker="o",
        label="Stress range",
    )
    plt.plot(
        cycle_numbers,
        comparison.external_work_relative_differences,
        marker="o",
        label="External work",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Relative difference")
    plt.title("Cycle-by-cycle damage contribution")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.figure()
    plt.plot(
        cycle_numbers,
        comparison.coupled_damage_ends,
        marker="o",
        label="Coupled damage",
    )
    plt.plot(
        cycle_numbers,
        comparison.disabled_damage_ends,
        marker="o",
        label="Damage disabled",
    )
    plt.xlabel("Cycle number")
    plt.ylabel("Maximum damage at cycle end")
    plt.title("Damage-switch verification")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare coupled and damage-disabled multi-cycle "
            "tower responses."
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
        default=4,
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
        default=8,
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
            "Dimensionless threshold used inside each paired "
            "multi-cycle result."
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
    """Run the paired damage mechanism-separation study."""
    arguments = parse_arguments()

    configuration = TowerConfiguration(
        horizontal_force=arguments.force_amplitude,
        n_elements=arguments.elements,
        n_gauss=arguments.gauss,
        n_circumferential=arguments.circumferential,
        n_radial=arguments.radial,
    )
    material = MaterialParameters()

    comparison = run_damage_mechanism_comparison(
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
    print_comparison(comparison)

    if not arguments.no_plot:
        plot_comparison(comparison)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""
Quantitative Bauschinger-effect probe for the nonlinear NREL 5 MW tower.

The probe reuses the fully reversed tower analysis and examines one fixed
critical fiber. It distinguishes two related quantities:

1. Observed yield-onset stress:
   the stress at which the evaluated yield function crosses zero along the
   initial and reverse loading branches.

2. Kinematic yield-surface translation:
   the reduction of the reverse-direction yield limit caused by the
   backstress beta, relative to a counterfactual state with the same
   isotropic hardening and damage but beta = 0.

The second quantity isolates the Bauschinger contribution from simultaneous
isotropic hardening and damage evolution.

Use:
    python -m examples.nonlinear_tower_bauschinger_probe
    python -m examples.nonlinear_tower_bauschinger_probe --no-plot
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_reversed_response import (
    NonlinearReversedResponse,
    run_nonlinear_reversed_analysis,
)
from fem.tower_loading import create_reversed_top_force_history
from material.viscoplastic_damage_1d import (
    MaterialParameters,
    isotropic_hardening_force,
)


@dataclass(frozen=True)
class YieldOnset:
    """Interpolated yield-function zero crossing on one loading branch."""

    lower_index: int
    upper_index: int
    normalized_time: float
    force: float
    stress: float
    plastic_strain: float
    backstress: float
    yield_function: float


@dataclass(frozen=True)
class BauschingerMetrics:
    """Observed and constitutive measures of reverse-yield translation."""

    initial_onset: YieldOnset
    reverse_onset: YieldOnset
    virgin_yield_stress: float
    reversal_backstress: float
    reversal_isotropic_force: float
    reversal_damage: float
    reverse_limit_with_kinematic: float
    reverse_limit_without_kinematic: float

    @property
    def observed_initial_yield_magnitude(self) -> float:
        """Return absolute initial-branch yield-onset stress."""
        return abs(float(self.initial_onset.stress))

    @property
    def observed_reverse_yield_magnitude(self) -> float:
        """Return absolute reverse-branch yield-onset stress."""
        return abs(float(self.reverse_onset.stress))

    @property
    def observed_change_from_initial(self) -> float:
        """
        Return initial minus reverse observed yield magnitude.

        A positive value means that the reverse onset occurs at a lower
        stress magnitude than the initial onset. A negative value can occur
        when isotropic hardening exceeds the kinematic translation.
        """
        return float(
            self.observed_initial_yield_magnitude
            - self.observed_reverse_yield_magnitude
        )

    @property
    def kinematic_reverse_yield_reduction(self) -> float:
        """Return the reverse-limit reduction caused by backstress."""
        return float(
            self.reverse_limit_without_kinematic
            - self.reverse_limit_with_kinematic
        )

    @property
    def kinematic_reduction_ratio(self) -> float:
        """Return the kinematic reverse-limit reduction as a fraction."""
        denominator = abs(
            float(self.reverse_limit_without_kinematic)
        )
        if denominator <= np.finfo(np.float64).eps:
            return 0.0
        return float(
            self.kinematic_reverse_yield_reduction
            / denominator
        )


def interpolate_zero_crossing(
    response: NonlinearReversedResponse,
    lower_index: int,
    upper_index: int,
) -> YieldOnset:
    """Interpolate one yield-function zero crossing."""
    if lower_index < 0:
        raise ValueError("lower_index must be non-negative.")
    if upper_index != lower_index + 1:
        raise ValueError(
            "upper_index must immediately follow lower_index."
        )

    yield_functions = (
        response.critical_fiber_yield_functions
    )
    lower_value = float(yield_functions[lower_index])
    upper_value = float(yield_functions[upper_index])

    if lower_value > 0.0 or upper_value <= 0.0:
        raise ValueError(
            "Indices must bracket a negative-to-positive crossing."
        )

    denominator = upper_value - lower_value
    if denominator <= 0.0:
        raise ValueError(
            "Yield-function crossing must have positive slope."
        )

    fraction = float(-lower_value / denominator)
    fraction = float(np.clip(fraction, 0.0, 1.0))

    def interpolate(values: np.ndarray) -> float:
        lower = float(values[lower_index])
        upper = float(values[upper_index])
        return float(lower + fraction * (upper - lower))

    normalized_times = (
        response.loading.times / response.loading.period
    )

    return YieldOnset(
        lower_index=int(lower_index),
        upper_index=int(upper_index),
        normalized_time=interpolate(normalized_times),
        force=interpolate(response.loading.forces),
        stress=interpolate(
            response.critical_fiber_stresses
        ),
        plastic_strain=interpolate(
            response.critical_fiber_plastic_strains
        ),
        backstress=interpolate(
            response.critical_fiber_backstresses
        ),
        yield_function=0.0,
    )


def find_yield_onset(
    response: NonlinearReversedResponse,
    start_index: int,
    end_index: int,
) -> YieldOnset:
    """Find the first negative-to-positive yield crossing in a range."""
    yield_functions = (
        response.critical_fiber_yield_functions
    )

    if start_index < 0:
        raise ValueError("start_index must be non-negative.")
    if end_index >= yield_functions.size:
        raise ValueError("end_index exceeds the history length.")
    if end_index <= start_index:
        raise ValueError(
            "end_index must be greater than start_index."
        )

    for upper_index in range(start_index + 1, end_index + 1):
        lower_index = upper_index - 1
        if (
            float(yield_functions[lower_index]) <= 0.0
            and float(yield_functions[upper_index]) > 0.0
        ):
            return interpolate_zero_crossing(
                response=response,
                lower_index=lower_index,
                upper_index=upper_index,
            )

    raise RuntimeError(
        "No negative-to-positive yield-function crossing "
        "was found in the requested branch."
    )


def directional_reverse_yield_limits(
    material: MaterialParameters,
    state: np.ndarray,
) -> Tuple[float, float, float, float, float]:
    """
    Return positive reverse-yield limits at one unloading state.

    Returns:
        beta,
        isotropic force R,
        damage,
        limit with kinematic hardening,
        counterfactual limit with beta = 0.
    """
    values = np.asarray(state, dtype=np.float64)
    if values.shape != (4,):
        raise ValueError("state must have shape (4,).")
    if np.any(~np.isfinite(values)):
        raise ValueError("state must contain finite values.")

    alpha = float(values[1])
    r_bar = float(values[2])
    damage = float(values[3])
    beta = float(material.C * alpha)
    _r_bar_force, isotropic_force = (
        isotropic_hardening_force(
            r_bar=r_bar,
            material=material,
        )
    )

    nonlinear_kinematic_term = float(
        material.a * beta**2 / (2.0 * material.C)
    )
    undamaged_directional_limit = float(
        material.sigma_y
        + isotropic_force
        - nonlinear_kinematic_term
        + beta
    )
    limit_with_kinematic = float(
        (1.0 - damage) * undamaged_directional_limit
    )
    limit_without_kinematic = float(
        (1.0 - damage)
        * (material.sigma_y + isotropic_force)
    )

    return (
        beta,
        float(isotropic_force),
        damage,
        limit_with_kinematic,
        limit_without_kinematic,
    )


def evaluate_bauschinger_metrics(
    response: NonlinearReversedResponse,
) -> BauschingerMetrics:
    """Evaluate observed and constitutive reverse-yield metrics."""
    loading = response.loading
    if loading.n_cycles < 1:
        raise ValueError("At least one cycle is required.")

    quarter = loading.increments_per_cycle // 4
    half = 2 * quarter
    three_quarters = 3 * quarter

    initial_onset = find_yield_onset(
        response=response,
        start_index=0,
        end_index=quarter,
    )
    reverse_onset = find_yield_onset(
        response=response,
        start_index=half,
        end_index=three_quarters,
    )

    reversal_state = response.critical_fiber_states[half]
    (
        beta,
        isotropic_force,
        damage,
        limit_with_kinematic,
        limit_without_kinematic,
    ) = directional_reverse_yield_limits(
        material=response.material,
        state=reversal_state,
    )

    return BauschingerMetrics(
        initial_onset=initial_onset,
        reverse_onset=reverse_onset,
        virgin_yield_stress=float(
            response.material.sigma_y
        ),
        reversal_backstress=beta,
        reversal_isotropic_force=isotropic_force,
        reversal_damage=damage,
        reverse_limit_with_kinematic=(
            limit_with_kinematic
        ),
        reverse_limit_without_kinematic=(
            limit_without_kinematic
        ),
    )


def print_yield_onset(
    label: str,
    onset: YieldOnset,
) -> None:
    """Print one interpolated yield-onset state."""
    print(
        f"{label:<18}"
        f"{onset.normalized_time:12.6f}"
        f"{onset.force / 1.0e6:14.6f}"
        f"{onset.stress:16.6f}"
        f"{onset.plastic_strain:16.6e}"
        f"{onset.backstress:14.6f}"
    )


def print_summary(
    response: NonlinearReversedResponse,
    metrics: BauschingerMetrics,
) -> None:
    """Print the quantitative Bauschinger-effect diagnostics."""
    element_index, gauss_index, fiber_index = (
        response.critical_location
    )

    print("=" * 106)
    print(
        "Quantitative Bauschinger-effect probe: "
        "NREL 5 MW tower"
    )
    print("=" * 106)
    print(
        f"Force amplitude="
        f"{response.loading.force_amplitude / 1.0e6:.6g} MN, "
        f"increments/cycle="
        f"{response.loading.increments_per_cycle}"
    )
    print(
        "Fixed critical fiber: "
        f"element {element_index + 1}, "
        f"Gauss point {gauss_index + 1}, "
        f"fiber {fiber_index + 1}, "
        f"height={response.critical_height:.6e} m, "
        f"y={response.critical_y_coordinate:.6e} m"
    )
    print("-" * 106)
    print(
        f"{'onset':<18}"
        f"{'t/T':>12}"
        f"{'force(MN)':>14}"
        f"{'stress(MPa)':>16}"
        f"{'eps_p':>16}"
        f"{'beta(MPa)':>14}"
    )
    print_yield_onset(
        "initial loading",
        metrics.initial_onset,
    )
    print_yield_onset(
        "reverse loading",
        metrics.reverse_onset,
    )
    print("-" * 106)
    print(
        "Observed initial yield magnitude: "
        f"{metrics.observed_initial_yield_magnitude:.6f} MPa"
    )
    print(
        "Observed reverse yield magnitude: "
        f"{metrics.observed_reverse_yield_magnitude:.6f} MPa"
    )
    print(
        "Observed initial-minus-reverse change: "
        f"{metrics.observed_change_from_initial:.6f} MPa"
    )
    print(
        "State at first zero-force reversal point: "
        f"beta={metrics.reversal_backstress:.6f} MPa, "
        f"R={metrics.reversal_isotropic_force:.6f} MPa, "
        f"D={metrics.reversal_damage:.6e}"
    )
    print(
        "Predicted positive reverse-yield limit "
        "with kinematic hardening: "
        f"{metrics.reverse_limit_with_kinematic:.6f} MPa"
    )
    print(
        "Counterfactual reverse-yield limit "
        "with the same R and D but beta=0: "
        f"{metrics.reverse_limit_without_kinematic:.6f} MPa"
    )
    print(
        "Kinematic reverse-yield reduction: "
        f"{metrics.kinematic_reverse_yield_reduction:.6f} MPa"
    )
    print(
        "Kinematic reduction ratio: "
        f"{100.0 * metrics.kinematic_reduction_ratio:.6f}%"
    )
    print("-" * 106)
    print(
        "Interpretation: the kinematic reduction isolates "
        "yield-surface translation. The observed initial-to-reverse "
        "change also contains isotropic hardening, damage, time-step, "
        "and viscoplastic overstress effects."
    )
    print("=" * 106)


def plot_results(
    response: NonlinearReversedResponse,
    metrics: BauschingerMetrics,
) -> None:
    """Plot fixed-fiber stress-strain and yield-function histories."""
    normalized_time = (
        response.loading.times / response.loading.period
    )

    plt.figure()
    plt.plot(
        response.critical_fiber_strains,
        response.critical_fiber_stresses,
        marker="o",
        markersize=2,
    )
    plt.scatter(
        [
            np.interp(
                metrics.initial_onset.normalized_time,
                normalized_time,
                response.critical_fiber_strains,
            ),
            np.interp(
                metrics.reverse_onset.normalized_time,
                normalized_time,
                response.critical_fiber_strains,
            ),
        ],
        [
            metrics.initial_onset.stress,
            metrics.reverse_onset.stress,
        ],
    )
    plt.xlabel("Critical-fiber total strain")
    plt.ylabel("Critical-fiber stress (MPa)")
    plt.title("Initial and reverse yield onsets")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        response.critical_fiber_yield_functions,
    )
    plt.axhline(0.0)
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Critical-fiber yield function (MPa)")
    plt.title("Yield-function zero crossings")
    plt.grid(True)
    plt.tight_layout()

    plt.figure()
    plt.plot(
        normalized_time,
        response.critical_fiber_backstresses,
    )
    plt.xlabel("Normalized time, t/T")
    plt.ylabel("Critical-fiber backstress, beta (MPa)")
    plt.title("Kinematic-hardening translation")
    plt.grid(True)
    plt.tight_layout()

    plt.show()


def parse_arguments() -> argparse.Namespace:
    """Read command-line options."""
    parser = argparse.ArgumentParser(
        description=(
            "Quantify initial and reverse yield at one tower fiber."
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
        "--increments",
        type=int,
        default=160,
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
    """Run the Bauschinger-effect diagnostic."""
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
        n_cycles=1,
        increments_per_cycle=arguments.increments,
    )
    response = run_nonlinear_reversed_analysis(
        configuration=configuration,
        material=material,
        loading=loading,
    )
    metrics = evaluate_bauschinger_metrics(response)

    print_summary(
        response=response,
        metrics=metrics,
    )

    if not arguments.no_plot:
        plot_results(
            response=response,
            metrics=metrics,
        )


if __name__ == "__main__":
    main()

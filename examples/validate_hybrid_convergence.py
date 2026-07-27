# -*- coding: utf-8 -*-
"""
Validate the hybrid LATIN-PGD stopping rule against a 100-iteration reference.

Run from the project root:
    python examples/validate_hybrid_convergence.py
"""

from __future__ import annotations

import csv
import time
from pathlib import Path
from typing import Any, Dict, List, Sequence, Tuple, cast

import numpy as np
from numpy.typing import NDArray

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
    prescribed_displacement,
)
from fem.bar_1d import BarMesh1D, create_uniform_bar_mesh
from latin.initialization import compute_elastic_initialization
from latin.pgd_solver import PGDLatinResult, solve_latin_pgd
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]
FieldMetric = Tuple[float, float, float]
RegionDamageRow = Dict[str, float]

OUTPUT_DIR = Path("results") / "convergence_validation"

STATE_FIELDS = (
    "stress",
    "elastic_strain",
    "plastic_strain",
    "plastic_strain_rate",
    "damage",
    "damage_rate",
    "alpha",
    "alpha_rate",
    "beta",
    "r_bar",
    "r_bar_rate",
    "R_bar",
    "energy_release_rate",
)


def termination_value(result: PGDLatinResult) -> str:
    """Return the user-facing termination reason."""
    return str(result.termination_reason.value)


def relative_l2_error(
    approximation: FloatArray,
    reference: FloatArray,
) -> float:
    """Return ||approximation-reference||_2 / ||reference||_2."""
    numerator = float(np.linalg.norm(approximation - reference))
    denominator = float(np.linalg.norm(reference))

    if denominator <= np.finfo(float).eps:
        return 0.0 if numerator <= np.finfo(float).eps else np.inf

    return numerator / denominator


def relative_max_error(
    approximation: FloatArray,
    reference: FloatArray,
) -> float:
    """Return the maximum absolute error normalised by max|reference|."""
    numerator = float(np.max(np.abs(approximation - reference)))
    denominator = float(np.max(np.abs(reference)))

    if denominator <= np.finfo(float).eps:
        return 0.0 if numerator <= np.finfo(float).eps else np.inf

    return numerator / denominator


def field_metrics(
    hybrid_state: LatinState,
    reference_state: LatinState,
) -> Dict[str, FieldMetric]:
    """Compare all principal state fields over the complete space-time grid."""
    metrics: Dict[str, FieldMetric] = {}

    for field_name in STATE_FIELDS:
        hybrid = np.asarray(
            getattr(hybrid_state, field_name),
            dtype=np.float64,
        )
        reference = np.asarray(
            getattr(reference_state, field_name),
            dtype=np.float64,
        )

        if hybrid.shape != reference.shape:
            raise ValueError(
                f"Shape mismatch for {field_name}: "
                f"{hybrid.shape} != {reference.shape}"
            )

        metrics[field_name] = (
            relative_l2_error(hybrid, reference),
            relative_max_error(hybrid, reference),
            float(np.max(np.abs(hybrid - reference))),
        )

    hybrid_total_strain = (
        hybrid_state.elastic_strain + hybrid_state.plastic_strain
    )
    reference_total_strain = (
        reference_state.elastic_strain
        + reference_state.plastic_strain
    )

    metrics["total_strain"] = (
        relative_l2_error(
            hybrid_total_strain,
            reference_total_strain,
        ),
        relative_max_error(
            hybrid_total_strain,
            reference_total_strain,
        ),
        float(
            np.max(
                np.abs(
                    hybrid_total_strain
                    - reference_total_strain
                )
            )
        ),
    )

    return metrics


def region_damage_rows(
    hybrid_state: LatinState,
    reference_state: LatinState,
) -> List[RegionDamageRow]:
    """Compare final mean damage in the three equal material regions."""
    n_elements = hybrid_state.n_elements

    if n_elements % 3 != 0:
        raise ValueError(
            "The validation script expects three equal material regions."
        )

    region_size = n_elements // 3
    rows: List[RegionDamageRow] = []

    for region_index in range(3):
        start = region_index * region_size
        end = (region_index + 1) * region_size

        hybrid_mean = float(
            np.mean(hybrid_state.damage[-1, start:end])
        )
        reference_mean = float(
            np.mean(reference_state.damage[-1, start:end])
        )
        absolute_difference = abs(hybrid_mean - reference_mean)

        if abs(reference_mean) <= np.finfo(float).eps:
            relative_difference = (
                0.0
                if absolute_difference <= np.finfo(float).eps
                else np.inf
            )
        else:
            relative_difference = (
                absolute_difference / abs(reference_mean)
            )

        rows.append(
            {
                "region": float(region_index + 1),
                "element_start": float(start),
                "element_end": float(end - 1),
                "hybrid_final_mean_damage": hybrid_mean,
                "reference_final_mean_damage": reference_mean,
                "absolute_difference": absolute_difference,
                "relative_difference": relative_difference,
            }
        )

    return rows


def save_state_npz(
    path: Path,
    result: PGDLatinResult,
    elapsed_seconds: float,
) -> None:
    """Save a compact numerical archive for later plotting or checking."""
    state = result.state

    arrays = {
        "time": state.time,
        "indicator_history": result.indicator_history,
        "saturation_history": result.saturation_history,
        "modes_added_history": result.modes_added_history,
        "elapsed_seconds": np.array(
            [elapsed_seconds],
            dtype=np.float64,
        ),
        "iterations": np.array(
            [result.iterations],
            dtype=np.int64,
        ),
        "basis_modes": np.array(
            [result.basis.n_modes],
            dtype=np.int64,
        ),
    }

    for field_name in STATE_FIELDS:
        arrays[field_name] = getattr(state, field_name)

    arrays["total_strain"] = (
        state.elastic_strain + state.plastic_strain
    )

    np.savez_compressed(path, **arrays)


def run_case(
    *,
    initial_state: LatinState,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    reference: bool,
) -> PGDLatinResult:
    """Run either the hybrid calculation or the 100-iteration reference."""
    common = dict(
        initial_state=initial_state,
        mesh=mesh,
        area=area,
        materials=materials,
        relaxation=0.8,
        reduced_tolerance=1.0e-3,
        max_enrichments_per_iteration=3,
    )

    if reference:
        # The stagnation counter cannot reach 101 within 100 iterations.
        # This disables the stagnation stop without changing the solver.
        return solve_latin_pgd(
            **common,
            tolerance=1.0e-12,
            max_iterations=100,
            stagnation_indicator_threshold=1.0e-12,
            stagnation_absolute_tolerance=1.0e-16,
            stagnation_required_iterations=101,
        )

    return solve_latin_pgd(
        **common,
        tolerance=1.0e-4,
        max_iterations=50,
        stagnation_indicator_threshold=1.0e-3,
        stagnation_absolute_tolerance=1.0e-6,
        stagnation_required_iterations=3,
    )


def write_summary(
    path: Path,
    hybrid_result: PGDLatinResult,
    reference_result: PGDLatinResult,
    metrics: Dict[str, FieldMetric],
    damage_rows: List[RegionDamageRow],
    hybrid_seconds: float,
    reference_seconds: float,
) -> str:
    """Build and save a plain-text validation summary."""
    stress_l2 = metrics["stress"][0]
    total_strain_l2 = metrics["total_strain"][0]
    damage_l2 = metrics["damage"][0]
    max_region_damage_error = max(
        row["relative_difference"] for row in damage_rows
    )

    speedup = (
        reference_seconds / hybrid_seconds
        if hybrid_seconds > 0.0
        else np.inf
    )

    summary_lines = [
        "LATIN-PGD hybrid convergence validation",
        "=" * 44,
        "",
        "[Hybrid calculation]",
        f"termination: {termination_value(hybrid_result)}",
        f"converged: {hybrid_result.converged}",
        f"iterations: {hybrid_result.iterations}",
        f"basis modes: {hybrid_result.basis.n_modes}",
        f"final indicator: {hybrid_result.final_indicator:.12e}",
        f"elapsed seconds: {hybrid_seconds:.6f}",
        "",
        "[100-iteration reference]",
        f"termination: {termination_value(reference_result)}",
        f"converged: {reference_result.converged}",
        f"iterations: {reference_result.iterations}",
        f"basis modes: {reference_result.basis.n_modes}",
        f"final indicator: {reference_result.final_indicator:.12e}",
        f"elapsed seconds: {reference_seconds:.6f}",
        "",
        "[Main comparison]",
        f"stress relative L2 error: {stress_l2:.12e}",
        (
            "total strain relative L2 error: "
            f"{total_strain_l2:.12e}"
        ),
        f"damage relative L2 error: {damage_l2:.12e}",
        (
            "maximum regional final-damage relative error: "
            f"{max_region_damage_error:.12e}"
        ),
        f"measured speed-up: {speedup:.6f}",
        "",
        "[Suggested acceptance limits]",
        "stress relative L2 error <= 1.0e-3",
        "total strain relative L2 error <= 1.0e-3",
        "damage relative L2 error <= 5.0e-3",
        "regional final-damage relative error <= 5.0e-3",
    ]

    summary = "\n".join(summary_lines)

    with path.open("w", encoding="utf-8", newline="\n") as file:
        file.write(summary)

    return summary


def main() -> None:
    """Run both calculations, compare them and save all outputs."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    configuration = BenchmarkConfiguration()
    mesh = create_uniform_bar_mesh(
        configuration.length,
        configuration.n_elements,
    )
    materials = create_three_material_distribution(
        configuration.n_elements
    )
    time_grid = create_time_grid(
        configuration.total_time,
        configuration.time_step,
    )
    right_displacement = prescribed_displacement(
        time_grid,
        configuration.displacement_amplitude,
        configuration.period,
    )

    elastic = compute_elastic_initialization(
        mesh=mesh,
        area=configuration.area,
        materials=materials,
        time=time_grid,
        right_displacement=right_displacement,
    )
    initial_state = LatinState.from_elastic_initialization(
        elastic,
        materials,
    )

    print("Running hybrid stopping-rule calculation...")
    start = time.perf_counter()
    hybrid_result = run_case(
        initial_state=initial_state,
        mesh=mesh,
        area=configuration.area,
        materials=materials,
        reference=False,
    )
    hybrid_seconds = time.perf_counter() - start

    print("Running 100-iteration reference calculation...")
    start = time.perf_counter()
    reference_result = run_case(
        initial_state=initial_state,
        mesh=mesh,
        area=configuration.area,
        materials=materials,
        reference=True,
    )
    reference_seconds = time.perf_counter() - start

    metrics = field_metrics(
        hybrid_result.state,
        reference_result.state,
    )
    damage_rows = region_damage_rows(
        hybrid_result.state,
        reference_result.state,
    )

    with (OUTPUT_DIR / "field_error_metrics.csv").open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "field",
                "relative_l2_error",
                "relative_max_error",
                "maximum_absolute_error",
            ]
        )
        for field_name, values in metrics.items():
            writer.writerow([field_name, *values])

    with (OUTPUT_DIR / "region_damage_comparison.csv").open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as file:
        fieldnames = list(damage_rows[0].keys())
        writer = csv.DictWriter(
            cast(Any, file),
            fieldnames=fieldnames,
        )
        writer.writeheader()
        writer.writerows(damage_rows)

    save_state_npz(
        OUTPUT_DIR / "hybrid_result.npz",
        hybrid_result,
        hybrid_seconds,
    )
    save_state_npz(
        OUTPUT_DIR / "reference_result.npz",
        reference_result,
        reference_seconds,
    )

    summary = write_summary(
        OUTPUT_DIR / "validation_summary.txt",
        hybrid_result,
        reference_result,
        metrics,
        damage_rows,
        hybrid_seconds,
        reference_seconds,
    )

    print()
    print(summary)
    print()
    print(f"Files written to: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    main()

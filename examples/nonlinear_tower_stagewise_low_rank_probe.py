# -*- coding: utf-8 -*-
"""
Stage-wise low-rank diagnostics for the frozen 100-cycle tower FOM dataset.

This script performs no full-order solve. It loads the frozen reference dataset

    outputs/tower_100cycle_fom_reference_v1.npz

and partitions the global cycle axis into the three mechanism-based stages

    Stage I   : cycles 1-20
    Stage II  : cycles 21-46
    Stage III : cycles 47-100

without renumbering the cycles locally.

For each stage, the existing HOSVD-style mode-wise SVD diagnostics are applied
to

    u, sigma, eps_p, D

in both raw and cycle-start-referenced increment form.

The main comparison uses the 99.99% squared-Frobenius energy threshold. The
script also prints selected singular-value ratios for the cycle-increment
plastic-strain and damage fields, because these irreversible fields carry the
most important stage-dependent complexity identified by the 100-cycle study.

This remains an empirical low-rank diagnostic. It is not a CP decomposition,
not a PGD approximation, and not a LATIN-PGD solve.

Run from the repository root:

    python -m examples.nonlinear_tower_stagewise_low_rank_probe
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Tuple

from examples.nonlinear_tower_low_rank_diagnostics import (
    FieldLowRankDiagnostics,
    SvdSpectrum,
    TowerLowRankDiagnostics,
    analyze_tower_low_rank,
)
from examples.nonlinear_tower_snapshot_tensor import (
    TowerCyclePhaseSnapshots,
    load_tower_cycle_phase_snapshots,
    select_tower_cycle_range,
)


REFERENCE_FILENAME = "tower_100cycle_fom_reference_v1.npz"

ENERGY_TARGET = 0.9999
DETAILED_RATIO_INDICES = (1, 2, 4)

FIELD_ORDER = (
    "u",
    "sigma",
    "eps_p",
    "D",
)

STAGE_DEFINITIONS = (
    ("Stage I", 1, 20),
    ("Stage II", 21, 46),
    ("Stage III", 47, 100),
)


@dataclass(frozen=True)
class StageResult:
    """One globally numbered stage and its low-rank diagnostics."""

    name: str
    first_cycle: int
    last_cycle: int
    snapshots: TowerCyclePhaseSnapshots
    diagnostics: TowerLowRankDiagnostics


def _reference_path() -> Path:
    """Return the local frozen 100-cycle FOM reference path."""
    repository_root = Path(__file__).resolve().parents[1]
    return repository_root / "outputs" / REFERENCE_FILENAME


def _validate_reference(
    snapshots: TowerCyclePhaseSnapshots,
) -> None:
    """Check that the loaded dataset covers the formal 100-cycle study."""
    first_cycle = int(snapshots.cycle_numbers[0])
    last_cycle = int(snapshots.cycle_numbers[-1])

    if first_cycle != 1 or last_cycle != 100:
        raise ValueError(
            "The stage-wise probe requires a frozen reference covering "
            "global cycles 1 through 100. Available range is [{}, {}].".format(
                first_cycle,
                last_cycle,
            )
        )

    if snapshots.n_cycles != 100:
        raise ValueError(
            "The stage-wise probe requires exactly 100 stored cycles."
        )

    if snapshots.n_phase_points != 41:
        raise ValueError(
            "The formal 100-cycle reference is expected to contain "
            "41 phase points per cycle."
        )


def _build_stage_results(
    reference: TowerCyclePhaseSnapshots,
) -> Tuple[StageResult, ...]:
    """Select the three global cycle stages and analyze each one."""
    results = []

    for name, first_cycle, last_cycle in STAGE_DEFINITIONS:
        stage_snapshots = select_tower_cycle_range(
            snapshots=reference,
            first_cycle=first_cycle,
            last_cycle=last_cycle,
        )
        diagnostics = analyze_tower_low_rank(
            stage_snapshots
        )

        results.append(
            StageResult(
                name=name,
                first_cycle=first_cycle,
                last_cycle=last_cycle,
                snapshots=stage_snapshots,
                diagnostics=diagnostics,
            )
        )

    return tuple(results)


def _representation(
    field: FieldLowRankDiagnostics,
    name: str,
):
    """Return raw or cycle-increment tensor diagnostics."""
    if name == "raw":
        return field.raw
    if name == "cycle_increment":
        return field.cycle_increment

    raise ValueError(
        "name must be 'raw' or 'cycle_increment'."
    )


def _rank_triplet(
    field: FieldLowRankDiagnostics,
    representation_name: str,
    target_energy: float = ENERGY_TARGET,
) -> Tuple[int, int, int]:
    """Return (cycle, phase, space) ranks at one energy threshold."""
    tensor = _representation(
        field=field,
        name=representation_name,
    )

    return (
        tensor.cycle_mode.rank_for_energy(target_energy),
        tensor.phase_mode.rank_for_energy(target_energy),
        tensor.space_mode.rank_for_energy(target_energy),
    )


def _ratio_at(
    spectrum: SvdSpectrum,
    index: int,
) -> float:
    """Return sigma[index] / sigma[0], handling short or zero spectra."""
    if spectrum.is_zero:
        return 0.0

    singular_values = spectrum.singular_values

    if index >= singular_values.size:
        return float("nan")

    leading = float(singular_values[0])

    if leading == 0.0:
        return 0.0

    return float(
        singular_values[index] / leading
    )


def _ratio_text(
    spectrum: SvdSpectrum,
    indices: Iterable[int] = DETAILED_RATIO_INDICES,
) -> str:
    """Format selected normalized singular values."""
    entries = []

    for index in indices:
        ratio = _ratio_at(
            spectrum=spectrum,
            index=index,
        )
        entries.append(
            "s{}/s1={:.3e}".format(
                index + 1,
                ratio,
            )
        )

    return "  ".join(entries)


def _field_dict(
    diagnostics: TowerLowRankDiagnostics,
) -> Dict[str, FieldLowRankDiagnostics]:
    """Return diagnostics in the fixed physical field order."""
    fields = diagnostics.as_dict()
    return {
        name: fields[name]
        for name in FIELD_ORDER
    }


def _print_reference_summary(
    path: Path,
    reference: TowerCyclePhaseSnapshots,
) -> None:
    """Print the frozen reference dataset anchors."""
    print("=" * 118)
    print("Frozen 100-cycle tower FOM: stage-wise low-rank diagnostics")
    print("=" * 118)
    print("Reference: {}".format(path))
    print(
        "Global cycle range: {}-{}  |  cycles={}  |  phase points={}".format(
            int(reference.cycle_numbers[0]),
            int(reference.cycle_numbers[-1]),
            reference.n_cycles,
            reference.n_phase_points,
        )
    )
    print(
        "u shape={}  |  sigma shape={}  |  state shape={}".format(
            reference.nodal_displacements.shape,
            reference.fiber_stresses.shape,
            reference.fiber_states.shape,
        )
    )
    print(
        "Energy threshold used for the comparison table: {:.4f}%".format(
            100.0 * ENERGY_TARGET
        )
    )


def _print_stage_summary(
    stage: StageResult,
) -> None:
    """Print stage size and preserved global timing anchors."""
    snapshots = stage.snapshots

    print("\n" + "-" * 118)
    print(
        "{} | global cycles {}-{} | {} cycles".format(
            stage.name,
            stage.first_cycle,
            stage.last_cycle,
            snapshots.n_cycles,
        )
    )
    print(
        "analysis time range: {:.6g} -> {:.6g}".format(
            float(snapshots.analysis_times[0, 0]),
            float(snapshots.analysis_times[-1, -1]),
        )
    )


def _print_rank_comparison(
    stages: Tuple[StageResult, ...],
) -> None:
    """Print 99.99%-energy multilinear ranks for all stages and fields."""
    print("\n" + "=" * 118)
    print("99.99%-energy mode-rank comparison")
    print("=" * 118)
    print(
        "Each tuple is (cycle rank, phase rank, space rank). "
        "These are HOSVD-style mode ranks, not CP/PGD ranks."
    )

    for representation_name in (
        "raw",
        "cycle_increment",
    ):
        print(
            "\nRepresentation: {}".format(
                representation_name
            )
        )
        print("-" * 118)
        print(
            "{:<10s} {:<20s} {:<20s} {:<20s}".format(
                "field",
                "Stage I (1-20)",
                "Stage II (21-46)",
                "Stage III (47-100)",
            )
        )

        stage_fields = [
            _field_dict(stage.diagnostics)
            for stage in stages
        ]

        for field_name in FIELD_ORDER:
            rank_values = [
                _rank_triplet(
                    field=fields[field_name],
                    representation_name=representation_name,
                )
                for fields in stage_fields
            ]

            print(
                "{:<10s} {:<20s} {:<20s} {:<20s}".format(
                    field_name,
                    str(rank_values[0]),
                    str(rank_values[1]),
                    str(rank_values[2]),
                )
            )


def _print_irreversible_increment_details(
    stages: Tuple[StageResult, ...],
) -> None:
    """
    Print singular-value decay for Delta eps_p and Delta D.

    These are the fields for which the full 100-cycle study showed the largest
    cycle/space complexity in the cycle-increment representation.
    """
    print("\n" + "=" * 118)
    print("Detailed singular-value decay: cycle-increment irreversible fields")
    print("=" * 118)

    for field_name in ("eps_p", "D"):
        print(
            "\n{} | cycle_increment".format(
                field_name
            )
        )
        print("-" * 118)

        for stage in stages:
            field = _field_dict(
                stage.diagnostics
            )[field_name]
            tensor = field.cycle_increment

            cycle_rank = tensor.cycle_mode.rank_for_energy(
                ENERGY_TARGET
            )
            phase_rank = tensor.phase_mode.rank_for_energy(
                ENERGY_TARGET
            )
            space_rank = tensor.space_mode.rank_for_energy(
                ENERGY_TARGET
            )

            print(
                "{} | cycles {:>3d}-{:>3d}".format(
                    stage.name,
                    stage.first_cycle,
                    stage.last_cycle,
                )
            )
            print(
                "  cycle: rank={:<2d}  {}".format(
                    cycle_rank,
                    _ratio_text(tensor.cycle_mode),
                )
            )
            print(
                "  phase: rank={:<2d}  {}".format(
                    phase_rank,
                    _ratio_text(tensor.phase_mode),
                )
            )
            print(
                "  space: rank={:<2d}  {}".format(
                    space_rank,
                    _ratio_text(tensor.space_mode),
                )
            )


def main() -> None:
    """Load the frozen FOM and compare Stage I/II/III low-rank structure."""
    reference_path = _reference_path()

    if not reference_path.is_file():
        raise FileNotFoundError(
            "Frozen 100-cycle FOM reference not found: {}".format(
                reference_path
            )
        )

    reference = load_tower_cycle_phase_snapshots(
        reference_path
    )
    _validate_reference(reference)

    stages = _build_stage_results(reference)

    _print_reference_summary(
        path=reference_path,
        reference=reference,
    )

    for stage in stages:
        _print_stage_summary(stage)

    _print_rank_comparison(stages)
    _print_irreversible_increment_details(stages)

    print("\n" + "=" * 118)
    print(
        "Stage-wise analysis complete. No FOM solve was performed."
    )
    print("=" * 118)


if __name__ == "__main__":
    main()

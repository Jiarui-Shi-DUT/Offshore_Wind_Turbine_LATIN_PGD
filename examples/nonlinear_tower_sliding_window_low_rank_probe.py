# -*- coding: utf-8 -*-
"""
Sliding-window low-rank diagnostics for the frozen 100-cycle tower FOM.

Purpose
-------
The equal-window control showed that the low-rank structure changes across the
cycle-46 transition even when the compared windows have identical lengths.

This script resolves that transition more finely by moving a fixed 20-cycle
window in 5-cycle steps:

    1-20, 6-25, 11-30, ..., 81-100

For each window, the script analyzes the cycle-increment irreversible fields

    Delta eps_p
    Delta D

and records, for the cycle / phase / space unfoldings:

    - rank at 99.99% squared-Frobenius energy
    - s2 / s1

The resulting sequence provides a modal-complexity indicator versus cycle
number and is intended to support later design of an adaptive PGD enrichment
criterion.

The script performs no full-order solve. It loads

    outputs/tower_100cycle_fom_reference_v1.npz

and writes a compact CSV table to

    outputs/tower_100cycle_sliding_window_low_rank_v1.csv

The reported ranks remain HOSVD-style mode ranks; they are not CP/PGD
separation ranks.

Run from the repository root:

    python -m examples.nonlinear_tower_sliding_window_low_rank_probe
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple

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
OUTPUT_FILENAME = "tower_100cycle_sliding_window_low_rank_v1.csv"

WINDOW_LENGTH = 20
WINDOW_STEP = 5
ENERGY_TARGET = 0.9999

TARGET_FIELDS = (
    "eps_p",
    "D",
)

MODE_NAMES = (
    "cycle",
    "phase",
    "space",
)


@dataclass(frozen=True)
class SlidingWindowResult:
    """One fixed-length sliding cycle window and its diagnostics."""

    first_cycle: int
    last_cycle: int
    snapshots: TowerCyclePhaseSnapshots
    diagnostics: TowerLowRankDiagnostics

    @property
    def center_cycle(self) -> float:
        """Return the geometric center of the inclusive cycle window."""
        return 0.5 * (
            float(self.first_cycle)
            + float(self.last_cycle)
        )


def _repository_root() -> Path:
    """Return the repository root from this example module."""
    return Path(__file__).resolve().parents[1]


def _reference_path() -> Path:
    """Return the frozen FOM reference path."""
    return (
        _repository_root()
        / "outputs"
        / REFERENCE_FILENAME
    )


def _output_path() -> Path:
    """Return the CSV output path."""
    return (
        _repository_root()
        / "outputs"
        / OUTPUT_FILENAME
    )


def _validate_reference(
    snapshots: TowerCyclePhaseSnapshots,
) -> None:
    """Validate the formal frozen 100-cycle reference."""
    if snapshots.n_cycles != 100:
        raise ValueError(
            "Sliding-window diagnostics require exactly 100 stored cycles."
        )

    if int(snapshots.cycle_numbers[0]) != 1:
        raise ValueError(
            "The frozen reference must start at global cycle 1."
        )

    if int(snapshots.cycle_numbers[-1]) != 100:
        raise ValueError(
            "The frozen reference must end at global cycle 100."
        )

    if snapshots.n_phase_points != 41:
        raise ValueError(
            "The formal frozen reference must contain 41 phase points "
            "per cycle."
        )


def _window_ranges() -> Tuple[Tuple[int, int], ...]:
    """Return all fixed 20-cycle windows shifted by 5 cycles."""
    last_start = (
        100
        - WINDOW_LENGTH
        + 1
    )

    starts = range(
        1,
        last_start + 1,
        WINDOW_STEP,
    )

    windows = tuple(
        (
            start,
            start + WINDOW_LENGTH - 1,
        )
        for start in starts
    )

    if windows[-1] != (
        81,
        100,
    ):
        raise RuntimeError(
            "Unexpected sliding-window definition."
        )

    return windows


def _build_results(
    reference: TowerCyclePhaseSnapshots,
) -> Tuple[SlidingWindowResult, ...]:
    """Select and analyze every fixed-length sliding window."""
    results: List[SlidingWindowResult] = []

    for first_cycle, last_cycle in _window_ranges():
        snapshots = select_tower_cycle_range(
            snapshots=reference,
            first_cycle=first_cycle,
            last_cycle=last_cycle,
        )

        diagnostics = analyze_tower_low_rank(
            snapshots
        )

        results.append(
            SlidingWindowResult(
                first_cycle=first_cycle,
                last_cycle=last_cycle,
                snapshots=snapshots,
                diagnostics=diagnostics,
            )
        )

    return tuple(results)


def _field_dict(
    diagnostics: TowerLowRankDiagnostics,
) -> Dict[str, FieldLowRankDiagnostics]:
    """Return low-rank fields under compact physical names."""
    return diagnostics.as_dict()


def _ratio_at(
    spectrum: SvdSpectrum,
    index: int,
) -> float:
    """Return sigma[index] / sigma[0], handling zero spectra."""
    if spectrum.is_zero:
        return 0.0

    singular_values = spectrum.singular_values

    if index >= singular_values.size:
        return float("nan")

    leading = float(
        singular_values[0]
    )

    if leading == 0.0:
        return 0.0

    return float(
        singular_values[index]
        / leading
    )


def _window_metrics(
    result: SlidingWindowResult,
) -> Dict[str, float]:
    """
    Return one flat record of cycle-increment ranks and s2/s1 values.
    """
    record: Dict[str, float] = {
        "first_cycle": float(result.first_cycle),
        "last_cycle": float(result.last_cycle),
        "center_cycle": result.center_cycle,
    }

    fields = _field_dict(
        result.diagnostics
    )

    for field_name in TARGET_FIELDS:
        tensor = fields[
            field_name
        ].cycle_increment

        for mode_name in MODE_NAMES:
            spectrum = tensor.mode(
                mode_name
            )

            rank_key = (
                field_name
                + "_"
                + mode_name
                + "_rank_9999"
            )
            ratio_key = (
                field_name
                + "_"
                + mode_name
                + "_s2_s1"
            )

            record[rank_key] = float(
                spectrum.rank_for_energy(
                    ENERGY_TARGET
                )
            )
            record[ratio_key] = _ratio_at(
                spectrum=spectrum,
                index=1,
            )

    return record


def _all_records(
    results: Tuple[SlidingWindowResult, ...],
) -> Tuple[Dict[str, float], ...]:
    """Convert all window diagnostics to flat records."""
    return tuple(
        _window_metrics(result)
        for result in results
    )


def _write_csv(
    records: Tuple[Dict[str, float], ...],
    path: Path,
) -> None:
    """Write all sliding-window metrics to a CSV file."""
    if not records:
        raise ValueError(
            "At least one sliding-window record is required."
        )

    path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    fieldnames = list(
        records[0].keys()
    )

    with path.open(
        "w",
        encoding="utf-8",
        newline="",
    ) as file_object:
        writer = csv.DictWriter(
            file_object,
            fieldnames=fieldnames,
        )
        writer.writeheader()

        for record in records:
            writer.writerow(
                record
            )


def _print_reference_summary(
    reference_path: Path,
    output_path: Path,
) -> None:
    """Print the sliding-window setup."""
    print("=" * 118)
    print("Frozen 100-cycle tower FOM: 20-cycle sliding-window low-rank diagnostics")
    print("=" * 118)
    print(
        "Reference: {}".format(
            reference_path
        )
    )
    print(
        "Window length={} cycles | step={} cycles | number of windows={}".format(
            WINDOW_LENGTH,
            WINDOW_STEP,
            len(_window_ranges()),
        )
    )
    print(
        "Energy threshold: {:.4f}%".format(
            100.0 * ENERGY_TARGET
        )
    )
    print(
        "CSV output: {}".format(
            output_path
        )
    )


def _print_compact_table(
    records: Tuple[Dict[str, float], ...],
    field_name: str,
) -> None:
    """Print one compact modal-complexity table for one field."""
    print("\n" + "=" * 118)
    print(
        "{} | cycle_increment | sliding-window modal complexity".format(
            field_name
        )
    )
    print("=" * 118)

    print(
        "{:<9s} {:>6s} | "
        "{:>5s} {:>10s} | "
        "{:>5s} {:>10s} | "
        "{:>5s} {:>10s}".format(
            "window",
            "center",
            "r_n",
            "s2/s1_n",
            "r_tau",
            "s2/s1_tau",
            "r_x",
            "s2/s1_x",
        )
    )
    print("-" * 118)

    for record in records:
        first_cycle = int(
            record["first_cycle"]
        )
        last_cycle = int(
            record["last_cycle"]
        )

        print(
            "{:>3d}-{:>3d} {:>6.1f} | "
            "{:>5d} {:>10.3e} | "
            "{:>5d} {:>10.3e} | "
            "{:>5d} {:>10.3e}".format(
                first_cycle,
                last_cycle,
                record["center_cycle"],
                int(
                    record[
                        field_name
                        + "_cycle_rank_9999"
                    ]
                ),
                record[
                    field_name
                    + "_cycle_s2_s1"
                ],
                int(
                    record[
                        field_name
                        + "_phase_rank_9999"
                    ]
                ),
                record[
                    field_name
                    + "_phase_s2_s1"
                ],
                int(
                    record[
                        field_name
                        + "_space_rank_9999"
                    ]
                ),
                record[
                    field_name
                    + "_space_s2_s1"
                ],
            )
        )


def _print_transition_neighborhood(
    records: Tuple[Dict[str, float], ...],
) -> None:
    """
    Print the windows that straddle or closely surround the cycle-46 transition.
    """
    selected = tuple(
        record
        for record in records
        if 26 <= int(record["first_cycle"]) <= 46
    )

    print("\n" + "=" * 118)
    print("Transition neighborhood around cycle 46")
    print("=" * 118)
    print(
        "Focused windows: 26-45, 31-50, 36-55, 41-60, 46-65"
    )

    for field_name in TARGET_FIELDS:
        print(
            "\n{}".format(
                field_name
            )
        )
        print(
            "{:<9s} {:>10s} {:>10s} {:>10s}".format(
                "window",
                "s2/s1_n",
                "s2/s1_tau",
                "s2/s1_x",
            )
        )
        print("-" * 56)

        for record in selected:
            first_cycle = int(
                record["first_cycle"]
            )
            last_cycle = int(
                record["last_cycle"]
            )

            print(
                "{:>3d}-{:>3d} {:>10.3e} {:>10.3e} {:>10.3e}".format(
                    first_cycle,
                    last_cycle,
                    record[
                        field_name
                        + "_cycle_s2_s1"
                    ],
                    record[
                        field_name
                        + "_phase_s2_s1"
                    ],
                    record[
                        field_name
                        + "_space_s2_s1"
                    ],
                )
            )


def main() -> None:
    """Run the 20-cycle sliding-window offline low-rank diagnostics."""
    reference_path = _reference_path()
    output_path = _output_path()

    if not reference_path.is_file():
        raise FileNotFoundError(
            "Frozen 100-cycle FOM reference not found: {}".format(
                reference_path
            )
        )

    reference = load_tower_cycle_phase_snapshots(
        reference_path
    )

    _validate_reference(
        reference
    )

    results = _build_results(
        reference
    )

    records = _all_records(
        results
    )

    _write_csv(
        records=records,
        path=output_path,
    )

    _print_reference_summary(
        reference_path=reference_path,
        output_path=output_path,
    )

    for field_name in TARGET_FIELDS:
        _print_compact_table(
            records=records,
            field_name=field_name,
        )

    _print_transition_neighborhood(
        records
    )

    print("\n" + "=" * 118)
    print(
        "Sliding-window analysis complete. No FOM solve was performed."
    )
    print("=" * 118)


if __name__ == "__main__":
    main()

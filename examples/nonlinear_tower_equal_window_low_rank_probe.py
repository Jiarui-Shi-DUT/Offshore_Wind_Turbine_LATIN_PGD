# -*- coding: utf-8 -*-
"""
Equal-window low-rank control diagnostics for the frozen 100-cycle tower FOM.

Purpose
-------
The previous Stage I / II / III comparison used windows with different cycle
counts:

    Stage I   : 1-20   -> 20 cycles
    Stage II  : 21-46  -> 26 cycles
    Stage III : 47-100 -> 54 cycles

Because a longer cycle window can itself increase the observed cycle-mode rank,
this script performs a controlled comparison using equal 20-cycle windows:

    W1 : cycles 1-20
    W2 : cycles 27-46
    W3 : cycles 47-66
    W4 : cycles 81-100

All four windows contain exactly 20 cycles.

The comparison is focused on the cycle-increment irreversible fields

    Delta eps_p
    Delta D

because the previous 100-cycle diagnostics showed that their cycle and spatial
directions carry the strongest stage-dependent low-rank complexity.

The script loads the frozen reference dataset

    outputs/tower_100cycle_fom_reference_v1.npz

and performs no full-order solve.

This remains an empirical HOSVD-style diagnostic. The reported cycle/phase/
space ranks are not CP/PGD separation ranks.

Run from the repository root:

    python -m examples.nonlinear_tower_equal_window_low_rank_probe
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

ENERGY_TARGETS = (
    0.99,
    0.999,
    0.9999,
)

RATIO_INDICES = (
    1,
    2,
    4,
)

WINDOW_DEFINITIONS = (
    (
        "W1",
        "early adjustment",
        1,
        20,
    ),
    (
        "W2",
        "late Stage II / pre-transition",
        27,
        46,
    ),
    (
        "W3",
        "early Stage III / post-transition",
        47,
        66,
    ),
    (
        "W4",
        "late Stage III",
        81,
        100,
    ),
)

TARGET_FIELDS = (
    "eps_p",
    "D",
)


@dataclass(frozen=True)
class WindowResult:
    """One fixed-length global-cycle window and its diagnostics."""

    name: str
    description: str
    first_cycle: int
    last_cycle: int
    snapshots: TowerCyclePhaseSnapshots
    diagnostics: TowerLowRankDiagnostics

    @property
    def n_cycles(self) -> int:
        """Return the number of cycles in the window."""
        return self.last_cycle - self.first_cycle + 1


def _reference_path() -> Path:
    """Return the frozen 100-cycle FOM reference path."""
    repository_root = Path(__file__).resolve().parents[1]
    return repository_root / "outputs" / REFERENCE_FILENAME


def _validate_reference(
    snapshots: TowerCyclePhaseSnapshots,
) -> None:
    """Validate the formal frozen 100-cycle reference dataset."""
    if snapshots.n_cycles != 100:
        raise ValueError(
            "Equal-window diagnostics require exactly 100 stored cycles."
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
            "The formal reference must contain 41 phase points per cycle."
        )


def _validate_window_definitions() -> None:
    """Require every diagnostic window to have exactly the same length."""
    lengths = tuple(
        last_cycle - first_cycle + 1
        for _, _, first_cycle, last_cycle
        in WINDOW_DEFINITIONS
    )

    if len(set(lengths)) != 1:
        raise ValueError(
            "All equal-window definitions must contain the same number "
            "of cycles."
        )

    if lengths[0] != 20:
        raise ValueError(
            "The formal equal-window control uses exactly 20 cycles per "
            "window."
        )


def _build_window_results(
    reference: TowerCyclePhaseSnapshots,
) -> Tuple[WindowResult, ...]:
    """Select all equal-length windows and compute their diagnostics."""
    results = []

    for (
        name,
        description,
        first_cycle,
        last_cycle,
    ) in WINDOW_DEFINITIONS:
        snapshots = select_tower_cycle_range(
            snapshots=reference,
            first_cycle=first_cycle,
            last_cycle=last_cycle,
        )

        diagnostics = analyze_tower_low_rank(
            snapshots
        )

        results.append(
            WindowResult(
                name=name,
                description=description,
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
    """Return sigma[index] / sigma[0]."""
    if spectrum.is_zero:
        return 0.0

    singular_values = spectrum.singular_values

    if index >= singular_values.size:
        return float("nan")

    leading = float(singular_values[0])

    if leading == 0.0:
        return 0.0

    return float(
        singular_values[index]
        / leading
    )


def _ratio_text(
    spectrum: SvdSpectrum,
    indices: Iterable[int] = RATIO_INDICES,
) -> str:
    """Format selected normalized singular values."""
    entries = []

    for index in indices:
        entries.append(
            "s{}/s1={:.3e}".format(
                index + 1,
                _ratio_at(
                    spectrum=spectrum,
                    index=index,
                ),
            )
        )

    return "  ".join(entries)


def _rank_text(
    spectrum: SvdSpectrum,
) -> str:
    """Format ranks at 99%, 99.9%, and 99.99% energy."""
    entries = []

    for target in ENERGY_TARGETS:
        entries.append(
            "{:.2f}%:{}".format(
                100.0 * target,
                spectrum.rank_for_energy(target),
            )
        )

    return "  ".join(entries)


def _print_reference_summary(
    path: Path,
    reference: TowerCyclePhaseSnapshots,
) -> None:
    """Print the frozen FOM anchors and window definitions."""
    print("=" * 118)
    print("Frozen 100-cycle tower FOM: equal-window low-rank control")
    print("=" * 118)
    print("Reference: {}".format(path))
    print(
        "Global cycle range: {}-{} | cycles={} | phase points={}".format(
            int(reference.cycle_numbers[0]),
            int(reference.cycle_numbers[-1]),
            reference.n_cycles,
            reference.n_phase_points,
        )
    )
    print(
        "All comparison windows contain exactly 20 cycles."
    )

    for (
        name,
        description,
        first_cycle,
        last_cycle,
    ) in WINDOW_DEFINITIONS:
        print(
            "  {}: cycles {}-{} | {}".format(
                name,
                first_cycle,
                last_cycle,
                description,
            )
        )


def _print_window_time_ranges(
    windows: Tuple[WindowResult, ...],
) -> None:
    """Print preserved absolute analysis-time ranges."""
    print("\n" + "=" * 118)
    print("Window timing anchors")
    print("=" * 118)

    for window in windows:
        print(
            "{} | cycles {:>3d}-{:>3d} | "
            "analysis time {:.6g} -> {:.6g}".format(
                window.name,
                window.first_cycle,
                window.last_cycle,
                float(
                    window.snapshots.analysis_times[0, 0]
                ),
                float(
                    window.snapshots.analysis_times[-1, -1]
                ),
            )
        )


def _print_rank_table(
    windows: Tuple[WindowResult, ...],
) -> None:
    """
    Print 99.99%-energy cycle/phase/space ranks for Delta eps_p and Delta D.
    """
    print("\n" + "=" * 118)
    print("99.99%-energy mode-rank comparison: cycle-increment fields")
    print("=" * 118)
    print(
        "Each tuple is (cycle rank, phase rank, space rank). "
        "All windows have equal cycle length."
    )

    header = "{:<10s}".format("field")

    for window in windows:
        header += " {:<22s}".format(
            "{} ({}-{})".format(
                window.name,
                window.first_cycle,
                window.last_cycle,
            )
        )

    print("\n" + header)
    print("-" * 118)

    for field_name in TARGET_FIELDS:
        row = "{:<10s}".format(field_name)

        for window in windows:
            field = _field_dict(
                window.diagnostics
            )[field_name]

            tensor = field.cycle_increment

            ranks = (
                tensor.cycle_mode.rank_for_energy(0.9999),
                tensor.phase_mode.rank_for_energy(0.9999),
                tensor.space_mode.rank_for_energy(0.9999),
            )

            row += " {:<22s}".format(
                str(ranks)
            )

        print(row)


def _print_detailed_spectra(
    windows: Tuple[WindowResult, ...],
) -> None:
    """Print detailed singular-value decay for each irreversible field."""
    print("\n" + "=" * 118)
    print("Detailed equal-window singular-value decay")
    print("=" * 118)

    for field_name in TARGET_FIELDS:
        print(
            "\n{} | cycle_increment".format(
                field_name
            )
        )
        print("-" * 118)

        for window in windows:
            field = _field_dict(
                window.diagnostics
            )[field_name]
            tensor = field.cycle_increment

            print(
                "{} | cycles {:>3d}-{:>3d} | {}".format(
                    window.name,
                    window.first_cycle,
                    window.last_cycle,
                    window.description,
                )
            )

            for mode_name in (
                "cycle",
                "phase",
                "space",
            ):
                spectrum = tensor.mode(
                    mode_name
                )

                print(
                    "  {:<6s} | {:<34s} | {}".format(
                        mode_name,
                        _rank_text(spectrum),
                        _ratio_text(spectrum),
                    )
                )


def _print_transition_control(
    windows: Tuple[WindowResult, ...],
) -> None:
    """
    Print the direct W2 -> W3 comparison around the cycle-46 transition.

    W2 ends at cycle 46 and W3 starts at cycle 47, so these two equal-length
    windows form the cleanest control for testing whether post-transition
    complexity increases independently of window length.
    """
    lookup = {
        window.name: window
        for window in windows
    }

    w2 = lookup["W2"]
    w3 = lookup["W3"]

    print("\n" + "=" * 118)
    print("Direct transition control: W2 (27-46) -> W3 (47-66)")
    print("=" * 118)

    for field_name in TARGET_FIELDS:
        field_w2 = _field_dict(
            w2.diagnostics
        )[field_name].cycle_increment

        field_w3 = _field_dict(
            w3.diagnostics
        )[field_name].cycle_increment

        print(
            "\n{}".format(
                field_name
            )
        )

        for mode_name in (
            "cycle",
            "phase",
            "space",
        ):
            spectrum_w2 = field_w2.mode(
                mode_name
            )
            spectrum_w3 = field_w3.mode(
                mode_name
            )

            rank_w2 = spectrum_w2.rank_for_energy(
                0.9999
            )
            rank_w3 = spectrum_w3.rank_for_energy(
                0.9999
            )

            second_ratio_w2 = _ratio_at(
                spectrum=spectrum_w2,
                index=1,
            )
            second_ratio_w3 = _ratio_at(
                spectrum=spectrum_w3,
                index=1,
            )

            if second_ratio_w2 == 0.0:
                ratio_change = float("inf")
            else:
                ratio_change = (
                    second_ratio_w3
                    / second_ratio_w2
                )

            print(
                "  {:<6s} | rank {:>2d} -> {:>2d} | "
                "s2/s1 {:.3e} -> {:.3e} | factor={:.3f}".format(
                    mode_name,
                    rank_w2,
                    rank_w3,
                    second_ratio_w2,
                    second_ratio_w3,
                    ratio_change,
                )
            )


def main() -> None:
    """Run the equal-window offline low-rank control analysis."""
    _validate_window_definitions()

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

    _validate_reference(
        reference
    )

    windows = _build_window_results(
        reference
    )

    _print_reference_summary(
        path=reference_path,
        reference=reference,
    )
    _print_window_time_ranges(
        windows
    )
    _print_rank_table(
        windows
    )
    _print_detailed_spectra(
        windows
    )
    _print_transition_control(
        windows
    )

    print("\n" + "=" * 118)
    print(
        "Equal-window control complete. No FOM solve was performed."
    )
    print("=" * 118)


if __name__ == "__main__":
    main()

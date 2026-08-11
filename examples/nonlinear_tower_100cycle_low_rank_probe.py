# -*- coding: utf-8 -*-
"""
100-cycle FOM low-rank probe for the NREL 5 MW tower.

This script runs the current reduced full-order tower model under the formal
asymmetric cyclic loading used in the 100-cycle transition study, tensorizes
the response into slow-cycle x fast-phase coordinates, freezes the complete
cycle-phase snapshots as a compressed NPZ reference dataset, and prints
mode-wise SVD diagnostics for

    u, sigma, eps_p, D.

The analysis is diagnostic only. It does not construct or solve a LATIN-PGD
reduced model.

Run from the repository root:

    python -m examples.nonlinear_tower_100cycle_low_rank_probe

Formal probe discretisation
---------------------------
Tower:
    10 beam elements
    2 Gauss points / element
    16 circumferential fibers
    1 radial fiber layer

Loading:
    F_max = 1.0 MN
    R_F = F_min / F_max = -0.5
    T = 10
    100 cycles
    40 increments / cycle

The preload used by ``run_nonlinear_asymmetric_analysis`` is retained exactly
as in the previous asymmetric-cycle studies and is not counted as a cycle.

Frozen reference dataset
------------------------
The complete ``TowerCyclePhaseSnapshots`` object is saved to

    outputs/tower_100cycle_fom_reference_v1.npz

The ``outputs`` directory is already excluded by the repository ``.gitignore``.
The binary FOM reference is therefore kept locally and is not intended to be
committed to Git.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Iterable, Tuple

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_response import (
    run_nonlinear_asymmetric_analysis,
)
from examples.nonlinear_tower_low_rank_diagnostics import (
    FieldLowRankDiagnostics,
    SvdSpectrum,
    analyze_tower_low_rank,
)
from examples.nonlinear_tower_snapshot_tensor import (
    build_tower_cycle_phase_snapshots,
    save_tower_cycle_phase_snapshots,
)
from fem.tower_loading import (
    create_asymmetric_cyclic_top_force_history,
)
from material.viscoplastic_damage_1d import MaterialParameters


ENERGY_TARGETS = (0.99, 0.999, 0.9999)
SINGULAR_RATIO_INDICES = (1, 2, 4, 9)

FROZEN_SNAPSHOT_FILENAME = "tower_100cycle_fom_reference_v1.npz"


def _default_snapshot_path() -> Path:
    """Return the local path for the frozen 100-cycle FOM reference dataset."""
    repository_root = Path(__file__).resolve().parents[1]
    return repository_root / "outputs" / FROZEN_SNAPSHOT_FILENAME


def _ratio_at(
    spectrum: SvdSpectrum,
    index: int,
) -> float:
    """Return sigma[index] / sigma[0], or zero for a zero spectrum."""
    singular_values = spectrum.singular_values

    if spectrum.is_zero:
        return 0.0
    if index >= singular_values.size:
        return float("nan")

    leading = float(singular_values[0])
    if leading == 0.0:
        return 0.0

    return float(singular_values[index] / leading)


def _rank_text(
    spectrum: SvdSpectrum,
    targets: Iterable[float],
) -> str:
    """Format energy ranks compactly."""
    entries = []
    for target in targets:
        rank = spectrum.rank_for_energy(target)
        percentage = 100.0 * target
        entries.append(
            "{:.2f}%:{}".format(percentage, rank)
        )
    return "  ".join(entries)


def _singular_ratio_text(
    spectrum: SvdSpectrum,
) -> str:
    """Format selected normalized singular values."""
    entries = []
    for index in SINGULAR_RATIO_INDICES:
        ratio = _ratio_at(spectrum, index)
        entries.append(
            "s{}/s1={:.3e}".format(index + 1, ratio)
        )
    return "  ".join(entries)


def _print_tensor_diagnostics(
    field_name: str,
    representation_name: str,
    field: FieldLowRankDiagnostics,
) -> None:
    """Print all three unfolding diagnostics for one field representation."""
    diagnostics = (
        field.raw
        if representation_name == "raw"
        else field.cycle_increment
    )

    print(
        "\n{} | {} | tensor shape = {}".format(
            field_name,
            representation_name,
            diagnostics.tensor_shape,
        )
    )
    print("-" * 118)

    for mode_name in ("cycle", "phase", "space"):
        spectrum = diagnostics.mode(mode_name)
        print(
            "{:<6s} | {:<34s} | {}".format(
                mode_name,
                _rank_text(spectrum, ENERGY_TARGETS),
                _singular_ratio_text(spectrum),
            )
        )


def _print_global_summary(
    response,
    snapshots,
    elapsed_fom: float,
    elapsed_tensorization: float,
) -> None:
    """Print FOM anchors before the persistence and SVD stages."""
    loading = response.loading
    maximum_damage = float(response.maximum_damages[-1])
    maximum_plastic = float(
        response.maximum_absolute_plastic_strains[-1]
    )

    print("=" * 118)
    print("100-cycle asymmetric tower FOM: low-rank diagnostic probe")
    print("=" * 118)
    print(
        "Discretisation: {} elements, {} Gauss/element, "
        "{} x {} fibers/section".format(
            response.fiber_strains.shape[1],
            response.fiber_strains.shape[2],
            response.fiber_strains.shape[3],
            1,
        )
    )
    print(
        "Loading: Fmax={:.6g} MN, Fmin={:.6g} MN, "
        "R_F={:.6g}, cycles={}, increments/cycle={}".format(
            loading.maximum_force / 1.0e6,
            loading.minimum_force / 1.0e6,
            loading.force_ratio,
            loading.n_cycles,
            loading.increments_per_cycle,
        )
    )
    print(
        "Periodic stored points: {}".format(
            response.analysis_times.size
        )
    )
    print(
        "Cycle-phase grid: {} cycles x {} phase points".format(
            snapshots.n_cycles,
            snapshots.n_phase_points,
        )
    )
    print(
        "Full nodal displacement tensor: {}".format(
            snapshots.nodal_displacements.shape
        )
    )
    print(
        "Fiber stress tensor: {}".format(
            snapshots.fiber_stresses.shape
        )
    )
    print(
        "Fiber state tensor: {}".format(
            snapshots.fiber_states.shape
        )
    )
    print(
        "Final max |eps_p|={:.9e}, final max D={:.9e}".format(
            maximum_plastic,
            maximum_damage,
        )
    )
    print(
        "Critical fiber: element={}, Gauss={}, fiber={}".format(
            response.critical_location[0] + 1,
            response.critical_location[1] + 1,
            response.critical_location[2] + 1,
        )
    )
    print(
        "Critical height={:.9e} m, y={:.9e} m".format(
            response.critical_height,
            response.critical_y_coordinate,
        )
    )
    print(
        "Max Newton iterations={}, max free-DOF residual={:.9e} N".format(
            int(np.max(response.iterations)),
            float(np.max(response.residual_norms)),
        )
    )
    print(
        "Elapsed: FOM={:.3f} s, tensorization={:.3f} s".format(
            elapsed_fom,
            elapsed_tensorization,
        )
    )


def _print_snapshot_summary(
    snapshot_path: Path,
    elapsed_save: float,
) -> None:
    """Print the location, size, and write time of the frozen dataset."""
    size_mib = snapshot_path.stat().st_size / (1024.0 ** 2)

    print("\n" + "=" * 118)
    print("Frozen 100-cycle FOM snapshot dataset")
    print("=" * 118)
    print("Path: {}".format(snapshot_path))
    print("File size: {:.3f} MiB".format(size_mib))
    print("Elapsed snapshot save: {:.3f} s".format(elapsed_save))
    print(
        "This NPZ is the frozen full-order reference for subsequent "
        "offline Stage I/II/III and reduced-order diagnostics."
    )


def main() -> None:
    """Run, freeze, and diagnose the formal 100-cycle tower FOM."""
    configuration = TowerConfiguration(
        horizontal_force=1.0e6,
        n_elements=10,
        n_gauss=2,
        n_circumferential=16,
        n_radial=1,
    )
    material = MaterialParameters()
    loading = create_asymmetric_cyclic_top_force_history(
        maximum_force=1.0e6,
        force_ratio=-0.5,
        period=10.0,
        n_cycles=100,
        increments_per_cycle=40,
    )

    start = time.perf_counter()
    response = run_nonlinear_asymmetric_analysis(
        configuration=configuration,
        material=material,
        loading=loading,
        max_iterations=40,
    )
    after_fom = time.perf_counter()

    snapshots = build_tower_cycle_phase_snapshots(response)
    after_tensorization = time.perf_counter()

    _print_global_summary(
        response=response,
        snapshots=snapshots,
        elapsed_fom=after_fom - start,
        elapsed_tensorization=(
            after_tensorization - after_fom
        ),
    )

    snapshot_path = _default_snapshot_path()
    snapshot_path.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    save_start = time.perf_counter()
    saved_snapshot_path = save_tower_cycle_phase_snapshots(
        snapshots=snapshots,
        file_path=snapshot_path,
    )
    after_snapshot_save = time.perf_counter()

    _print_snapshot_summary(
        snapshot_path=saved_snapshot_path,
        elapsed_save=after_snapshot_save - save_start,
    )

    diagnostics = analyze_tower_low_rank(snapshots)
    after_svd = time.perf_counter()

    print("\n" + "=" * 118)
    print("Mode-wise SVD diagnostics")
    print("=" * 118)
    print(
        "Ranks are the minimum ranks capturing "
        "99%, 99.9%, and 99.99% squared-Frobenius energy."
    )
    print(
        "Singular-value ratios show s2/s1, s3/s1, s5/s1, and s10/s1."
    )
    print(
        "These mode ranks diagnose multilinear compressibility; "
        "they are not direct CP/PGD separation ranks."
    )

    fields: Tuple[
        Tuple[str, FieldLowRankDiagnostics],
        ...,
    ] = (
        ("u", diagnostics.nodal_displacements),
        ("sigma", diagnostics.fiber_stresses),
        ("eps_p", diagnostics.fiber_plastic_strains),
        ("D", diagnostics.fiber_damages),
    )

    for field_name, field in fields:
        _print_tensor_diagnostics(
            field_name=field_name,
            representation_name="raw",
            field=field,
        )
        _print_tensor_diagnostics(
            field_name=field_name,
            representation_name="cycle_increment",
            field=field,
        )

    print("\n" + "=" * 118)
    print(
        "Elapsed SVD diagnostics={:.3f} s; total={:.3f} s".format(
            after_svd - after_snapshot_save,
            after_svd - start,
        )
    )
    print("=" * 118)


if __name__ == "__main__":
    main()

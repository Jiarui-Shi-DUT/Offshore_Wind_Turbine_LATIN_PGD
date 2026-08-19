# -*- coding: utf-8 -*-
"""Multi-iteration convergence probe for the nonzero tower LATIN-PGD solve."""

from __future__ import annotations

from time import perf_counter

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration, create_tower_geometry
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import ViscoplasticDamageTowerSystem2D
from latin.tower_equilibrium_operator import build_tower_equilibrium_operator
from latin.tower_initialization import compute_tower_elastic_initialization
import latin.tower_pgd_enrichment as enrichment_module
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from material.viscoplastic_damage_1d import MaterialParameters


def main() -> None:
    configuration = TowerConfiguration(
        horizontal_force=1.0e6,
        n_elements=10,
        n_gauss=2,
        n_circumferential=16,
        n_radial=1,
    )
    material = MaterialParameters()
    geometry = create_tower_geometry(configuration)
    mesh = create_uniform_vertical_tower_mesh(
        height=configuration.height,
        n_elements=configuration.n_elements,
    )
    system = ViscoplasticDamageTowerSystem2D(
        mesh=mesh,
        tower_geometry=geometry,
        material=material,
        n_gauss=configuration.n_gauss,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )
    operator = build_tower_equilibrium_operator(system)

    loading = create_reversed_top_force_history(
        force_amplitude=1.0e6,
        period=10.0,
        n_cycles=1,
        increments_per_cycle=40,
    )
    load_vectors = np.stack(
        [
            top_horizontal_load_vector(
                mesh=mesh,
                horizontal_force=float(force),
            )
            for force in loading.forces
        ],
        axis=0,
    )
    initialization = compute_tower_elastic_initialization(
        time=loading.times,
        load_vectors=load_vectors,
        materials=material,
        equilibrium_operator=operator,
        stress_to_force_factor=system.stress_to_force_factor,
    )

    print("=" * 112)
    print("Tower LATIN-PGD reversed multi-iteration convergence probe")
    print("=" * 112)
    print(
        "Grid: Nt={}, Nq={} ({} elements x {} Gauss x {} fibers)".format(
            initialization.state.n_time,
            initialization.state.n_material_points,
            configuration.n_elements,
            configuration.n_gauss,
            configuration.n_circumferential * configuration.n_radial,
        )
    )
    print(
        "Elastic init max |stress| = {:.9e} MPa; max equilibrium residual = {:.9e} N".format(
            float(np.max(np.abs(initialization.state.stress))),
            initialization.maximum_free_equilibrium_residual,
        )
    )
    print(
        "Diagnostic controls: mode_significance_tolerance=0, acceptance_tolerance=0, max_iterations=20, max_fixed_point_iterations=120"
    )

    raw_fixed_point_records = []
    original_raw_fixed_point = enrichment_module._raw_fixed_point
    original_pair_change = enrichment_module._pair_change

    def instrumented_raw_fixed_point(*args, **kwargs):
        modes = []

        def recording_pair_change(previous, current, time, H_sigma, metric):
            if not modes:
                modes.append(previous)
            modes.append(current)
            return original_pair_change(
                previous,
                current,
                time,
                H_sigma,
                metric,
            )

        saved_pair_change = enrichment_module._pair_change
        enrichment_module._pair_change = recording_pair_change
        try:
            raw_mode, history, converged = original_raw_fixed_point(
                *args,
                **kwargs,
            )
        finally:
            enrichment_module._pair_change = saved_pair_change

        raw_fixed_point_records.append(
            {
                "modes": tuple(modes),
                "time": np.asarray(args[1], dtype=np.float64),
                "H_sigma": np.asarray(args[3], dtype=np.float64),
                "metric": args[4],
            }
        )
        return raw_mode, history, converged

    enrichment_module._raw_fixed_point = instrumented_raw_fixed_point
    start = perf_counter()
    try:
        result = solve_tower_latin_pgd(
            initial_state=initialization.state,
            materials=material,
            metric=operator.metric,
            equilibrium_operator=operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            max_iterations=20,
            max_fixed_point_iterations=120,
        )
    finally:
        enrichment_module._raw_fixed_point = original_raw_fixed_point
    elapsed = perf_counter() - start

    print("-" * 112)
    print("Solver summary")
    print("-" * 112)
    print("termination_reason     =", result.termination_reason.value)
    print("converged              =", result.converged)
    print("failure_reason         =", result.failure_reason)
    print("attempted iterations   =", result.attempted_iterations)
    print("committed iterations   =", result.iterations)
    print("trial evaluations      =", result.trial_evaluations)
    print("final basis modes      =", result.basis.n_modes)
    print("total modes added      =", result.total_modes_added)
    print("final accepted xi      = {:.9e}".format(result.final_indicator))
    print("elapsed                = {:.6f} s".format(elapsed))

    print("-" * 112)
    print("Trial history")
    print("-" * 112)
    print(
        "{:>6s} {:>6s} {:>8s} {:>16s} {:>16s} {:>18s}".format(
            "trial", "kind", "rank", "xi", "zeta", "reduced_residual"
        )
    )
    for i, (
        kind,
        rank,
        xi,
        zeta,
        reduced,
    ) in enumerate(
        zip(
            result.trial_kind_history,
            result.trial_basis_size_history,
            result.trial_indicator_history,
            result.saturation_history,
            result.trial_reduced_residual_history,
        ),
        start=1,
    ):
        print(
            "{:6d} {:>6s} {:8d} {:16.9e} {:16.9e} {:18.9e}".format(
                i,
                kind,
                int(rank),
                float(xi),
                float(zeta),
                float(reduced),
            )
        )

    print("-" * 112)
    print("Persistent commit history")
    print("-" * 112)
    print(
        "{:>8s} {:>8s} {:>12s} {:>16s}".format(
            "commit", "kind", "modes_added", "accepted_xi"
        )
    )
    for i, (kind, added, xi) in enumerate(
        zip(
            result.commit_kind_history,
            result.modes_added_history,
            result.indicator_history,
        ),
        start=1,
    ):
        print(
            "{:8d} {:>8s} {:12d} {:16.9e}".format(
                i,
                kind,
                int(added),
                float(xi),
            )
        )

    enrichment = result.last_enrichment_result
    if enrichment is not None:
        print("-" * 112)
        print("Last enrichment diagnostics")
        print("-" * 112)
        print("accepted               =", enrichment.accepted)
        print("failure_reason          =", enrichment.failure_reason)
        print("fixed-point converged   =", enrichment.fixed_point_converged)
        print("fixed-point iterations  =", enrichment.fixed_point_iterations)
        print("spatial novelty         = {:.9e}".format(enrichment.spatial_novelty))
        print(
            "temporal significance   = {:.9e}".format(
                enrichment.temporal_significance
            )
        )
        print(
            "orthogonality error     = {:.9e}".format(
                enrichment.orthogonality_error
            )
        )
        print(
            "residual before/after   = {:.9e} / {:.9e}".format(
                enrichment.residual_norm_before,
                enrichment.residual_norm_after,
            )
        )
        print(
            "residual benefit        = {:.9e}".format(
                enrichment.residual_benefit
            )
        )
        print(
            "fixed-point history     = {}".format(
                np.array2string(
                    enrichment.fixed_point_history,
                    precision=6,
                    suppress_small=False,
                    max_line_width=160,
                )
            )
        )

    if raw_fixed_point_records:
        record = raw_fixed_point_records[-1]
        modes = record["modes"]
        if len(modes) >= 4:
            print("-" * 112)
            print("Last raw fixed-point pair lag diagnostics")
            print("-" * 112)
            print(
                "{:>8s} {:>16s} {:>16s} {:>16s}".format(
                    "sweep",
                    "lag-1 distance",
                    "lag-2 distance",
                    "lag-3 distance",
                )
            )
            start_index = max(3, len(modes) - 12)
            for k in range(start_index, len(modes)):
                d1 = original_pair_change(
                    modes[k - 1],
                    modes[k],
                    record["time"],
                    record["H_sigma"],
                    record["metric"],
                )
                d2 = original_pair_change(
                    modes[k - 2],
                    modes[k],
                    record["time"],
                    record["H_sigma"],
                    record["metric"],
                )
                d3 = original_pair_change(
                    modes[k - 3],
                    modes[k],
                    record["time"],
                    record["H_sigma"],
                    record["metric"],
                )
                print(
                    "{:8d} {:16.9e} {:16.9e} {:16.9e}".format(
                        k,
                        d1,
                        d2,
                        d3,
                    )
                )

    print("-" * 112)
    print("Returned persistent state")
    print("-" * 112)
    print(
        "max |plastic strain|    = {:.9e}".format(
            float(np.max(np.abs(result.state.plastic_strain)))
        )
    )
    print(
        "max damage              = {:.9e}".format(
            float(np.max(result.state.damage))
        )
    )
    print(
        "max |stress|            = {:.9e} MPa".format(
            float(np.max(np.abs(result.state.stress)))
        )
    )
    print("=" * 112)


if __name__ == "__main__":
    main()

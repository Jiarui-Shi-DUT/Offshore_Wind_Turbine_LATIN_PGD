# -*- coding: utf-8 -*-
"""Diagnose usefulness of the three late phases of the failing fourth PGD mode.

This is a diagnostic-only example.  It does NOT relax the formal requirement
that a raw enrichment pair must converge before Gram-Schmidt and acceptance.

The script reproduces the common seven-commit tower LATIN-PGD baseline, rebuilds
the same failing Trial-A reduced problem, runs the original unrelaxed
Eq. (70)-(72) alternating map for 120 sweeps, and stores the last three raw
pairs of the already-observed period-3 orbit.

Each late raw pair is then *diagnostically* passed through the existing
post-fixed-point operations:
    M-weighted Gram-Schmidt,
    exact temporal-coordinate transformation,
    all-mode Eq. (58)-(59) temporal re-optimisation with full forcing,
    full mechanical residual benefit.

This answers a narrow question:
Is the fourth-mode period-3 orbit merely cycling among candidates that would
all be insignificant/useless after basis management, or does it contain a
genuinely useful new spatial direction whose fixed-point map simply fails to
settle?

No persistent LATIN-PGD source file is modified.
"""

from __future__ import annotations

import numpy as np

from examples.elastic_tapered_tower import (
    TowerConfiguration,
    create_tower_geometry,
)
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.pgd_basis import PGDModeTower
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_global_finishing import prepare_frozen_global_data
from latin.tower_initialization import (
    compute_tower_elastic_initialization,
)
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_pgd_time_update import (
    update_tower_pgd_time_functions,
)
from latin.tower_search_directions import (
    compute_tower_descent_search_directions,
)
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


N_RAW_SWEEPS = 120


def _build_problem():
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
    return material, operator, initialization


def _rebuild_failing_trial_a(material, operator, initialization):
    baseline = solve_tower_latin_pgd(
        initial_state=initialization.state,
        materials=material,
        metric=operator.metric,
        equilibrium_operator=operator,
        mode_significance_tolerance=0.0,
        acceptance_tolerance=0.0,
        max_iterations=20,
        max_fixed_point_iterations=120,
    )

    materials = (material,) * baseline.state.n_material_points
    local_state = solve_tower_local_stage(
        global_state=baseline.state,
        materials=materials,
    )
    directions = compute_tower_descent_search_directions(
        local_state=local_state,
        materials=materials,
        regularization=0.15,
    )
    frozen_data = prepare_frozen_global_data(
        baseline_state=baseline.state,
        local_state=local_state,
        directions=directions,
        equilibrium_operator=operator,
    )
    fixed_a = update_tower_pgd_time_functions(
        basis=baseline.basis,
        time=baseline.state.time,
        forcing=frozen_data.full_plastic_forcing,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        equilibrium_operator=operator,
        reduced_tolerance=1.0e-4,
        rcond=1.0e-12,
    )
    return baseline, directions, frozen_data, fixed_a


def _late_raw_modes(
    fixed_a,
    time,
    H_sigma,
    metric,
    operator,
    iteration_added,
):
    basis = fixed_a.basis
    defect = fixed_a.mechanical_residual

    p = enrichment_module._seed(
        defect,
        H_sigma,
        metric,
        1.0e-14,
    )
    s = operator.apply_spatial(p).stress
    lam, rate = enrichment_module._temporal_solve(
        p,
        s,
        time,
        defect,
        H_sigma,
        metric,
    )
    current = PGDModeTower(
        p,
        s,
        lam,
        rate,
        iteration_added,
    )

    modes = []
    chis = []

    for sweep in range(1, N_RAW_SWEEPS + 1):
        p, s = enrichment_module._spatial_solve(
            current.temporal_amplitude,
            current.temporal_rate,
            time,
            defect,
            H_sigma,
            metric,
            operator,
            1.0e-14,
        )
        lam, rate = enrichment_module._temporal_solve(
            p,
            s,
            time,
            defect,
            H_sigma,
            metric,
        )
        candidate = PGDModeTower(
            p,
            s,
            lam,
            rate,
            iteration_added,
        )
        chi = enrichment_module._pair_change(
            current,
            candidate,
            time,
            H_sigma,
            metric,
        )
        current = candidate
        if sweep > N_RAW_SWEEPS - 6:
            modes.append((sweep, candidate))
            chis.append((sweep, float(chi)))

    return modes, chis


def main() -> None:
    material, operator, initialization = _build_problem()
    baseline, directions, frozen_data, fixed_a = (
        _rebuild_failing_trial_a(
            material,
            operator,
            initialization,
        )
    )

    print("=" * 154)
    print("Tower PGD fourth-mode late-cycle usefulness diagnostic")
    print("=" * 154)
    print(
        "baseline: termination={}, committed={}, rank={}, xi={:.9e}".format(
            baseline.termination_reason.value,
            baseline.iterations,
            baseline.basis.n_modes,
            baseline.final_indicator,
        )
    )
    print(
        "failing Trial-A relative residual = {:.9e}".format(
            fixed_a.relative_residual
        )
    )

    modes, chis = _late_raw_modes(
        fixed_a=fixed_a,
        time=baseline.state.time,
        H_sigma=directions.H_sigma,
        metric=operator.metric,
        operator=operator,
        iteration_added=baseline.attempted_iterations,
    )

    print("-" * 154)
    print("Last six unrelaxed fixed-point sweep changes")
    print("-" * 154)
    for sweep, chi in chis:
        print("sweep {:3d}: chi = {:.9e}".format(sweep, chi))

    # The last three modes are one representative of each phase of the
    # established period-3 orbit.
    representatives = modes[-3:]

    before = enrichment_module._be_residual_norm(
        fixed_a.mechanical_residual,
        baseline.state.time,
        directions.H_sigma,
        operator.metric,
    )

    print("-" * 154)
    print("Diagnostic post-processing of the three late period-3 phases")
    print("-" * 154)
    print(
        "{:>7s} {:>12s} {:>12s} {:>12s} {:>13s} {:>14s} {:>14s} {:>16s}".format(
            "sweep",
            "gamma_sp",
            "orth_scale",
            "gamma_lambda",
            "ortho_error",
            "resid_after",
            "resid_benefit",
            "relative_resid",
        )
    )

    for sweep, raw_mode in representatives:
        try:
            (
                tentative,
                _coeff,
                c,
                gamma_sp,
                gamma_lambda,
                ortho_err,
                _ep_err,
                _rate_err,
                _stress_err,
            ) = enrichment_module._post_fixed_point_transform(
                fixed_a.basis,
                raw_mode,
                baseline.state.time,
                operator.metric,
                operator,
                1.0e-14,
                1.0e-12,
                0.0,
                1.0e-8,
                1.0e-10,
                2,
            )

            candidate = update_tower_pgd_time_functions(
                basis=tentative,
                time=baseline.state.time,
                forcing=frozen_data.full_plastic_forcing,
                H_sigma=directions.H_sigma,
                metric=operator.metric,
                equilibrium_operator=operator,
                reduced_tolerance=1.0e-4,
                rcond=1.0e-12,
            )

            after = enrichment_module._be_residual_norm(
                candidate.mechanical_residual,
                baseline.state.time,
                directions.H_sigma,
                operator.metric,
            )
            benefit = (
                0.0
                if before <= np.finfo(np.float64).eps
                else 1.0 - after / before
            )

            print(
                "{:7d} {:12.5e} {:12.5e} {:12.5e} {:13.5e} "
                "{:14.7e} {:14.7e} {:16.7e}".format(
                    sweep,
                    gamma_sp,
                    c,
                    gamma_lambda,
                    ortho_err,
                    after,
                    benefit,
                    candidate.relative_residual,
                )
            )
        except Exception as exc:
            print(
                "{:7d} diagnostic post-processing failed: {}: {}".format(
                    sweep,
                    type(exc).__name__,
                    exc,
                )
            )

    print("-" * 154)
    print(
        "Important: these three raw pairs are NOT formally accepted modes; "
        "they did not satisfy the fixed-point convergence gate."
    )
    print(
        "The table is diagnostic only and is used to decide whether the "
        "periodic orbit is cycling among useless or potentially useful "
        "fourth-mode candidates."
    )
    print("=" * 154)


if __name__ == "__main__":
    main()

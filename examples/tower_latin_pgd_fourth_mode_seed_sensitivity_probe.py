# -*- coding: utf-8 -*-
"""Seed-sensitivity diagnostic for the failing fourth tower PGD enrichment.

The persistent LATIN-PGD implementation is not modified.  The script reproduces
the common seven-commit baseline, rebuilds the same failing Trial-A reduced
problem, and reruns the ORIGINAL unrelaxed Eq. (70)-(72) raw fixed-point map
from several deterministic residual-row seeds.

The tested seeds are the highest-energy time rows of the same shifted defect
used by the current residual-driven seed strategy.  This isolates whether the
observed period-3 orbit is primarily a basin-of-attraction artefact of the
single argmax seed or a more intrinsic feature of the failing fourth-mode map.

For every seed the script reports:
    convergence,
    final complete-pair change,
    lag-3 complete-pair distance,
    diagnostic post-Gram-Schmidt spatial/temporal significance,
    diagnostic full residual benefit.

Passing an unconverged final pair through basis management is diagnostic only
and does not alter the formal fixed-point convergence gate.
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
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_search_directions import (
    compute_tower_descent_search_directions,
)
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


N_SEEDS = 10
MAX_FIXED_POINT_ITERATIONS = 200
FIXED_POINT_TOLERANCE = 1.0e-6


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


def _seed_row_energy(defect, H_sigma, metric):
    return np.sum(
        defect**2
        * metric.weights[None, :]
        / H_sigma,
        axis=1,
    )


def _run_raw_from_row_seed(
    seed_row,
    fixed_a,
    time,
    H_sigma,
    metric,
    operator,
    iteration_added,
):
    defect = fixed_a.mechanical_residual

    p = -np.asarray(defect[seed_row, :], dtype=np.float64).copy()
    p_norm = metric.norm(p)
    if (
        not np.isfinite(p_norm)
        or p_norm <= 1.0e-14
    ):
        raise RuntimeError("selected residual-row seed is degenerate")
    p /= p_norm

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

    recent_modes = [current]
    history = []
    converged = False

    for _ in range(MAX_FIXED_POINT_ITERATIONS):
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
        history.append(float(chi))
        current = candidate
        recent_modes.append(current)
        if len(recent_modes) > 4:
            recent_modes.pop(0)

        if chi <= FIXED_POINT_TOLERANCE:
            converged = True
            break

    lag3 = float("nan")
    if len(recent_modes) == 4:
        lag3 = enrichment_module._pair_change(
            recent_modes[0],
            recent_modes[3],
            time,
            H_sigma,
            metric,
        )

    return (
        current,
        np.asarray(history, dtype=np.float64),
        converged,
        lag3,
    )


def _diagnostic_usefulness(
    raw_mode,
    fixed_a,
    time,
    H_sigma,
    metric,
    operator,
    full_forcing,
):
    before = enrichment_module._be_residual_norm(
        fixed_a.mechanical_residual,
        time,
        H_sigma,
        metric,
    )

    (
        tentative,
        _coeff,
        _c,
        gamma_sp,
        gamma_lambda,
        _ortho_err,
        _ep_err,
        _rate_err,
        _stress_err,
    ) = enrichment_module._post_fixed_point_transform(
        fixed_a.basis,
        raw_mode,
        time,
        metric,
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
        time=time,
        forcing=full_forcing,
        H_sigma=H_sigma,
        metric=metric,
        equilibrium_operator=operator,
        reduced_tolerance=1.0e-4,
        rcond=1.0e-12,
    )
    after = enrichment_module._be_residual_norm(
        candidate.mechanical_residual,
        time,
        H_sigma,
        metric,
    )
    benefit = (
        0.0
        if before <= np.finfo(np.float64).eps
        else 1.0 - after / before
    )
    return (
        float(gamma_sp),
        float(gamma_lambda),
        float(benefit),
        float(candidate.relative_residual),
    )


def main() -> None:
    material, operator, initialization = _build_problem()
    baseline, directions, frozen_data, fixed_a = (
        _rebuild_failing_trial_a(
            material,
            operator,
            initialization,
        )
    )

    energy = _seed_row_energy(
        fixed_a.mechanical_residual,
        directions.H_sigma,
        operator.metric,
    )
    order = np.argsort(energy)[::-1]
    selected = order[: min(N_SEEDS, order.size)]
    max_energy = float(energy[selected[0]])

    print("=" * 178)
    print("Tower PGD fourth-mode deterministic residual-row seed sensitivity")
    print("=" * 178)
    print(
        "baseline: termination={}, committed={}, rank={}, xi={:.9e}".format(
            baseline.termination_reason.value,
            baseline.iterations,
            baseline.basis.n_modes,
            baseline.final_indicator,
        )
    )
    print(
        "failing Trial-A relative residual = {:.9e}; "
        "testing top {} residual-energy rows".format(
            fixed_a.relative_residual,
            selected.size,
        )
    )
    print("-" * 178)
    print(
        "{:>5s} {:>9s} {:>11s} {:>9s} {:>14s} {:>14s} "
        "{:>12s} {:>12s} {:>14s} {:>16s}".format(
            "rank",
            "time_idx",
            "energy/max",
            "conv",
            "last_chi",
            "lag3",
            "gamma_sp",
            "gamma_lam",
            "resid_benefit",
            "relative_resid",
        )
    )

    for rank, seed_row in enumerate(selected, start=1):
        try:
            raw_mode, history, converged, lag3 = (
                _run_raw_from_row_seed(
                    seed_row=int(seed_row),
                    fixed_a=fixed_a,
                    time=baseline.state.time,
                    H_sigma=directions.H_sigma,
                    metric=operator.metric,
                    operator=operator,
                    iteration_added=baseline.attempted_iterations,
                )
            )
            (
                gamma_sp,
                gamma_lambda,
                benefit,
                relative_residual,
            ) = _diagnostic_usefulness(
                raw_mode=raw_mode,
                fixed_a=fixed_a,
                time=baseline.state.time,
                H_sigma=directions.H_sigma,
                metric=operator.metric,
                operator=operator,
                full_forcing=frozen_data.full_plastic_forcing,
            )

            last_chi = (
                float(history[-1])
                if history.size
                else float("nan")
            )
            print(
                "{:5d} {:9d} {:11.5e} {:>9s} {:14.7e} {:14.7e} "
                "{:12.5e} {:12.5e} {:14.7e} {:16.7e}".format(
                    rank,
                    int(seed_row),
                    float(energy[seed_row] / max_energy),
                    str(converged),
                    last_chi,
                    lag3,
                    gamma_sp,
                    gamma_lambda,
                    benefit,
                    relative_residual,
                )
            )
        except Exception as exc:
            print(
                "{:5d} {:9d} {:11.5e} FAILED: {}: {}".format(
                    rank,
                    int(seed_row),
                    float(energy[seed_row] / max_energy),
                    type(exc).__name__,
                    exc,
                )
            )

    print("-" * 178)
    print(
        "Interpretation rule: if alternative high-energy residual-row seeds "
        "converge, the failure is seed/basin sensitive.  If they all remain "
        "nonconvergent and show small lag-3 distances, the period-3 attractor "
        "is robust to this deterministic restart family."
    )
    print(
        "Post-basis-management metrics for unconverged raw pairs are "
        "diagnostic only; they do not bypass the formal fixed-point gate."
    )
    print("=" * 178)


if __name__ == "__main__":
    main()

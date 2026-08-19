# -*- coding: utf-8 -*-
"""Sweep spatial under-relaxation for the exact failing fourth PGD mode.

The script first reproduces the original tower LATIN-PGD solve until the
fourth enrichment fails.  The returned persistent state/basis are therefore
the common seven-commit baseline.  It then rebuilds the same Trial-A reduced
problem once and tests several spatial relaxation factors on exactly that same
shifted defect.

For every relaxed iterate, the script also evaluates the defect of the
ORIGINAL unrelaxed alternating map.  This prevents a small relaxation factor
from appearing converged merely because the damped step length is small.

No persistent LATIN-PGD source file is modified.
"""

from __future__ import annotations

from time import perf_counter

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration, create_tower_geometry
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import ViscoplasticDamageTowerSystem2D
from latin.pgd_basis import PGDModeTower
from latin.tower_equilibrium_operator import build_tower_equilibrium_operator
from latin.tower_global_finishing import prepare_frozen_global_data
from latin.tower_initialization import compute_tower_elastic_initialization
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_search_directions import compute_tower_descent_search_directions
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


OMEGAS = (1.0, 0.75, 0.5, 0.35, 0.25, 0.20, 0.15, 0.10, 0.05)
MAX_FIXED_POINT_ITERATIONS = 400


def _make_relaxed_raw_fixed_point(
    omega: float,
    raw_history_sink: list[float],
):
    def relaxed_raw_fixed_point(
        basis,
        time,
        defect,
        H_sigma,
        metric,
        operator,
        iteration_added,
        fixed_point_tolerance,
        max_fixed_point_iterations,
        minimum_spatial_norm,
    ):
        p = enrichment_module._seed(
            defect,
            H_sigma,
            metric,
            minimum_spatial_norm,
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

        history = []
        converged = False

        for _ in range(max_fixed_point_iterations):
            p_star, s_star = enrichment_module._spatial_solve(
                current.temporal_amplitude,
                current.temporal_rate,
                time,
                defect,
                H_sigma,
                metric,
                operator,
                minimum_spatial_norm,
            )

            if metric.inner_product(
                current.spatial_plastic_strain,
                p_star,
            ) < 0.0:
                p_star = -p_star
                s_star = -s_star

            # Measure the defect of the ORIGINAL unrelaxed alternating map
            # F(z)-z at the current relaxed iterate.  This is the quantity
            # that must vanish at a true fixed point.  The relaxed step
            # distance alone scales approximately with omega and can become
            # artificially small as omega -> 0.
            lam_star, rate_star = enrichment_module._temporal_solve(
                p_star,
                s_star,
                time,
                defect,
                H_sigma,
                metric,
            )
            raw_candidate = PGDModeTower(
                p_star,
                s_star,
                lam_star,
                rate_star,
                iteration_added,
            )
            raw_chi = enrichment_module._pair_change(
                current,
                raw_candidate,
                time,
                H_sigma,
                metric,
            )
            raw_history_sink.append(float(raw_chi))

            p_mix = (
                (1.0 - omega) * current.spatial_plastic_strain
                + omega * p_star
            )
            p_norm = metric.norm(p_mix)
            if (
                not np.isfinite(p_norm)
                or p_norm <= minimum_spatial_norm
            ):
                raise enrichment_module._CandidateFailure(
                    "relaxed_spatial_mode_degenerate"
                )

            p = p_mix / p_norm
            s = operator.apply_spatial(p).stress
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
            history.append(chi)
            current = candidate

            # Convergence is diagnosed with the unrelaxed-map defect,
            # not with the damped step length.
            if raw_chi <= fixed_point_tolerance:
                converged = True
                break

        return (
            current,
            np.asarray(history, dtype=np.float64),
            converged,
        )

    return relaxed_raw_fixed_point


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

    print("=" * 132)
    print("Tower PGD fourth-mode spatial-relaxation sweep")
    print("=" * 132)
    print("Reproducing common seven-commit persistent baseline...")

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

    print(
        "baseline termination={}, committed={}, rank={}, xi={:.9e}".format(
            baseline.termination_reason.value,
            baseline.iterations,
            baseline.basis.n_modes,
            baseline.final_indicator,
        )
    )

    local_state = solve_tower_local_stage(
        global_state=baseline.state,
        materials=(material,) * baseline.state.n_material_points,
    )
    directions = compute_tower_descent_search_directions(
        local_state=local_state,
        materials=(material,) * baseline.state.n_material_points,
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

    print(
        "common failing Trial-A reduced residual = {:.9e}".format(
            fixed_a.relative_residual
        )
    )
    print("-" * 132)
    print(
        "{:>8s} {:>11s} {:>9s} {:>14s} {:>14s} {:>14s} {:>14s} {:>20s}".format(
            "omega",
            "converged",
            "fp_iters",
            "relaxed_last",
            "raw_last",
            "raw_min_tail",
            "raw_mean_tail",
            "failure_reason",
        )
    )

    original_raw_fixed_point = enrichment_module._raw_fixed_point

    for omega in OMEGAS:
        raw_history: list[float] = []
        enrichment_module._raw_fixed_point = _make_relaxed_raw_fixed_point(
            float(omega),
            raw_history,
        )
        start = perf_counter()
        try:
            result = enrichment_module.enrich_tower_pgd_basis_once(
                fixed_basis_result=fixed_a,
                time=baseline.state.time,
                full_forcing=frozen_data.full_plastic_forcing,
                shifted_defect=fixed_a.mechanical_residual,
                H_sigma=directions.H_sigma,
                metric=operator.metric,
                equilibrium_operator=operator,
                mode_significance_tolerance=0.0,
                acceptance_tolerance=0.0,
                iteration_added=baseline.attempted_iterations,
                fixed_point_tolerance=1.0e-6,
                max_fixed_point_iterations=MAX_FIXED_POINT_ITERATIONS,
                minimum_spatial_norm=1.0e-14,
                minimum_spatial_novelty=1.0e-12,
                basis_health_tolerance=1.0e-8,
                field_invariance_tolerance=1.0e-10,
                reorthogonalization_passes=2,
                reduced_tolerance=1.0e-4,
                rcond=1.0e-12,
            )
        finally:
            enrichment_module._raw_fixed_point = original_raw_fixed_point
        elapsed = perf_counter() - start

        hist = result.fixed_point_history
        relaxed_last = float(hist[-1]) if hist.size else float("nan")
        raw = np.asarray(raw_history, dtype=np.float64)
        raw_tail = raw[-min(50, raw.size):]
        raw_last = float(raw[-1]) if raw.size else float("nan")
        raw_min_tail = float(np.min(raw_tail)) if raw_tail.size else float("nan")
        raw_mean_tail = float(np.mean(raw_tail)) if raw_tail.size else float("nan")

        print(
            "{:8.3f} {:>11s} {:9d} {:14.7e} {:14.7e} {:14.7e} {:14.7e} {:>20s}  ({:.3f}s)".format(
                float(omega),
                str(result.fixed_point_converged),
                int(result.fixed_point_iterations),
                relaxed_last,
                raw_last,
                raw_min_tail,
                raw_mean_tail,
                str(result.failure_reason),
                elapsed,
            )
        )

        if result.fixed_point_converged:
            print(
                "         gamma_sp={:.6e}, gamma_lambda={:.6e}, "
                "residual_benefit={:.6e}".format(
                    result.spatial_novelty,
                    result.temporal_significance,
                    result.residual_benefit,
                )
            )

    print("=" * 132)


if __name__ == "__main__":
    main()

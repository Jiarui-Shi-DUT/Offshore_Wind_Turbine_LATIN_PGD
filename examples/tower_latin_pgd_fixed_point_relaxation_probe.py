# -*- coding: utf-8 -*-
"""Diagnostic spatial under-relaxation for tower PGD enrichment fixed points.

This file does not modify the persistent LATIN-PGD implementation.  During this
single process it temporarily replaces the raw one-mode alternating fixed-point
map by a spatially under-relaxed version with omega = 0.5.

The relaxed spatial update is

    p_mix = (1 - omega) p_old + omega p_star,

after metric-sign alignment, followed by M-normalisation, equilibrium stress
reconstruction, and an exact temporal half-step.  Thus every iterate remains
one separated rank-one pair and every fixed point of the original spatial map
remains a fixed point of the relaxed iteration.
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
from latin.tower_initialization import compute_tower_elastic_initialization
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


SPATIAL_RELAXATION = 0.25


def _relaxed_raw_fixed_point(
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
    """Run the tower raw fixed point with spatial under-relaxation only."""
    omega = float(SPATIAL_RELAXATION)

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
        p_star, _ = enrichment_module._spatial_solve(
            current.temporal_amplitude,
            current.temporal_rate,
            time,
            defect,
            H_sigma,
            metric,
            operator,
            minimum_spatial_norm,
        )

        # The separated representation has a sign gauge.  Align the new
        # spatial direction before convex mixing so opposite but equivalent
        # representations are not averaged toward zero.
        overlap = metric.inner_product(
            current.spatial_plastic_strain,
            p_star,
        )
        if overlap < 0.0:
            p_star = -p_star

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

        if chi <= fixed_point_tolerance:
            converged = True
            break

    return (
        current,
        np.asarray(history, dtype=np.float64),
        converged,
    )


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
    print("Tower LATIN-PGD fixed-point spatial under-relaxation probe")
    print("=" * 112)
    print(
        "omega={}, Nt={}, Nq={}, max_iterations=20, "
        "max_fixed_point_iterations=240".format(
            SPATIAL_RELAXATION,
            initialization.state.n_time,
            initialization.state.n_material_points,
        )
    )

    original_raw_fixed_point = enrichment_module._raw_fixed_point
    enrichment_module._raw_fixed_point = _relaxed_raw_fixed_point
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
            max_fixed_point_iterations=240,
        )
    finally:
        enrichment_module._raw_fixed_point = original_raw_fixed_point
    elapsed = perf_counter() - start

    print("-" * 112)
    print("Solver summary")
    print("-" * 112)
    print("termination_reason   =", result.termination_reason.value)
    print("converged            =", result.converged)
    print("failure_reason       =", result.failure_reason)
    print("attempted iterations =", result.attempted_iterations)
    print("committed iterations =", result.iterations)
    print("final basis modes    =", result.basis.n_modes)
    print("total modes added    =", result.total_modes_added)
    print("final accepted xi    = {:.9e}".format(result.final_indicator))
    print("elapsed              = {:.6f} s".format(elapsed))
    print("commit kinds         =", result.commit_kind_history)
    print(
        "indicator history    =",
        np.array2string(
            result.indicator_history,
            precision=6,
            suppress_small=False,
            max_line_width=160,
        ),
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
        print(
            "spatial novelty         = {:.9e}".format(
                enrichment.spatial_novelty
            )
        )
        print(
            "temporal significance   = {:.9e}".format(
                enrichment.temporal_significance
            )
        )
        print(
            "residual benefit        = {:.9e}".format(
                enrichment.residual_benefit
            )
        )
        history = enrichment.fixed_point_history
        print(
            "fixed-point history tail =",
            np.array2string(
                history[-min(30, history.size):],
                precision=6,
                suppress_small=False,
                max_line_width=160,
            ),
        )

    print("-" * 112)
    print("Returned persistent state")
    print("-" * 112)
    print(
        "max |plastic strain| = {:.9e}".format(
            float(np.max(np.abs(result.state.plastic_strain)))
        )
    )
    print(
        "max damage           = {:.9e}".format(
            float(np.max(result.state.damage))
        )
    )
    print("=" * 112)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""Diagnose the failing fourth tower PGD mode with a whole-time temporal solve.

The persistent LATIN-PGD implementation is not modified.  The script first
reproduces the common seven-commit baseline with the current implementation,
rebuilds the same failing Trial-A reduced problem, and then temporarily replaces
only the single-new-mode temporal half-step inside tower PGD enrichment.

The diagnostic temporal half-step keeps the current tower-v1 backward-Euler
residual and lambda(0)=0 convention, but minimises the complete discrete
space-time residual jointly over lambda_1,...,lambda_N rather than performing a
causal sequence of one-step least-squares solves.

This is NOT claimed to reproduce the paper's exact DG0 algebra.  It isolates
whether the sequential temporal optimisation scope is responsible for the
observed fourth-mode fixed-point cycle.
"""

from __future__ import annotations

from time import perf_counter

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration, create_tower_geometry
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import ViscoplasticDamageTowerSystem2D
from latin.tower_equilibrium_operator import build_tower_equilibrium_operator
from latin.tower_global_finishing import prepare_frozen_global_data
from latin.tower_initialization import compute_tower_elastic_initialization
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_search_directions import compute_tower_descent_search_directions
import latin.tower_pgd_enrichment as enrichment_module
from material.viscoplastic_damage_1d import MaterialParameters


MAX_FIXED_POINT_ITERATIONS = 200


def _whole_time_temporal_solve(
    p,
    s,
    time,
    defect,
    H_sigma,
    metric,
):
    """Jointly minimise the BE-discretised Eq. (72) residual over all times."""
    t = np.asarray(time, dtype=np.float64)
    d = np.asarray(defect, dtype=np.float64)
    Hs = np.asarray(H_sigma, dtype=np.float64)

    nt = t.size
    nq = p.size
    n_unknown = nt - 1

    lam = np.zeros(nt, dtype=np.float64)
    rate = np.zeros(nt, dtype=np.float64)

    # Keep the same current tower-v1 t0 rate treatment.  The fixed-point graph
    # norm and the spatial half-step use the positive-time slabs.
    w0 = metric.weights / Hs[0, :]
    den0 = float(np.dot(w0 * p, p))
    if den0 <= np.finfo(np.float64).tiny:
        raise enrichment_module._CandidateFailure(
            "whole_time_temporal_initial_denominator_degenerate"
        )
    rate[0] = -float(np.dot(w0 * p, d[0, :])) / den0

    if n_unknown == 0:
        return lam, rate

    # Each time slab contributes Nq weighted residual equations:
    #
    #   r_n = (p/dt - H_n s) lambda_n
    #         - (p/dt) lambda_{n-1}
    #         + defect_n .
    #
    # lambda_0 is fixed to zero.  Stack all slabs and solve one global
    # weighted least-squares problem for lambda_1,...,lambda_N.
    matrix = np.zeros(
        (n_unknown * nq, n_unknown),
        dtype=np.float64,
    )
    rhs = np.zeros(n_unknown * nq, dtype=np.float64)

    for n in range(1, nt):
        dt = float(t[n] - t[n - 1])
        Hn = Hs[n, :]
        sqrt_weights = np.sqrt(
            dt * metric.weights / Hn
        )

        g = p / dt - Hn * s
        h = -p / dt

        row0 = (n - 1) * nq
        row1 = n * nq

        matrix[row0:row1, n - 1] = sqrt_weights * g
        if n > 1:
            matrix[row0:row1, n - 2] = sqrt_weights * h

        rhs[row0:row1] = -sqrt_weights * d[n, :]

    solution, _, _, singular_values = np.linalg.lstsq(
        matrix,
        rhs,
        rcond=1.0e-12,
    )
    if np.any(~np.isfinite(solution)):
        raise enrichment_module._CandidateFailure(
            "whole_time_temporal_solution_non_finite"
        )

    lam[1:] = solution
    rate[1:] = np.diff(lam) / np.diff(t)

    if singular_values.size:
        smallest = float(singular_values[-1])
        largest = float(singular_values[0])
        condition = (
            np.inf
            if smallest <= np.finfo(np.float64).eps
            else largest / smallest
        )
    else:
        condition = 0.0

    _whole_time_temporal_solve.last_condition = float(condition)
    return lam, rate


_whole_time_temporal_solve.last_condition = float("nan")


def main() -> None:
    configuration = TowerConfiguration(
        horizontal_force=1.0e6,
        n_elements=10,
        n_gauss=2,
        n_circumferential=16,
        n_radial=1,
    )
    material = MaterialParameters()
    materials = None

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

    print("=" * 118)
    print("Tower PGD fourth-mode whole-time temporal-minimisation probe")
    print("=" * 118)
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

    print(
        "common failing Trial-A reduced residual = {:.9e}".format(
            fixed_a.relative_residual
        )
    )
    print(
        "diagnostic: whole-time BE least squares for NEW mode only; "
        "max_fixed_point_iterations={}".format(
            MAX_FIXED_POINT_ITERATIONS
        )
    )

    original_temporal_solve = enrichment_module._temporal_solve
    enrichment_module._temporal_solve = _whole_time_temporal_solve
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
        enrichment_module._temporal_solve = original_temporal_solve
    elapsed = perf_counter() - start

    print("-" * 118)
    print("Enrichment result")
    print("-" * 118)
    print("accepted               =", result.accepted)
    print("failure_reason          =", result.failure_reason)
    print("fixed-point converged   =", result.fixed_point_converged)
    print("fixed-point iterations  =", result.fixed_point_iterations)
    print(
        "last temporal condition = {:.9e}".format(
            _whole_time_temporal_solve.last_condition
        )
    )
    history = result.fixed_point_history
    print(
        "fixed-point history tail =",
        np.array2string(
            history[-min(40, history.size):],
            precision=6,
            suppress_small=False,
            max_line_width=170,
        ),
    )
    print(
        "spatial novelty         = {:.9e}".format(
            result.spatial_novelty
        )
    )
    print(
        "temporal significance   = {:.9e}".format(
            result.temporal_significance
        )
    )
    print(
        "orthogonality error     = {:.9e}".format(
            result.orthogonality_error
        )
    )
    print(
        "residual before/after   = {:.9e} / {:.9e}".format(
            result.residual_norm_before,
            result.residual_norm_after,
        )
    )
    print(
        "residual benefit        = {:.9e}".format(
            result.residual_benefit
        )
    )
    print("elapsed                 = {:.6f} s".format(elapsed))
    print("=" * 118)


if __name__ == "__main__":
    main()

# -*- coding: utf-8 -*-
"""I-4B production-path probe on the known failing fourth tower PGD mode.

This probe compares the two spatial strategies through the real
enrich_tower_pgd_basis_once() transaction:

1) paper_galerkin  -> expected fourth-mode fixed-point failure;
2) residual_ls     -> expected fixed-point convergence and accepted rank growth.

It also verifies that the provisional Trial-A object is not mutated.
"""

from __future__ import annotations

import numpy as np

from examples.tower_latin_pgd_fourth_mode_direct_ls_spatial_probe import (
    _build_problem,
    _rebuild_failing_trial_a,
)
from latin.tower_pgd_enrichment import enrich_tower_pgd_basis_once


FIXED_POINT_TOLERANCE = 1.0e-6
MAX_FIXED_POINT_ITERATIONS = 200


def _equilibrium_relative(mode, operator) -> tuple[float, float]:
    residual = operator.equilibrium_residual(mode.spatial_stress)
    absolute = float(np.linalg.norm(residual))
    scale_vector = (
        np.abs(operator.compatibility_matrix).T
        @ (
            operator.metric.weights
            * np.abs(mode.spatial_stress)
        )
    )
    scale = float(np.linalg.norm(scale_vector))
    relative = absolute / max(
        scale,
        np.finfo(np.float64).eps,
    )
    return absolute, float(relative)


def main() -> None:
    material, operator, initialization = _build_problem()
    baseline, directions, fixed_a = _rebuild_failing_trial_a(
        material,
        operator,
        initialization,
    )

    time = baseline.state.time
    H_sigma = directions.H_sigma
    defect = fixed_a.mechanical_residual
    basis = fixed_a.basis

    # From the fixed-basis residual definition
    #
    #   r = P lambda_dot - H_sigma S lambda - f,
    #
    # recover exactly the frozen full forcing f used to construct Trial-A.
    full_forcing = (
        fixed_a.plastic_strain_rate_correction
        - H_sigma * basis.stress_correction()
        - defect
    )

    before_basis = fixed_a.basis.copy()
    before_residual = fixed_a.mechanical_residual.copy()
    before_rank = fixed_a.basis.n_modes

    common = dict(
        fixed_basis_result=fixed_a,
        time=time,
        full_forcing=full_forcing,
        shifted_defect=defect,
        H_sigma=H_sigma,
        metric=operator.metric,
        equilibrium_operator=operator,
        mode_significance_tolerance=0.0,
        acceptance_tolerance=0.0,
        iteration_added=baseline.attempted_iterations,
        fixed_point_tolerance=FIXED_POINT_TOLERANCE,
        max_fixed_point_iterations=MAX_FIXED_POINT_ITERATIONS,
    )

    print("=" * 110)
    print("I-4B production enrichment: known fourth-mode benchmark")
    print("=" * 110)
    print(
        "baseline: termination={}, committed={}, attempted={}, rank={}, xi={:.12e}".format(
            baseline.termination_reason.value,
            baseline.iterations,
            baseline.attempted_iterations,
            baseline.basis.n_modes,
            baseline.final_indicator,
        )
    )
    print(
        "Trial-A: relative residual={:.12e}, Nq={}, Nt={}".format(
            fixed_a.relative_residual,
            operator.n_material_points,
            time.size,
        )
    )

    print("-" * 110)
    print("A) Existing paper_galerkin production path")
    paper = enrich_tower_pgd_basis_once(
        **common,
        spatial_strategy="paper_galerkin",
    )
    print("  accepted                 =", paper.accepted)
    print("  failure reason           =", paper.failure_reason)
    print("  fixed-point converged    =", paper.fixed_point_converged)
    print("  fixed-point iterations   =", paper.fixed_point_iterations)
    if paper.fixed_point_history.size:
        print(
            "  final chi                = {:.12e}".format(
                float(paper.fixed_point_history[-1])
            )
        )

    print("-" * 110)
    print("B) New residual_ls production path")
    residual_ls = enrich_tower_pgd_basis_once(
        **common,
        spatial_strategy="residual_ls",
    )
    print("  accepted                 =", residual_ls.accepted)
    print("  failure reason           =", residual_ls.failure_reason)
    print("  fixed-point converged    =", residual_ls.fixed_point_converged)
    print("  fixed-point iterations   =", residual_ls.fixed_point_iterations)
    if residual_ls.fixed_point_history.size:
        print(
            "  final chi                = {:.12e}".format(
                float(residual_ls.fixed_point_history[-1])
            )
        )
    print("  rank before              =", before_rank)
    print("  rank after               =", residual_ls.n_modes)
    print(
        "  residual before          = {:.12e}".format(
            residual_ls.residual_norm_before
        )
    )
    print(
        "  residual after           = {:.12e}".format(
            residual_ls.residual_norm_after
        )
    )
    print(
        "  full residual benefit    = {:.6%}".format(
            residual_ls.residual_benefit
        )
    )
    print(
        "  spatial novelty          = {:.12e}".format(
            residual_ls.spatial_novelty
        )
    )
    print(
        "  temporal significance    = {:.12e}".format(
            residual_ls.temporal_significance
        )
    )
    print(
        "  orthogonality error      = {:.12e}".format(
            residual_ls.orthogonality_error
        )
    )

    if residual_ls.raw_mode is not None:
        eq_abs, eq_rel = _equilibrium_relative(
            residual_ls.raw_mode,
            operator,
        )
        print(
            "  raw-mode equilibrium abs = {:.12e}".format(
                eq_abs
            )
        )
        print(
            "  raw-mode equilibrium rel = {:.12e}".format(
                eq_rel
            )
        )

    basis_unchanged = (
        fixed_a.basis.n_modes == before_basis.n_modes
        and np.array_equal(
            fixed_a.basis.plastic_strain_correction(),
            before_basis.plastic_strain_correction(),
        )
        and np.array_equal(
            fixed_a.basis.stress_correction(),
            before_basis.stress_correction(),
        )
    )
    residual_unchanged = np.array_equal(
        fixed_a.mechanical_residual,
        before_residual,
    )

    print("-" * 110)
    print("Transaction audit")
    print("  Trial-A basis unchanged  =", basis_unchanged)
    print("  Trial-A residual unchanged =", residual_unchanged)

    strong_pass = (
        (not paper.accepted)
        and paper.failure_reason == "fixed_point_not_converged"
        and residual_ls.accepted
        and residual_ls.fixed_point_converged
        and residual_ls.n_modes == before_rank + 1
        and residual_ls.residual_benefit > 0.0
        and residual_ls.orthogonality_error < 1.0e-8
        and basis_unchanged
        and residual_unchanged
    )

    print("-" * 110)
    print("I-4B STRONG PASS =", strong_pass)
    print("=" * 110)

    if not strong_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()

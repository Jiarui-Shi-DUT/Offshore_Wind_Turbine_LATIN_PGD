# -*- coding: utf-8 -*-
"""I-5 outer LATIN convergence probe using residual-LS PGD enrichment.

This diagnostic runs the known nonlinear fully reversed tower benchmark through
the complete transactional solve_tower_latin_pgd() path with

    spatial_strategy="residual_ls".

It does not change production code.  The purpose is to determine whether the
I-4B fourth-mode recovery persists when the outer LATIN iterations continue
beyond that enrichment event.
"""

from __future__ import annotations

import numpy as np

from examples.tower_latin_pgd_fourth_mode_direct_ls_spatial_probe import (
    _build_problem,
)
from latin.tower_latin_pgd_solver import solve_tower_latin_pgd


def _history_text(values, fmt=".6e") -> str:
    arr = np.asarray(values)
    if arr.size == 0:
        return "<empty>"
    return " ".join(format(float(x), fmt) for x in arr)


def main() -> None:
    material, operator, initialization = _build_problem()

    result = solve_tower_latin_pgd(
        initial_state=initialization.state,
        materials=material,
        metric=operator.metric,
        equilibrium_operator=operator,
        mode_significance_tolerance=0.0,
        acceptance_tolerance=0.0,
        max_iterations=30,
        max_fixed_point_iterations=200,
        spatial_strategy="residual_ls",
    )

    indicators = np.asarray(result.indicator_history, dtype=np.float64)
    baselines = np.asarray(
        result.baseline_indicator_history,
        dtype=np.float64,
    )
    trials = np.asarray(result.trial_indicator_history, dtype=np.float64)
    reduced = np.asarray(
        result.trial_reduced_residual_history,
        dtype=np.float64,
    )
    modes_added = np.asarray(result.modes_added_history, dtype=np.int64)

    if indicators.size >= 2:
        decreases = int(np.sum(np.diff(indicators) < 0.0))
        increases = int(np.sum(np.diff(indicators) > 0.0))
        unchanged = int(indicators.size - 1 - decreases - increases)
        overall_ratio = float(indicators[-1] / indicators[0])
        best = float(np.min(indicators))
    elif indicators.size == 1:
        decreases = 0
        increases = 0
        unchanged = 0
        overall_ratio = 1.0
        best = float(indicators[0])
    else:
        decreases = 0
        increases = 0
        unchanged = 0
        overall_ratio = float("nan")
        best = float("nan")

    print("=" * 116)
    print("I-5 outer LATIN convergence probe: residual-LS enrichment")
    print("=" * 116)
    print("termination reason          =", result.termination_reason.value)
    print("solver converged            =", result.converged)
    print("failure reason              =", result.failure_reason)
    print("committed iterations        =", result.iterations)
    print("attempted iterations        =", result.attempted_iterations)
    print("trial evaluations           =", result.trial_evaluations)
    print("final PGD rank              =", result.basis.n_modes)
    print("total modes added           =", result.total_modes_added)
    print("final LATIN indicator       = {:.12e}".format(result.final_indicator))
    print("best committed indicator    = {:.12e}".format(best))
    print("commit kinds                =", " ".join(result.commit_kind_history))
    print("modes added per commit      =", " ".join(str(int(x)) for x in modes_added))
    print()
    print("Committed LATIN indicator history:")
    print("  ", _history_text(indicators, ".9e"))
    print()
    print("Outer-iteration baseline indicator history:")
    print("  ", _history_text(baselines, ".9e"))
    print()
    print("Trial kinds:")
    print("  ", " ".join(result.trial_kind_history))
    print()
    print("Trial LATIN indicator history:")
    print("  ", _history_text(trials, ".9e"))
    print()
    print("Trial reduced-residual history:")
    print("  ", _history_text(reduced, ".9e"))

    print("-" * 116)
    print("Trend audit")
    print("  committed decreases       =", decreases)
    print("  committed increases       =", increases)
    print("  committed unchanged       =", unchanged)
    print("  final / first indicator   = {:.12e}".format(overall_ratio))

    enrichment = result.last_enrichment_result
    if enrichment is not None:
        print("-" * 116)
        print("Last enrichment diagnostic")
        print("  accepted                  =", enrichment.accepted)
        print("  failure reason            =", enrichment.failure_reason)
        print("  fixed-point converged     =", enrichment.fixed_point_converged)
        print("  fixed-point iterations    =", enrichment.fixed_point_iterations)
        if enrichment.fixed_point_history.size:
            print(
                "  final fixed-point chi     = {:.12e}".format(
                    float(enrichment.fixed_point_history[-1])
                )
            )
        print(
            "  residual benefit          = {:.6%}".format(
                enrichment.residual_benefit
            )
        )
        print(
            "  orthogonality error       = {:.12e}".format(
                enrichment.orthogonality_error
            )
        )

    print("-" * 116)
    print("Interpretation:")
    print(
        "  I-5 PASS requires the solver to move beyond the former attempted "
        "iteration 8 / rank-3 enrichment failure, preserve valid transactions, "
        "and show a stable overall decrease of the committed LATIN indicator. "
        "Absolute CONVERGED or the solver's accepted STAGNATED termination is "
        "stronger evidence; MAX_ITERATIONS is diagnostic and must be judged "
        "from the indicator trend rather than treated automatically as failure."
    )
    print("=" * 116)


if __name__ == "__main__":
    main()

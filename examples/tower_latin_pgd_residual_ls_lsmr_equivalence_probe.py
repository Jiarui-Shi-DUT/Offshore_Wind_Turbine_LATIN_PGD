# -*- coding: utf-8 -*-
"""I-3 diagnostic: dense residual-LS versus matrix-free LSMR residual-LS.

Diagnostic only
---------------
This script does NOT modify latin/ production code.

Background
----------
The successful fourth-mode Case C uses the residual-driven spatial problem

    min_p  1/2 sum_n w_n
           || lambda_dot_n p
              - lambda_n H_sigma,n A_sigma p
              + Delta_n ||^2_(M H_sigma,n^-1),

where

    A_sigma p = s

is the equilibrated tower stress mode returned by
TowerEquilibriumOperator.apply_spatial(p).

The earlier coarse diagnostic formed A_sigma explicitly as an Nq x Nq dense
matrix.  That is useful for reference calculations but is not scalable.

This probe constructs the ORIGINAL weighted least-squares operator directly as
a scipy.sparse.linalg.LinearOperator and solves it with LSMR.  It therefore
avoids both

1) a dense material-point stress matrix in the operator-driven solve, and
2) the normal equations B^T W B, whose condition number is approximately the
   square of the original least-squares operator condition number.

Weighted least-squares operator
-------------------------------
Define, for each time node n,

    B_n p = lambda_dot_n p
            - lambda_n H_sigma,n A_sigma p,

and

    D_n = diag( sqrt(w_n M / H_sigma,n) ).

Then the spatial problem is

    min_p || C p - b ||_2,

with

    C_n = D_n B_n,
    b_n = -D_n Delta_n.

The forward action requires only ONE equilibrium-operator call:

    A_sigma p.

For the adjoint, the tower stress map satisfies

    A_sigma^T M = M A_sigma,

hence

    A_sigma^T q = M A_sigma(M^-1 q).

Therefore C^T y can also be evaluated with only ONE equilibrium-operator call,
without assembling A_sigma or A_sigma^T.

Checks performed
----------------
A) LinearOperator audit against an explicit dense stacked least-squares matrix
   for one representative fourth-mode half-step:
   - forward action error,
   - adjoint action error,
   - Euclidean adjoint identity error,
   - dense stacked-LS versus LSMR solution error,
   - old normal-equation Case-C solution versus stacked-LS solution error.

B) Reproduce the old dense Case-C fourth-mode fixed point.

C) Run the same fourth-mode fixed point with the matrix-free LSMR spatial
   half-step and compare the converged mode with the dense Case-C reference.

The dense matrix is used ONLY in sections A/B as a diagnostic reference.
Section C is the candidate scalable operator path.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse.linalg import LinearOperator, lsmr

from examples.tower_latin_pgd_fourth_mode_direct_ls_spatial_probe import (
    FIXED_POINT_TOLERANCE,
    MAX_FIXED_POINT_ITERATIONS,
    MINIMUM_SPATIAL_NORM,
    RCOND,
    _build_problem,
    _dense_stress_matrix,
    _direct_1d_style_spatial_solve,
    _orthogonalize_and_normalise_like_1d,
    _rebuild_failing_trial_a,
    _trapezoidal_weights,
)
from latin.pgd_basis import PGDModeTower
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
import latin.tower_pgd_enrichment as enrichment_module


LSMR_ATOL = 1.0e-12
LSMR_BTOL = 1.0e-12
LSMR_CONLIM = 1.0e12
LSMR_MAX_ITERATIONS = 2000


@dataclass(frozen=True)
class LSMRResult:
    solution: np.ndarray
    converged: bool
    istop: int
    iterations: int
    residual_norm: float
    normal_residual_norm: float
    operator_norm_estimate: float
    condition_estimate: float


@dataclass(frozen=True)
class FixedPointResult:
    mode: PGDModeTower
    converged: bool
    history: np.ndarray
    spatial_solver_iterations: np.ndarray
    spatial_solver_condition_estimates: np.ndarray
    spatial_solver_all_converged: bool


def _metric_relative_error(
    left: np.ndarray,
    right: np.ndarray,
    metric: MaterialPointMetric,
) -> float:
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    numerator = metric.norm(left_array - right_array)
    denominator = max(
        metric.norm(right_array),
        np.finfo(np.float64).eps,
    )
    return float(numerator / denominator)


def _euclidean_relative_error(
    left: np.ndarray,
    right: np.ndarray,
) -> float:
    left_array = np.asarray(left, dtype=np.float64)
    right_array = np.asarray(right, dtype=np.float64)
    numerator = float(np.linalg.norm(left_array - right_array))
    denominator = max(
        float(np.linalg.norm(right_array)),
        np.finfo(np.float64).eps,
    )
    return float(numerator / denominator)


def _stress_action(
    operator: TowerEquilibriumOperator,
    vector: np.ndarray,
) -> np.ndarray:
    """Apply A_sigma without forming the dense material-point matrix."""
    return operator.apply_spatial(
        np.asarray(vector, dtype=np.float64)
    ).stress


def _build_matrix_free_ls_operator(
    *,
    temporal_amplitude: np.ndarray,
    temporal_rate: np.ndarray,
    defect: np.ndarray,
    time: np.ndarray,
    H_sigma: np.ndarray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
):
    """Return scipy LinearOperator C and RHS b for the weighted residual-LS."""
    lam = np.asarray(temporal_amplitude, dtype=np.float64)
    rate = np.asarray(temporal_rate, dtype=np.float64)
    d = np.asarray(defect, dtype=np.float64)
    t = np.asarray(time, dtype=np.float64)
    Hs = np.asarray(H_sigma, dtype=np.float64)

    nt = t.size
    nq = metric.n_material_points

    if lam.shape != (nt,):
        raise ValueError("temporal_amplitude has incompatible shape")
    if rate.shape != (nt,):
        raise ValueError("temporal_rate has incompatible shape")
    if d.shape != (nt, nq):
        raise ValueError("defect has incompatible shape")
    if Hs.shape != (nt, nq):
        raise ValueError("H_sigma has incompatible shape")
    if np.any(Hs <= 0.0):
        raise ValueError("H_sigma must be strictly positive")

    time_weights = _trapezoidal_weights(t)
    sqrt_weight = np.sqrt(
        time_weights[:, None]
        * metric.weights[None, :]
        / Hs
    )

    rhs = (-sqrt_weight * d).reshape(nt * nq)

    def matvec(vector):
        p = np.asarray(vector, dtype=np.float64).reshape(nq)
        Ap = _stress_action(operator, p)

        field = (
            rate[:, None] * p[None, :]
            - lam[:, None] * Hs * Ap[None, :]
        )
        result = sqrt_weight * field

        if np.any(~np.isfinite(result)):
            raise FloatingPointError(
                "matrix-free residual-LS forward action became non-finite"
            )
        return result.reshape(nt * nq)

    def rmatvec(vector):
        y = np.asarray(vector, dtype=np.float64).reshape(nt, nq)

        # z_n = D_n y_n
        z = sqrt_weight * y

        # Sum_n lambda_dot_n z_n
        direct_part = np.sum(
            rate[:, None] * z,
            axis=0,
        )

        # Sum_n lambda_n H_n z_n
        adjoint_source = np.sum(
            lam[:, None] * Hs * z,
            axis=0,
        )

        # A_sigma^T q = M A_sigma(M^-1 q)
        transformed_source = (
            adjoint_source / metric.weights
        )
        adjoint_stress = _stress_action(
            operator,
            transformed_source,
        )
        stress_adjoint_part = (
            metric.weights * adjoint_stress
        )

        result = direct_part - stress_adjoint_part

        if np.any(~np.isfinite(result)):
            raise FloatingPointError(
                "matrix-free residual-LS adjoint action became non-finite"
            )
        return np.asarray(result, dtype=np.float64)

    linear_operator = LinearOperator(
        shape=(nt * nq, nq),
        matvec=matvec,
        rmatvec=rmatvec,
        dtype=np.float64,
    )

    return linear_operator, rhs, sqrt_weight


def _solve_matrix_free_lsmr(
    *,
    temporal_amplitude: np.ndarray,
    temporal_rate: np.ndarray,
    defect: np.ndarray,
    time: np.ndarray,
    H_sigma: np.ndarray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
) -> LSMRResult:
    """Solve the original weighted LS problem with matrix-free LSMR."""
    linear_operator, rhs, _ = _build_matrix_free_ls_operator(
        temporal_amplitude=temporal_amplitude,
        temporal_rate=temporal_rate,
        defect=defect,
        time=time,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
    )

    result = lsmr(
        linear_operator,
        rhs,
        damp=0.0,
        atol=LSMR_ATOL,
        btol=LSMR_BTOL,
        conlim=LSMR_CONLIM,
        maxiter=LSMR_MAX_ITERATIONS,
        show=False,
    )

    solution = np.asarray(result[0], dtype=np.float64)
    istop = int(result[1])
    iterations = int(result[2])
    residual_norm = float(result[3])
    normal_residual_norm = float(result[4])
    operator_norm_estimate = float(result[5])
    condition_estimate = float(result[6])

    # scipy LSMR stop codes:
    # 1,2 : requested compatible / least-squares accuracy reached
    # 4,5 : same tests reached at machine-precision level
    converged = istop in (1, 2, 4, 5)

    if np.any(~np.isfinite(solution)):
        raise RuntimeError(
            "matrix-free LSMR returned non-finite spatial values"
        )
    if metric.norm(solution) <= MINIMUM_SPATIAL_NORM:
        raise RuntimeError(
            "matrix-free LSMR returned a degenerate spatial vector"
        )

    return LSMRResult(
        solution=solution,
        converged=converged,
        istop=istop,
        iterations=iterations,
        residual_norm=residual_norm,
        normal_residual_norm=normal_residual_norm,
        operator_norm_estimate=operator_norm_estimate,
        condition_estimate=condition_estimate,
    )


def _initial_mode(
    *,
    basis,
    time,
    defect,
    H_sigma,
    metric,
    operator,
    iteration_added,
) -> PGDModeTower:
    raw_seed = enrichment_module._seed(
        defect,
        H_sigma,
        metric,
        MINIMUM_SPATIAL_NORM,
    )
    p, _, _ = _orthogonalize_and_normalise_like_1d(
        raw_seed,
        basis,
        metric,
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
    return PGDModeTower(
        p,
        s,
        lam,
        rate,
        iteration_added,
    )


def _explicit_stacked_ls_matrix(
    *,
    temporal_amplitude: np.ndarray,
    temporal_rate: np.ndarray,
    time: np.ndarray,
    H_sigma: np.ndarray,
    metric: MaterialPointMetric,
    stress_matrix: np.ndarray,
):
    """Build dense C only for the coarse reference audit."""
    lam = np.asarray(temporal_amplitude, dtype=np.float64)
    rate = np.asarray(temporal_rate, dtype=np.float64)
    t = np.asarray(time, dtype=np.float64)
    Hs = np.asarray(H_sigma, dtype=np.float64)

    nt = t.size
    nq = metric.n_material_points

    time_weights = _trapezoidal_weights(t)
    sqrt_weight = np.sqrt(
        time_weights[:, None]
        * metric.weights[None, :]
        / Hs
    )

    identity = np.eye(nq, dtype=np.float64)
    blocks = []

    for n in range(nt):
        Bn = (
            rate[n] * identity
            - lam[n]
            * Hs[n, :, None]
            * stress_matrix
        )
        blocks.append(
            sqrt_weight[n, :, None] * Bn
        )

    return np.vstack(blocks), sqrt_weight


def _audit_linear_operator(
    *,
    current: PGDModeTower,
    defect: np.ndarray,
    time: np.ndarray,
    H_sigma: np.ndarray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
    stress_matrix: np.ndarray,
    basis,
):
    """Audit matrix-free forward/adjoint actions and one spatial solve."""
    linear_operator, rhs, _ = _build_matrix_free_ls_operator(
        temporal_amplitude=current.temporal_amplitude,
        temporal_rate=current.temporal_rate,
        defect=defect,
        time=time,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
    )

    dense_C, sqrt_weight = _explicit_stacked_ls_matrix(
        temporal_amplitude=current.temporal_amplitude,
        temporal_rate=current.temporal_rate,
        time=time,
        H_sigma=H_sigma,
        metric=metric,
        stress_matrix=stress_matrix,
    )
    dense_rhs = (-sqrt_weight * defect).reshape(-1)

    rhs_error = _euclidean_relative_error(
        rhs,
        dense_rhs,
    )

    rng = np.random.default_rng(20260820)
    probe_p = rng.standard_normal(
        metric.n_material_points
    )
    probe_y = rng.standard_normal(
        dense_C.shape[0]
    )

    forward_dense = dense_C @ probe_p
    forward_operator = linear_operator.matvec(probe_p)
    forward_error = _euclidean_relative_error(
        forward_operator,
        forward_dense,
    )

    adjoint_dense = dense_C.T @ probe_y
    adjoint_operator = linear_operator.rmatvec(probe_y)
    adjoint_error = _euclidean_relative_error(
        adjoint_operator,
        adjoint_dense,
    )

    lhs = float(
        np.dot(
            linear_operator.matvec(probe_p),
            probe_y,
        )
    )
    rhs_identity = float(
        np.dot(
            probe_p,
            linear_operator.rmatvec(probe_y),
        )
    )
    adjoint_identity_error = (
        abs(lhs - rhs_identity)
        / max(
            1.0,
            abs(lhs),
            abs(rhs_identity),
        )
    )

    dense_stacked_solution = np.linalg.lstsq(
        dense_C,
        dense_rhs,
        rcond=RCOND,
    )[0]

    old_normal_equation_solution = (
        _direct_1d_style_spatial_solve(
            temporal_amplitude=current.temporal_amplitude,
            temporal_rate=current.temporal_rate,
            defect=defect,
            time=time,
            H_sigma=H_sigma,
            metric=metric,
            stress_matrix=stress_matrix,
        )
    )

    lsmr_result = _solve_matrix_free_lsmr(
        temporal_amplitude=current.temporal_amplitude,
        temporal_rate=current.temporal_rate,
        defect=defect,
        time=time,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
    )

    lsmr_vs_stacked_raw = _metric_relative_error(
        lsmr_result.solution,
        dense_stacked_solution,
        metric,
    )
    normal_vs_stacked_raw = _metric_relative_error(
        old_normal_equation_solution,
        dense_stacked_solution,
        metric,
    )

    p_lsmr, _, _ = _orthogonalize_and_normalise_like_1d(
        lsmr_result.solution,
        basis,
        metric,
    )
    p_stacked, _, _ = _orthogonalize_and_normalise_like_1d(
        dense_stacked_solution,
        basis,
        metric,
    )
    p_old, _, _ = _orthogonalize_and_normalise_like_1d(
        old_normal_equation_solution,
        basis,
        metric,
    )

    lsmr_vs_stacked_normalised = _metric_relative_error(
        p_lsmr,
        p_stacked,
        metric,
    )
    old_vs_stacked_normalised = _metric_relative_error(
        p_old,
        p_stacked,
        metric,
    )

    print("LinearOperator audit:")
    print(
        "  RHS relative error                 = {:.12e}".format(
            rhs_error
        )
    )
    print(
        "  forward-action relative error      = {:.12e}".format(
            forward_error
        )
    )
    print(
        "  adjoint-action relative error      = {:.12e}".format(
            adjoint_error
        )
    )
    print(
        "  <Cp,y>-<p,C^Ty> relative error     = {:.12e}".format(
            adjoint_identity_error
        )
    )
    print(
        "  LSMR converged / istop             = {} / {}".format(
            lsmr_result.converged,
            lsmr_result.istop,
        )
    )
    print(
        "  LSMR iterations                    = {}".format(
            lsmr_result.iterations
        )
    )
    print(
        "  LSMR residual norm                 = {:.12e}".format(
            lsmr_result.residual_norm
        )
    )
    print(
        "  LSMR normal residual norm          = {:.12e}".format(
            lsmr_result.normal_residual_norm
        )
    )
    print(
        "  LSMR condition estimate            = {:.12e}".format(
            lsmr_result.condition_estimate
        )
    )
    print(
        "  LSMR vs dense stacked-LS raw error = {:.12e}".format(
            lsmr_vs_stacked_raw
        )
    )
    print(
        "  old normal-eq vs stacked raw error = {:.12e}".format(
            normal_vs_stacked_raw
        )
    )
    print(
        "  LSMR vs stacked norm-p error       = {:.12e}".format(
            lsmr_vs_stacked_normalised
        )
    )
    print(
        "  old vs stacked norm-p error        = {:.12e}".format(
            old_vs_stacked_normalised
        )
    )

    return {
        "rhs_error": rhs_error,
        "forward_error": forward_error,
        "adjoint_error": adjoint_error,
        "adjoint_identity_error": adjoint_identity_error,
        "lsmr_result": lsmr_result,
        "lsmr_vs_stacked_raw": lsmr_vs_stacked_raw,
        "normal_vs_stacked_raw": normal_vs_stacked_raw,
        "lsmr_vs_stacked_normalised": lsmr_vs_stacked_normalised,
        "old_vs_stacked_normalised": old_vs_stacked_normalised,
    }


def _run_dense_case_c(
    *,
    basis,
    time,
    defect,
    H_sigma,
    metric,
    operator,
    iteration_added,
    stress_matrix,
) -> FixedPointResult:
    """Reproduce the old dense Case-C fixed point and return the final mode."""
    current = _initial_mode(
        basis=basis,
        time=time,
        defect=defect,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
        iteration_added=iteration_added,
    )

    history = []

    for _ in range(MAX_FIXED_POINT_ITERATIONS):
        p_raw = _direct_1d_style_spatial_solve(
            temporal_amplitude=current.temporal_amplitude,
            temporal_rate=current.temporal_rate,
            defect=defect,
            time=time,
            H_sigma=H_sigma,
            metric=metric,
            stress_matrix=stress_matrix,
        )

        p, _, _ = _orthogonalize_and_normalise_like_1d(
            p_raw,
            basis,
            metric,
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

        if chi <= FIXED_POINT_TOLERANCE:
            break

    history_array = np.asarray(
        history,
        dtype=np.float64,
    )
    converged = bool(
        history_array.size
        and history_array[-1] <= FIXED_POINT_TOLERANCE
    )

    return FixedPointResult(
        mode=current,
        converged=converged,
        history=history_array,
        spatial_solver_iterations=np.zeros(
            history_array.size,
            dtype=np.int64,
        ),
        spatial_solver_condition_estimates=np.zeros(
            history_array.size,
            dtype=np.float64,
        ),
        spatial_solver_all_converged=True,
    )


def _run_operator_lsmr_case(
    *,
    basis,
    time,
    defect,
    H_sigma,
    metric,
    operator,
    iteration_added,
) -> FixedPointResult:
    """Run the fourth-mode fixed point using only matrix-free LSMR spatial solves."""
    current = _initial_mode(
        basis=basis,
        time=time,
        defect=defect,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
        iteration_added=iteration_added,
    )

    history = []
    solver_iterations = []
    condition_estimates = []
    all_solver_converged = True

    for fp_iteration in range(
        1,
        MAX_FIXED_POINT_ITERATIONS + 1,
    ):
        spatial_result = _solve_matrix_free_lsmr(
            temporal_amplitude=current.temporal_amplitude,
            temporal_rate=current.temporal_rate,
            defect=defect,
            time=time,
            H_sigma=H_sigma,
            metric=metric,
            operator=operator,
        )

        all_solver_converged = (
            all_solver_converged
            and spatial_result.converged
        )
        solver_iterations.append(
            spatial_result.iterations
        )
        condition_estimates.append(
            spatial_result.condition_estimate
        )

        if not spatial_result.converged:
            print(
                "  operator FP {:3d}: LSMR FAILED "
                "(istop={}, iters={}, cond~{:.3e})".format(
                    fp_iteration,
                    spatial_result.istop,
                    spatial_result.iterations,
                    spatial_result.condition_estimate,
                )
            )
            break

        p, _, _ = _orthogonalize_and_normalise_like_1d(
            spatial_result.solution,
            basis,
            metric,
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

        print(
            "  operator FP {:3d}: chi={:.9e}  "
            "LSMR iters={:4d}  istop={}  cond~{:.3e}".format(
                fp_iteration,
                chi,
                spatial_result.iterations,
                spatial_result.istop,
                spatial_result.condition_estimate,
            )
        )

        current = candidate

        if chi <= FIXED_POINT_TOLERANCE:
            break

    history_array = np.asarray(
        history,
        dtype=np.float64,
    )
    converged = bool(
        all_solver_converged
        and history_array.size
        and history_array[-1] <= FIXED_POINT_TOLERANCE
    )

    return FixedPointResult(
        mode=current,
        converged=converged,
        history=history_array,
        spatial_solver_iterations=np.asarray(
            solver_iterations,
            dtype=np.int64,
        ),
        spatial_solver_condition_estimates=np.asarray(
            condition_estimates,
            dtype=np.float64,
        ),
        spatial_solver_all_converged=all_solver_converged,
    )


def _equilibrium_relative(
    stress: np.ndarray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
) -> tuple[float, float]:
    residual = operator.equilibrium_residual(stress)
    absolute = float(np.linalg.norm(residual))

    scale_vector = (
        np.abs(operator.compatibility_matrix).T
        @ (
            metric.weights
            * np.abs(stress)
        )
    )
    scale = float(np.linalg.norm(scale_vector))
    relative = absolute / max(
        scale,
        np.finfo(np.float64).eps,
    )
    return absolute, float(relative)


def _raw_mode_residual_benefit(
    *,
    mode: PGDModeTower,
    defect: np.ndarray,
    time: np.ndarray,
    H_sigma: np.ndarray,
    metric: MaterialPointMetric,
):
    before = enrichment_module._be_residual_norm(
        defect,
        time,
        H_sigma,
        metric,
    )

    after_field = (
        defect
        + (
            mode.temporal_rate[:, None]
            * mode.spatial_plastic_strain[None, :]
        )
        - (
            H_sigma
            * mode.temporal_amplitude[:, None]
            * mode.spatial_stress[None, :]
        )
    )

    after = enrichment_module._be_residual_norm(
        after_field,
        time,
        H_sigma,
        metric,
    )

    benefit = (
        0.0
        if before <= np.finfo(np.float64).eps
        else 1.0 - after / before
    )
    return float(before), float(after), float(benefit)


def _tail_text(
    values: np.ndarray,
    n: int = 12,
) -> str:
    if values.size == 0:
        return "<empty>"
    tail = values[-min(n, values.size):]
    return " ".join(
        "{:.9f}".format(float(x))
        for x in tail
    )


def main() -> None:
    material, operator, initialization = _build_problem()
    baseline, directions, fixed_a = _rebuild_failing_trial_a(
        material,
        operator,
        initialization,
    )

    time = baseline.state.time
    defect = fixed_a.mechanical_residual
    H_sigma = directions.H_sigma
    metric = operator.metric
    basis = fixed_a.basis
    iteration_added = baseline.attempted_iterations

    print("=" * 120)
    print(
        "I-3 residual-LS: dense reference versus matrix-free ORIGINAL-LS LSMR operator path"
    )
    print("=" * 120)
    print(
        "baseline: termination={}, committed={}, attempted={}, "
        "rank={}, xi={:.9e}".format(
            baseline.termination_reason.value,
            baseline.iterations,
            baseline.attempted_iterations,
            baseline.basis.n_modes,
            baseline.final_indicator,
        )
    )
    print(
        "Trial-A relative residual = {:.9e}; "
        "Nq={}; Nt={}; fixed-point tol={:.1e}".format(
            fixed_a.relative_residual,
            operator.n_material_points,
            time.size,
            FIXED_POINT_TOLERANCE,
        )
    )
    print(
        "LSMR: atol={:.1e}; btol={:.1e}; conlim={:.1e}; maxiter={}".format(
            LSMR_ATOL,
            LSMR_BTOL,
            LSMR_CONLIM,
            LSMR_MAX_ITERATIONS,
        )
    )

    print("-" * 120)
    print(
        "Building dense A_sigma for diagnostic reference sections A/B only "
        "(the operator LSMR path does not use it)..."
    )
    stress_matrix = _dense_stress_matrix(operator)

    initial = _initial_mode(
        basis=basis,
        time=time,
        defect=defect,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
        iteration_added=iteration_added,
    )

    print("-" * 120)
    print(
        "A) Matrix-free LinearOperator + LSMR audit "
        "against explicit stacked dense least squares"
    )
    audit = _audit_linear_operator(
        current=initial,
        defect=defect,
        time=time,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
        stress_matrix=stress_matrix,
        basis=basis,
    )

    print("-" * 120)
    print("B) Dense normal-equation Case-C reference fixed point")
    dense_result = _run_dense_case_c(
        basis=basis,
        time=time,
        defect=defect,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
        iteration_added=iteration_added,
        stress_matrix=stress_matrix,
    )
    print(
        "  dense FP converged      = {}".format(
            dense_result.converged
        )
    )
    print(
        "  dense FP iterations     = {}".format(
            dense_result.history.size
        )
    )
    print(
        "  dense final chi         = {:.12e}".format(
            float(dense_result.history[-1])
            if dense_result.history.size
            else float("nan")
        )
    )
    print(
        "  dense chi tail          = {}".format(
            _tail_text(dense_result.history)
        )
    )

    print("-" * 120)
    print(
        "C) Fully matrix-free ORIGINAL-LS LSMR fourth-mode fixed point "
        "(no dense A_sigma in this solve)"
    )
    operator_result = _run_operator_lsmr_case(
        basis=basis,
        time=time,
        defect=defect,
        H_sigma=H_sigma,
        metric=metric,
        operator=operator,
        iteration_added=iteration_added,
    )

    dense_mode = dense_result.mode
    operator_mode = operator_result.mode

    pair_difference = enrichment_module._pair_change(
        dense_mode,
        operator_mode,
        time,
        H_sigma,
        metric,
    )
    spatial_error = _metric_relative_error(
        operator_mode.spatial_plastic_strain,
        dense_mode.spatial_plastic_strain,
        metric,
    )
    stress_error = _metric_relative_error(
        operator_mode.spatial_stress,
        dense_mode.spatial_stress,
        metric,
    )

    eq_abs, eq_rel = _equilibrium_relative(
        operator_mode.spatial_stress,
        metric,
        operator,
    )

    residual_before, residual_after, benefit = (
        _raw_mode_residual_benefit(
            mode=operator_mode,
            defect=defect,
            time=time,
            H_sigma=H_sigma,
            metric=metric,
        )
    )

    print()
    print("Operator-LSMR fixed-point summary:")
    print(
        "  operator FP converged        = {}".format(
            operator_result.converged
        )
    )
    print(
        "  all LSMR solves converged    = {}".format(
            operator_result.spatial_solver_all_converged
        )
    )
    print(
        "  operator FP iterations       = {}".format(
            operator_result.history.size
        )
    )
    print(
        "  operator final chi           = {:.12e}".format(
            float(operator_result.history[-1])
            if operator_result.history.size
            else float("nan")
        )
    )
    if operator_result.spatial_solver_iterations.size:
        print(
            "  max LSMR iterations          = {}".format(
                int(
                    np.max(
                        operator_result.spatial_solver_iterations
                    )
                )
            )
        )
    if operator_result.spatial_solver_condition_estimates.size:
        print(
            "  max LSMR condition estimate  = {:.12e}".format(
                float(
                    np.max(
                        operator_result.spatial_solver_condition_estimates
                    )
                )
            )
        )
    print(
        "  dense/operator pair diff     = {:.12e}".format(
            pair_difference
        )
    )
    print(
        "  final spatial relative error = {:.12e}".format(
            spatial_error
        )
    )
    print(
        "  final stress relative error  = {:.12e}".format(
            stress_error
        )
    )
    print(
        "  equilibrium ||H^T M s||      = {:.12e}".format(
            eq_abs
        )
    )
    print(
        "  equilibrium relative         = {:.12e}".format(
            eq_rel
        )
    )
    print(
        "  BE residual before           = {:.12e}".format(
            residual_before
        )
    )
    print(
        "  BE residual after            = {:.12e}".format(
            residual_after
        )
    )
    print(
        "  raw-mode residual benefit    = {:.6%}".format(
            benefit
        )
    )
    print(
        "  operator chi tail            = {}".format(
            _tail_text(operator_result.history)
        )
    )

    print("-" * 120)
    print("Interpretation gate:")
    print(
        "  Strong PASS: forward/adjoint audit errors ~ machine precision; "
        "LSMR agrees with explicit stacked LS; operator fourth-mode fixed "
        "point converges; equilibrium remains satisfied; residual benefit is "
        "positive; and the final operator mode stays close to the dense Case-C "
        "reference."
    )
    print(
        "  If LSMR and stacked LS agree but the old normal-equation Case-C "
        "differs materially, the discrepancy belongs to normal-equation "
        "conditioning/rank truncation rather than to the matrix-free operator "
        "derivation."
    )
    print(
        "  Diagnostic only: no latin/ production enrichment or transaction "
        "semantics are modified here."
    )
    print("=" * 120)


if __name__ == "__main__":
    main()

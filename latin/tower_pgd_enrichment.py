# -*- coding: utf-8 -*-
"""
One-mode PGD enrichment transaction for the tower LATIN-PGD global stage.

Sequence:
    shifted defect R_A
    -> residual-driven spatial seed
    -> optional in-loop M-orthogonalisation for residual-LS
    -> Eq. (72) temporal solve
    -> selectable spatial half-step:
         * paper_galerkin: tower Eq. (70)-(71)
         * residual_ls: matrix-free weighted residual-LS + LSMR
    -> M-norm scale normalisation
    -> optional in-loop M-orthogonalisation for residual-LS
    -> Eq. (72) temporal solve
    -> complete-pair fixed-point convergence
    -> M-weighted modified Gram-Schmidt
    -> exact temporal-coordinate transformation
    -> field-invariance / significance checks
    -> enlarged-basis Eq. (58)-(59) temporal re-optimisation with full forcing
    -> full residual benefit
    -> accept or reject.

The input FixedBasisPGDResult is a provisional Trial-A value and is never
mutated.  This module does not update hardening, damage, relaxation, xi, zeta,
or persistent solver state.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import MaterialPointMetric, TowerEquilibriumOperator
from latin.tower_pgd_time_update import FixedBasisPGDResult, update_tower_pgd_time_functions

FloatArray = NDArray[np.float64]

_SPATIAL_STRATEGY_PAPER_GALERKIN = "paper_galerkin"
_SPATIAL_STRATEGY_RESIDUAL_LS = "residual_ls"
_VALID_SPATIAL_STRATEGIES = (
    _SPATIAL_STRATEGY_PAPER_GALERKIN,
    _SPATIAL_STRATEGY_RESIDUAL_LS,
)


class _CandidateFailure(RuntimeError):
    def __init__(self, reason: str) -> None:
        super().__init__(reason)
        self.reason = str(reason)


def _ro(values: FloatArray, name: str, ndim: int) -> FloatArray:
    a = np.array(values, dtype=np.float64, copy=True)
    if a.ndim != ndim:
        raise ValueError(f"{name} must have {ndim} dimension(s).")
    if np.any(~np.isfinite(a)):
        raise ValueError(f"{name} contains non-finite values.")
    a.setflags(write=False)
    return a


def _nonnegative_integer(value: int, name: str) -> int:
    """Return a validated non-negative integer, rejecting Boolean aliases."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(f"{name} must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(f"{name} must be an integer.")
    result = int(value)
    if result < 0:
        raise ValueError(f"{name} must be non-negative.")
    return result


def _positive_integer(value: int, name: str) -> int:
    """Return a validated strictly positive integer."""
    result = _nonnegative_integer(value, name)
    if result < 1:
        raise ValueError(f"{name} must be at least one.")
    return result


@dataclass(frozen=True)
class TowerEnrichmentResult:
    accepted: bool
    failure_reason: Optional[str]
    candidate_fixed_basis_result: Optional[FixedBasisPGDResult]
    raw_mode: Optional[PGDModeTower]
    fixed_point_history: FloatArray
    fixed_point_iterations: int
    fixed_point_converged: bool
    projection_coefficients: FloatArray
    orthogonal_scale: float
    spatial_novelty: float
    temporal_significance: float
    orthogonality_error: float
    plastic_field_invariance_error: float
    plastic_rate_field_invariance_error: float
    stress_field_invariance_error: float
    residual_norm_before: float
    residual_norm_after: float
    residual_benefit: float

    def __post_init__(self) -> None:
        hist = _ro(self.fixed_point_history, "fixed_point_history", 1)
        coeff = _ro(self.projection_coefficients, "projection_coefficients", 1)
        if self.accepted:
            if self.failure_reason is not None or self.candidate_fixed_basis_result is None:
                raise ValueError("Accepted enrichment requires a candidate and no failure reason.")
        else:
            if self.candidate_fixed_basis_result is not None:
                raise ValueError("Rejected enrichment must not expose a candidate result.")
        object.__setattr__(self, "fixed_point_history", hist)
        object.__setattr__(self, "projection_coefficients", coeff)
        object.__setattr__(self, "accepted", bool(self.accepted))
        object.__setattr__(self, "fixed_point_converged", bool(self.fixed_point_converged))
        object.__setattr__(self, "fixed_point_iterations", int(self.fixed_point_iterations))
        for name in (
            "orthogonal_scale", "spatial_novelty", "temporal_significance",
            "orthogonality_error", "plastic_field_invariance_error",
            "plastic_rate_field_invariance_error", "stress_field_invariance_error",
            "residual_norm_before", "residual_norm_after", "residual_benefit",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite.")
            object.__setattr__(self, name, value)

    @property
    def n_modes(self) -> Optional[int]:
        if self.candidate_fixed_basis_result is None:
            return None
        return self.candidate_fixed_basis_result.n_modes


def _validate_inputs(
    fixed_basis_result: FixedBasisPGDResult,
    time: FloatArray,
    full_forcing: FloatArray,
    shifted_defect: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
) -> Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
    if not isinstance(fixed_basis_result, FixedBasisPGDResult):
        raise TypeError("fixed_basis_result must be a FixedBasisPGDResult.")
    if not isinstance(metric, MaterialPointMetric):
        raise TypeError("metric must be a MaterialPointMetric.")
    if not isinstance(equilibrium_operator, TowerEquilibriumOperator):
        raise TypeError(
            "equilibrium_operator must be a TowerEquilibriumOperator."
        )
    basis = fixed_basis_result.basis
    t = np.asarray(time, dtype=np.float64)
    f = np.asarray(full_forcing, dtype=np.float64)
    d = np.asarray(shifted_defect, dtype=np.float64)
    Hs = np.asarray(H_sigma, dtype=np.float64)

    if t.shape != (basis.n_time,) or np.any(~np.isfinite(t)) or np.any(np.diff(t) <= 0.0):
        raise ValueError("time must be finite, one-dimensional, and strictly increasing.")
    for name, a in (("full_forcing", f), ("shifted_defect", d), ("H_sigma", Hs)):
        if a.shape != basis.field_shape or np.any(~np.isfinite(a)):
            raise ValueError(f"{name} must have finite shape {basis.field_shape}.")
    if np.any(Hs <= 0.0):
        raise ValueError("H_sigma must be strictly positive.")
    if metric.n_material_points != basis.n_material_points:
        raise ValueError("metric and basis material-point counts differ.")
    if equilibrium_operator.n_material_points != basis.n_material_points:
        raise ValueError("operator and basis material-point counts differ.")
    if not np.array_equal(metric.weights, equilibrium_operator.metric.weights):
        raise ValueError("metric must match the equilibrium operator metric.")

    scale = max(1.0, float(np.linalg.norm(d)))
    if not np.allclose(
        fixed_basis_result.mechanical_residual,
        d,
        rtol=0.0,
        atol=1.0e-12 * scale,
    ):
        raise ValueError(
            "shifted_defect must equal the current Trial-A fixed-basis mechanical residual."
        )
    return t, f, d, Hs


def _be_residual_norm(field: FloatArray, time: FloatArray, H_sigma: FloatArray, metric: MaterialPointMetric) -> float:
    dt = np.diff(time)
    density = field[1:, :] ** 2 / H_sigma[1:, :]
    squared = float(np.dot(dt, density @ metric.weights))
    if squared < -1.0e-12:
        raise FloatingPointError("Negative squared BE residual norm.")
    return float(np.sqrt(max(0.0, squared)))


def _seed(defect: FloatArray, H_sigma: FloatArray, metric: MaterialPointMetric, minimum_spatial_norm: float) -> FloatArray:
    energy = np.sum(defect**2 * metric.weights[None, :] / H_sigma, axis=1)
    k = int(np.argmax(energy))
    if energy[k] <= np.finfo(np.float64).eps:
        raise _CandidateFailure("shifted_defect_negligible")
    p = -defect[k, :].copy()
    n = metric.norm(p)
    if n <= minimum_spatial_norm:
        raise _CandidateFailure("initial_spatial_seed_degenerate")
    return p / n


def _temporal_solve(
    p: FloatArray,
    s: FloatArray,
    time: FloatArray,
    defect: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
) -> Tuple[FloatArray, FloatArray]:
    """Sequential backward-Euler discretisation of Eq. (72)."""
    nt = time.size
    lam = np.zeros(nt, dtype=np.float64)
    rate = np.zeros(nt, dtype=np.float64)

    w0 = metric.weights / H_sigma[0, :]
    den0 = float(np.dot(w0 * p, p))
    if den0 <= np.finfo(np.float64).tiny:
        raise _CandidateFailure("temporal_initial_denominator_degenerate")
    rate[0] = -float(np.dot(w0 * p, defect[0, :])) / den0

    for n in range(1, nt):
        dt = float(time[n] - time[n - 1])
        Hn = H_sigma[n, :]
        w = metric.weights / Hn
        g = p / dt - Hn * s
        b = p * (lam[n - 1] / dt) - defect[n, :]
        den = float(np.dot(w * g, g))
        if den <= np.finfo(np.float64).tiny or not np.isfinite(den):
            raise _CandidateFailure("temporal_controllability_degenerate")
        lam[n] = float(np.dot(w * g, b)) / den
        rate[n] = (lam[n] - lam[n - 1]) / dt
        if not np.isfinite(lam[n]) or not np.isfinite(rate[n]):
            raise _CandidateFailure("temporal_solution_non_finite")
    return lam, rate


def _spatial_solve(
    lam: FloatArray,
    rate: FloatArray,
    time: FloatArray,
    defect: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
    minimum_spatial_norm: float,
) -> Tuple[FloatArray, FloatArray]:
    """Tower Eq. (70)-(71), followed only by M-norm scale normalisation."""
    dt = np.diff(time)
    lr = lam[1:]
    rr = rate[1:]
    a_h = float(np.dot(dt, lr * rr))
    if a_h <= np.finfo(np.float64).tiny or not np.isfinite(a_h):
        raise _CandidateFailure("spatial_temporal_contraction_degenerate")

    Aq = np.sum(dt[:, None] * H_sigma[1:, :] * lr[:, None] ** 2, axis=0)
    delta_bar = np.sum(dt[:, None] * defect[1:, :] * lr[:, None], axis=0)
    invW = Aq + a_h / operator.reference_modulus
    if np.any(~np.isfinite(invW)) or np.any(invW <= 0.0):
        raise _CandidateFailure("spatial_effective_operator_nonpositive")
    W = 1.0 / invW

    H = operator.compatibility_matrix
    mw = metric.weights * W
    K = H.T @ (mw[:, None] * H)
    K = 0.5 * (K + K.T)
    rhs = -(H.T @ (mw * delta_bar))
    try:
        L = np.linalg.cholesky(K)
        u = np.linalg.solve(L.T, np.linalg.solve(L, rhs))
    except np.linalg.LinAlgError as exc:
        raise _CandidateFailure("spatial_effective_stiffness_not_positive_definite") from exc

    eps_tilde = H @ u
    sigma_raw = W * (eps_tilde + delta_bar)
    p_raw = eps_tilde / a_h - sigma_raw / operator.reference_modulus
    if np.any(~np.isfinite(p_raw)):
        raise _CandidateFailure("spatial_solution_non_finite")
    norm = metric.norm(p_raw)
    if norm <= minimum_spatial_norm:
        raise _CandidateFailure("spatial_mode_degenerate")
    p = p_raw / norm
    s = operator.apply_spatial(p).stress
    return p, s


def _pair_norm(mode: PGDModeTower, time: FloatArray, H_sigma: FloatArray, metric: MaterialPointMetric) -> float:
    dt = np.diff(time)
    rate_field = mode.temporal_rate[1:, None] * mode.spatial_plastic_strain[None, :]
    stress_field = mode.temporal_amplitude[1:, None] * mode.spatial_stress[None, :]
    density = rate_field**2 / H_sigma[1:, :] + stress_field**2 * H_sigma[1:, :]
    squared = float(np.dot(dt, density @ metric.weights))
    return float(np.sqrt(max(0.0, squared)))


def _pair_change(previous: PGDModeTower, current: PGDModeTower, time: FloatArray, H_sigma: FloatArray, metric: MaterialPointMetric) -> float:
    dt = np.diff(time)
    dr = (
        current.plastic_strain_rate_correction()[1:, :]
        - previous.plastic_strain_rate_correction()[1:, :]
    )
    ds = current.stress_correction()[1:, :] - previous.stress_correction()[1:, :]
    density = dr**2 / H_sigma[1:, :] + ds**2 * H_sigma[1:, :]
    num = float(np.sqrt(max(0.0, float(np.dot(dt, density @ metric.weights)))))
    den = _pair_norm(previous, time, H_sigma, metric) + _pair_norm(current, time, H_sigma, metric)
    if den <= np.finfo(np.float64).eps:
        raise _CandidateFailure("fixed_point_pair_negligible")
    return num / den


def _inloop_orthogonalize_spatial(
    basis: PGDBasisTower,
    raw: FloatArray,
    metric: MaterialPointMetric,
    minimum_spatial_norm: float,
    rcond: float,
) -> FloatArray:
    """Project one raw fixed-point spatial iterate off the existing basis.

    This is the production-candidate counterpart of the in-loop
    M-orthogonalisation used by the successful fourth-mode Case C diagnostic.
    It acts only on the *new raw iterate*.  Existing basis coordinates are not
    changed here; the exact post-fixed-point coordinate transformation remains
    the responsibility of _post_fixed_point_transform().
    """
    work = np.asarray(raw, dtype=np.float64).copy()
    raw_norm = metric.norm(work)
    if not np.isfinite(raw_norm) or raw_norm <= minimum_spatial_norm:
        raise _CandidateFailure("inloop_spatial_iterate_degenerate")

    if basis.n_modes > 0:
        P = basis.spatial_plastic_strain_matrix()
        weighted_P = metric.weights[:, None] * P
        gram = P.T @ weighted_P
        rhs = P.T @ (metric.weights * work)
        coeff = np.linalg.lstsq(
            gram,
            rhs,
            rcond=rcond,
        )[0]
        work -= P @ coeff

    norm = metric.norm(work)
    if not np.isfinite(norm) or norm <= minimum_spatial_norm:
        raise _CandidateFailure(
            "inloop_spatial_candidate_linearly_dependent"
        )
    return work / norm


def _raw_fixed_point(
    basis: PGDBasisTower,
    time: FloatArray,
    defect: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
    iteration_added: int,
    fixed_point_tolerance: float,
    max_fixed_point_iterations: int,
    minimum_spatial_norm: float,
    spatial_strategy: str = _SPATIAL_STRATEGY_PAPER_GALERKIN,
    rcond: float = 1.0e-12,
) -> Tuple[PGDModeTower, FloatArray, bool]:
    p = _seed(defect, H_sigma, metric, minimum_spatial_norm)
    if spatial_strategy == _SPATIAL_STRATEGY_RESIDUAL_LS:
        p = _inloop_orthogonalize_spatial(
            basis,
            p,
            metric,
            minimum_spatial_norm,
            rcond,
        )
    s = operator.apply_spatial(p).stress
    lam, rate = _temporal_solve(p, s, time, defect, H_sigma, metric)
    current = PGDModeTower(p, s, lam, rate, iteration_added)
    history = []
    converged = False

    for _ in range(max_fixed_point_iterations):
        if spatial_strategy == _SPATIAL_STRATEGY_PAPER_GALERKIN:
            p, s = _spatial_solve(
                current.temporal_amplitude,
                current.temporal_rate,
                time,
                defect,
                H_sigma,
                metric,
                operator,
                minimum_spatial_norm,
            )
        elif spatial_strategy == _SPATIAL_STRATEGY_RESIDUAL_LS:
            # Lazy import keeps the existing paper-Galerkin path free from a
            # mandatory SciPy import unless the residual-LS candidate is used.
            from latin.tower_residual_ls_spatial import (
                solve_tower_residual_ls_spatial,
            )

            try:
                residual_ls = solve_tower_residual_ls_spatial(
                    temporal_amplitude=current.temporal_amplitude,
                    temporal_rate=current.temporal_rate,
                    defect=defect,
                    time=time,
                    H_sigma=H_sigma,
                    metric=metric,
                    equilibrium_operator=operator,
                    minimum_spatial_norm=minimum_spatial_norm,
                )
            except (ValueError, FloatingPointError, np.linalg.LinAlgError) as exc:
                raise _CandidateFailure(
                    "residual_ls_spatial_solver_failed"
                ) from exc

            if not residual_ls.converged:
                raise _CandidateFailure(
                    "residual_ls_spatial_solver_not_converged"
                )

            p = _inloop_orthogonalize_spatial(
                basis,
                residual_ls.spatial_plastic_strain,
                metric,
                minimum_spatial_norm,
                rcond,
            )
            s = operator.apply_spatial(p).stress
        else:
            # Public validation should make this unreachable.
            raise ValueError(
                f"Unknown spatial_strategy: {spatial_strategy!r}."
            )

        lam, rate = _temporal_solve(p, s, time, defect, H_sigma, metric)
        candidate = PGDModeTower(p, s, lam, rate, iteration_added)
        chi = _pair_change(current, candidate, time, H_sigma, metric)
        history.append(chi)
        current = candidate
        if chi <= fixed_point_tolerance:
            converged = True
            break

    return current, np.asarray(history, dtype=np.float64), converged


def _weighted_mgs(
    basis: PGDBasisTower,
    raw: FloatArray,
    metric: MaterialPointMetric,
    minimum_spatial_norm: float,
    minimum_spatial_novelty: float,
    passes: int,
) -> Tuple[FloatArray, FloatArray, float, float]:
    raw_norm = metric.norm(raw)
    if raw_norm <= minimum_spatial_norm:
        raise _CandidateFailure("raw_spatial_mode_negligible")
    work = np.asarray(raw, dtype=np.float64).copy()
    coeff = np.zeros(basis.n_modes, dtype=np.float64)

    # Accumulating coefficients preserves raw = P a + work exactly at every pass.
    for _ in range(passes):
        for j, mode in enumerate(basis.modes):
            pj = mode.spatial_plastic_strain
            denom = metric.inner_product(pj, pj)
            if denom <= np.finfo(np.float64).tiny:
                raise _CandidateFailure("existing_basis_mode_degenerate")
            alpha = metric.inner_product(pj, work) / denom
            work -= alpha * pj
            coeff[j] += alpha

    c = metric.norm(work)
    gamma_sp = c / raw_norm
    if c <= minimum_spatial_norm or gamma_sp <= minimum_spatial_novelty:
        raise _CandidateFailure("spatial_candidate_linearly_dependent")
    return work / c, coeff, float(c), float(gamma_sp)


def _temporal_norm(amplitude: FloatArray, time: FloatArray) -> float:
    return float(np.sqrt(max(0.0, float(np.dot(np.diff(time), amplitude[1:] ** 2)))))


def _relative_error(a: FloatArray, b: FloatArray) -> float:
    return float(np.linalg.norm(a - b) / max(1.0, float(np.linalg.norm(b))))


def _post_fixed_point_transform(
    basis: PGDBasisTower,
    raw_mode: PGDModeTower,
    time: FloatArray,
    metric: MaterialPointMetric,
    operator: TowerEquilibriumOperator,
    minimum_spatial_norm: float,
    minimum_spatial_novelty: float,
    mode_significance_tolerance: float,
    basis_health_tolerance: float,
    field_invariance_tolerance: float,
    reorthogonalization_passes: int,
) -> Tuple[PGDBasisTower, FloatArray, float, float, float, float, float, float, float]:
    p_new, coeff, c, gamma_sp = _weighted_mgs(
        basis,
        raw_mode.spatial_plastic_strain,
        metric,
        minimum_spatial_norm,
        minimum_spatial_novelty,
        reorthogonalization_passes,
    )
    s_new = operator.apply_spatial(p_new).stress

    Aold = basis.temporal_amplitude_matrix()
    Rold = basis.temporal_rate_matrix()
    Aplus = Aold + raw_mode.temporal_amplitude[:, None] * coeff[None, :]
    Rplus = Rold + raw_mode.temporal_rate[:, None] * coeff[None, :]
    transformed = basis.with_temporal_coordinates(Aplus, Rplus)
    transformed = transformed.with_appended(
        PGDModeTower(
            p_new,
            s_new,
            c * raw_mode.temporal_amplitude,
            c * raw_mode.temporal_rate,
            raw_mode.iteration_added,
        )
    )

    raw_plastic = basis.plastic_strain_correction() + raw_mode.plastic_strain_correction()
    raw_rate = basis.plastic_strain_rate_correction() + raw_mode.plastic_strain_rate_correction()
    raw_stress = basis.stress_correction() + raw_mode.stress_correction()
    ep_err = _relative_error(transformed.plastic_strain_correction(), raw_plastic)
    rate_err = _relative_error(transformed.plastic_strain_rate_correction(), raw_rate)
    stress_err = _relative_error(transformed.stress_correction(), raw_stress)
    if max(ep_err, rate_err, stress_err) > field_invariance_tolerance:
        raise _CandidateFailure("gram_schmidt_field_invariance_failed")

    P = transformed.spatial_plastic_strain_matrix()
    gram = P.T @ (metric.weights[:, None] * P)
    ortho_err = float(np.linalg.norm(gram - np.eye(transformed.n_modes), ord="fro"))
    if ortho_err > basis_health_tolerance:
        raise _CandidateFailure("basis_orthogonality_health_failed")

    norms = np.array([_temporal_norm(m.temporal_amplitude, time) for m in transformed.modes])
    total = float(np.linalg.norm(norms))
    gamma_lambda = 0.0 if total <= np.finfo(np.float64).eps else float(norms[-1] / total)
    if gamma_lambda <= mode_significance_tolerance:
        raise _CandidateFailure("modified_time_function_insignificant")

    return transformed, coeff, c, gamma_sp, gamma_lambda, ortho_err, ep_err, rate_err, stress_err


def _reject(
    reason: str,
    basis: PGDBasisTower,
    residual_before: float,
    raw_mode: Optional[PGDModeTower] = None,
    history: Optional[FloatArray] = None,
    fixed_converged: bool = False,
    coeff: Optional[FloatArray] = None,
    c: float = 0.0,
    gamma_sp: float = 0.0,
    gamma_lambda: float = 0.0,
    ortho_err: float = 0.0,
    ep_err: float = 0.0,
    rate_err: float = 0.0,
    stress_err: float = 0.0,
    residual_after: Optional[float] = None,
) -> TowerEnrichmentResult:
    h = np.zeros(0, dtype=np.float64) if history is None else history
    a = np.zeros(basis.n_modes, dtype=np.float64) if coeff is None else coeff
    ra = residual_before if residual_after is None else residual_after
    benefit = 0.0 if residual_before <= np.finfo(np.float64).eps else 1.0 - ra / residual_before
    return TowerEnrichmentResult(
        accepted=False,
        failure_reason=reason,
        candidate_fixed_basis_result=None,
        raw_mode=raw_mode,
        fixed_point_history=h,
        fixed_point_iterations=int(h.size),
        fixed_point_converged=fixed_converged,
        projection_coefficients=a,
        orthogonal_scale=c,
        spatial_novelty=gamma_sp,
        temporal_significance=gamma_lambda,
        orthogonality_error=ortho_err,
        plastic_field_invariance_error=ep_err,
        plastic_rate_field_invariance_error=rate_err,
        stress_field_invariance_error=stress_err,
        residual_norm_before=residual_before,
        residual_norm_after=ra,
        residual_benefit=benefit,
    )


def enrich_tower_pgd_basis_once(
    fixed_basis_result: FixedBasisPGDResult,
    time: FloatArray,
    full_forcing: FloatArray,
    shifted_defect: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
    *,
    mode_significance_tolerance: float,
    acceptance_tolerance: float,
    iteration_added: int = 0,
    fixed_point_tolerance: float = 1.0e-6,
    max_fixed_point_iterations: int = 30,
    minimum_spatial_norm: float = 1.0e-14,
    minimum_spatial_novelty: float = 1.0e-12,
    basis_health_tolerance: float = 1.0e-8,
    field_invariance_tolerance: float = 1.0e-10,
    reorthogonalization_passes: int = 2,
    reduced_tolerance: float = 1.0e-4,
    rcond: float = 1.0e-12,
    spatial_strategy: str = _SPATIAL_STRATEGY_PAPER_GALERKIN,
) -> TowerEnrichmentResult:
    """
    Attempt exactly one Add-a-pair transaction from provisional B_m^A.

    The two tower-specific usefulness thresholds are required keyword
    arguments because their calibrated values are intentionally not frozen by
    the theory-stage documents.
    """
    t, forcing, defect, Hs = _validate_inputs(
        fixed_basis_result,
        time,
        full_forcing,
        shifted_defect,
        H_sigma,
        metric,
        equilibrium_operator,
    )
    for name, value, allow_zero in (
        ("mode_significance_tolerance", mode_significance_tolerance, True),
        ("acceptance_tolerance", acceptance_tolerance, True),
        ("fixed_point_tolerance", fixed_point_tolerance, False),
        ("minimum_spatial_norm", minimum_spatial_norm, False),
        ("minimum_spatial_novelty", minimum_spatial_novelty, True),
        ("basis_health_tolerance", basis_health_tolerance, False),
        ("field_invariance_tolerance", field_invariance_tolerance, False),
        ("reduced_tolerance", reduced_tolerance, False),
        ("rcond", rcond, False),
    ):
        x = float(value)
        if not np.isfinite(x) or (x < 0.0 if allow_zero else x <= 0.0):
            raise ValueError(f"Invalid {name}.")
    iteration_added = _nonnegative_integer(
        iteration_added,
        "iteration_added",
    )
    max_fixed_point_iterations = _positive_integer(
        max_fixed_point_iterations,
        "max_fixed_point_iterations",
    )
    reorthogonalization_passes = _positive_integer(
        reorthogonalization_passes,
        "reorthogonalization_passes",
    )
    if reorthogonalization_passes not in (1, 2):
        raise ValueError("reorthogonalization_passes must be 1 or 2.")
    if not isinstance(spatial_strategy, str):
        raise TypeError("spatial_strategy must be a string.")
    if spatial_strategy not in _VALID_SPATIAL_STRATEGIES:
        raise ValueError(
            "spatial_strategy must be one of "
            f"{_VALID_SPATIAL_STRATEGIES}."
        )

    basis = fixed_basis_result.basis
    before = _be_residual_norm(defect, t, Hs, metric)
    if before <= np.finfo(np.float64).eps:
        return _reject("shifted_defect_negligible", basis, before)

    try:
        raw_mode, history, converged = _raw_fixed_point(
            basis,
            t,
            defect,
            Hs,
            metric,
            equilibrium_operator,
            iteration_added,
            float(fixed_point_tolerance),
            max_fixed_point_iterations,
            float(minimum_spatial_norm),
            spatial_strategy,
            float(rcond),
        )
    except _CandidateFailure as failure:
        return _reject(failure.reason, basis, before)

    if not converged:
        return _reject(
            "fixed_point_not_converged",
            basis,
            before,
            raw_mode=raw_mode,
            history=history,
            fixed_converged=False,
        )

    try:
        (
            tentative,
            coeff,
            c,
            gamma_sp,
            gamma_lambda,
            ortho_err,
            ep_err,
            rate_err,
            stress_err,
        ) = _post_fixed_point_transform(
            basis,
            raw_mode,
            t,
            metric,
            equilibrium_operator,
            float(minimum_spatial_norm),
            float(minimum_spatial_novelty),
            float(mode_significance_tolerance),
            float(basis_health_tolerance),
            float(field_invariance_tolerance),
            reorthogonalization_passes,
        )
    except _CandidateFailure as failure:
        return _reject(
            failure.reason,
            basis,
            before,
            raw_mode=raw_mode,
            history=history,
            fixed_converged=True,
        )

    try:
        candidate = update_tower_pgd_time_functions(
            basis=tentative,
            time=t,
            forcing=forcing,
            H_sigma=Hs,
            metric=metric,
            equilibrium_operator=equilibrium_operator,
            reduced_tolerance=float(reduced_tolerance),
            rcond=float(rcond),
        )
    except (ValueError, FloatingPointError, np.linalg.LinAlgError):
        return _reject(
            "all_mode_temporal_reoptimization_failed",
            basis,
            before,
            raw_mode=raw_mode,
            history=history,
            fixed_converged=True,
            coeff=coeff,
            c=c,
            gamma_sp=gamma_sp,
            gamma_lambda=gamma_lambda,
            ortho_err=ortho_err,
            ep_err=ep_err,
            rate_err=rate_err,
            stress_err=stress_err,
        )

    after = _be_residual_norm(candidate.mechanical_residual, t, Hs, metric)
    benefit = 1.0 - after / before
    if benefit <= float(acceptance_tolerance):
        return _reject(
            "full_residual_benefit_insufficient",
            basis,
            before,
            raw_mode=raw_mode,
            history=history,
            fixed_converged=True,
            coeff=coeff,
            c=c,
            gamma_sp=gamma_sp,
            gamma_lambda=gamma_lambda,
            ortho_err=ortho_err,
            ep_err=ep_err,
            rate_err=rate_err,
            stress_err=stress_err,
            residual_after=after,
        )

    return TowerEnrichmentResult(
        accepted=True,
        failure_reason=None,
        candidate_fixed_basis_result=candidate,
        raw_mode=raw_mode,
        fixed_point_history=history,
        fixed_point_iterations=int(history.size),
        fixed_point_converged=True,
        projection_coefficients=coeff,
        orthogonal_scale=c,
        spatial_novelty=gamma_sp,
        temporal_significance=gamma_lambda,
        orthogonality_error=ortho_err,
        plastic_field_invariance_error=ep_err,
        plastic_rate_field_invariance_error=rate_err,
        stress_field_invariance_error=stress_err,
        residual_norm_before=before,
        residual_norm_after=after,
        residual_benefit=benefit,
    )

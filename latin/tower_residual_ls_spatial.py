# -*- coding: utf-8 -*-
"""Matrix-free residual-LS spatial solver for tower LATIN-PGD.

This module is a production-candidate component.  It does not alter the
existing enrichment transaction or select this solver automatically.

For fixed temporal functions lambda(t) and lambda_dot(t), solve

    min_p  1/2 sum_n w_n
           || lambda_dot_n p
              - lambda_n H_sigma,n A_sigma p
              + defect_n ||^2_(M H_sigma,n^-1),

where A_sigma p is the equilibrated stress mode returned by
TowerEquilibriumOperator.apply_spatial(p).

The least-squares operator is applied matrix-free through scipy LSMR; no dense
Nq x Nq material-point stress matrix and no normal equations are formed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import NDArray
from scipy.sparse.linalg import LinearOperator, lsmr

from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)

FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class TowerResidualLSSpatialResult:
    """Result of one matrix-free residual-LS spatial half-step."""

    spatial_plastic_strain: FloatArray
    spatial_stress: FloatArray
    raw_spatial_norm: float
    converged: bool
    istop: int
    iterations: int
    residual_norm: float
    normal_residual_norm: float
    operator_norm_estimate: float
    condition_estimate: float

    def __post_init__(self) -> None:
        p = np.array(
            self.spatial_plastic_strain,
            dtype=np.float64,
            copy=True,
        )
        s = np.array(
            self.spatial_stress,
            dtype=np.float64,
            copy=True,
        )
        if p.ndim != 1 or s.shape != p.shape:
            raise ValueError(
                "spatial_plastic_strain and spatial_stress must be "
                "one-dimensional arrays with identical shape."
            )
        if np.any(~np.isfinite(p)) or np.any(~np.isfinite(s)):
            raise ValueError("Spatial result contains non-finite values.")
        p.setflags(write=False)
        s.setflags(write=False)
        object.__setattr__(self, "spatial_plastic_strain", p)
        object.__setattr__(self, "spatial_stress", s)
        object.__setattr__(self, "converged", bool(self.converged))
        object.__setattr__(self, "istop", int(self.istop))
        object.__setattr__(self, "iterations", int(self.iterations))
        for name in (
            "raw_spatial_norm",
            "residual_norm",
            "normal_residual_norm",
            "operator_norm_estimate",
            "condition_estimate",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(f"{name} must be finite.")
            object.__setattr__(self, name, value)


def _trapezoidal_weights(time: FloatArray) -> FloatArray:
    t = np.asarray(time, dtype=np.float64)
    if (
        t.ndim != 1
        or t.size < 2
        or np.any(~np.isfinite(t))
        or np.any(np.diff(t) <= 0.0)
    ):
        raise ValueError(
            "time must be finite, one-dimensional, strictly increasing, "
            "and contain at least two points."
        )

    w = np.empty_like(t)
    w[0] = 0.5 * (t[1] - t[0])
    w[-1] = 0.5 * (t[-1] - t[-2])
    if t.size > 2:
        w[1:-1] = 0.5 * (t[2:] - t[:-2])
    return w


def _validate_inputs(
    temporal_amplitude: FloatArray,
    temporal_rate: FloatArray,
    defect: FloatArray,
    time: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
) -> Tuple[
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
]:
    if not isinstance(metric, MaterialPointMetric):
        raise TypeError("metric must be a MaterialPointMetric.")
    if not isinstance(
        equilibrium_operator,
        TowerEquilibriumOperator,
    ):
        raise TypeError(
            "equilibrium_operator must be a TowerEquilibriumOperator."
        )
    if (
        metric.n_material_points
        != equilibrium_operator.n_material_points
    ):
        raise ValueError(
            "metric and equilibrium operator material-point counts differ."
        )
    if not np.array_equal(
        metric.weights,
        equilibrium_operator.metric.weights,
    ):
        raise ValueError(
            "metric must match the equilibrium operator metric."
        )

    t = np.asarray(time, dtype=np.float64)
    _trapezoidal_weights(t)

    nt = t.size
    nq = metric.n_material_points

    lam = np.asarray(temporal_amplitude, dtype=np.float64)
    rate = np.asarray(temporal_rate, dtype=np.float64)
    d = np.asarray(defect, dtype=np.float64)
    Hs = np.asarray(H_sigma, dtype=np.float64)

    if lam.shape != (nt,) or np.any(~np.isfinite(lam)):
        raise ValueError(
            f"temporal_amplitude must have finite shape {(nt,)}."
        )
    if rate.shape != (nt,) or np.any(~np.isfinite(rate)):
        raise ValueError(
            f"temporal_rate must have finite shape {(nt,)}."
        )
    if d.shape != (nt, nq) or np.any(~np.isfinite(d)):
        raise ValueError(
            f"defect must have finite shape {(nt, nq)}."
        )
    if Hs.shape != (nt, nq) or np.any(~np.isfinite(Hs)):
        raise ValueError(
            f"H_sigma must have finite shape {(nt, nq)}."
        )
    if np.any(Hs <= 0.0):
        raise ValueError("H_sigma must be strictly positive.")

    return lam, rate, d, t, Hs


def build_tower_residual_ls_linear_operator(
    temporal_amplitude: FloatArray,
    temporal_rate: FloatArray,
    defect: FloatArray,
    time: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
) -> Tuple[LinearOperator, FloatArray]:
    """Build the matrix-free weighted least-squares operator C and RHS b."""
    lam, rate, d, t, Hs = _validate_inputs(
        temporal_amplitude,
        temporal_rate,
        defect,
        time,
        H_sigma,
        metric,
        equilibrium_operator,
    )

    nt = t.size
    nq = metric.n_material_points
    time_weights = _trapezoidal_weights(t)

    sqrt_weight = np.sqrt(
        time_weights[:, None]
        * metric.weights[None, :]
        / Hs
    )
    rhs = (-sqrt_weight * d).reshape(nt * nq)

    def stress_action(vector: FloatArray) -> FloatArray:
        return equilibrium_operator.apply_spatial(
            np.asarray(vector, dtype=np.float64)
        ).stress

    def matvec(vector: FloatArray) -> FloatArray:
        p = np.asarray(vector, dtype=np.float64).reshape(nq)
        Ap = stress_action(p)

        field = (
            rate[:, None] * p[None, :]
            - lam[:, None] * Hs * Ap[None, :]
        )
        result = sqrt_weight * field
        if np.any(~np.isfinite(result)):
            raise FloatingPointError(
                "Residual-LS forward action became non-finite."
            )
        return result.reshape(nt * nq)

    def rmatvec(vector: FloatArray) -> FloatArray:
        y = np.asarray(vector, dtype=np.float64).reshape(nt, nq)

        # z_n = D_n y_n.
        z = sqrt_weight * y

        direct_part = np.sum(
            rate[:, None] * z,
            axis=0,
        )

        # q = sum_n lambda_n H_sigma,n z_n.
        adjoint_source = np.sum(
            lam[:, None] * Hs * z,
            axis=0,
        )

        # A_sigma^T q = M A_sigma(M^-1 q), because
        # A_sigma^T M = M A_sigma.
        transformed_source = (
            adjoint_source / metric.weights
        )
        adjoint_stress = stress_action(
            transformed_source
        )
        stress_adjoint_part = (
            metric.weights * adjoint_stress
        )

        result = direct_part - stress_adjoint_part
        if np.any(~np.isfinite(result)):
            raise FloatingPointError(
                "Residual-LS adjoint action became non-finite."
            )
        return np.asarray(result, dtype=np.float64)

    operator = LinearOperator(
        shape=(nt * nq, nq),
        matvec=matvec,
        rmatvec=rmatvec,
        dtype=np.float64,
    )
    return operator, np.asarray(rhs, dtype=np.float64)


def solve_tower_residual_ls_spatial(
    temporal_amplitude: FloatArray,
    temporal_rate: FloatArray,
    defect: FloatArray,
    time: FloatArray,
    H_sigma: FloatArray,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
    *,
    minimum_spatial_norm: float = 1.0e-14,
    atol: float = 1.0e-12,
    btol: float = 1.0e-12,
    conlim: float = 1.0e12,
    max_iterations: int = 2000,
) -> TowerResidualLSSpatialResult:
    """Solve one matrix-free residual-LS spatial half-step with LSMR."""
    for name, value in (
        ("minimum_spatial_norm", minimum_spatial_norm),
        ("atol", atol),
        ("btol", btol),
        ("conlim", conlim),
    ):
        x = float(value)
        if not np.isfinite(x) or x <= 0.0:
            raise ValueError(f"{name} must be finite and strictly positive.")

    if isinstance(max_iterations, (bool, np.bool_)):
        raise TypeError("max_iterations must be an integer.")
    if not isinstance(max_iterations, (int, np.integer)):
        raise TypeError("max_iterations must be an integer.")
    maxiter = int(max_iterations)
    if maxiter < 1:
        raise ValueError("max_iterations must be at least one.")

    linear_operator, rhs = (
        build_tower_residual_ls_linear_operator(
            temporal_amplitude,
            temporal_rate,
            defect,
            time,
            H_sigma,
            metric,
            equilibrium_operator,
        )
    )

    output = lsmr(
        linear_operator,
        rhs,
        damp=0.0,
        atol=float(atol),
        btol=float(btol),
        conlim=float(conlim),
        maxiter=maxiter,
        show=False,
    )

    raw = np.asarray(output[0], dtype=np.float64)
    istop = int(output[1])
    iterations = int(output[2])
    residual_norm = float(output[3])
    normal_residual_norm = float(output[4])
    operator_norm_estimate = float(output[5])
    condition_estimate = float(output[6])

    if np.any(~np.isfinite(raw)):
        raise FloatingPointError(
            "Residual-LS LSMR solution became non-finite."
        )

    raw_norm = metric.norm(raw)
    if (
        not np.isfinite(raw_norm)
        or raw_norm <= float(minimum_spatial_norm)
    ):
        raise FloatingPointError(
            "Residual-LS spatial mode is degenerate."
        )

    p = raw / raw_norm
    s = equilibrium_operator.apply_spatial(p).stress

    # scipy LSMR stop codes:
    # 1,2 = requested compatible / least-squares accuracy reached;
    # 4,5 = equivalent tests reached at machine precision.
    converged = istop in (1, 2, 4, 5)

    return TowerResidualLSSpatialResult(
        spatial_plastic_strain=p,
        spatial_stress=s,
        raw_spatial_norm=float(raw_norm),
        converged=converged,
        istop=istop,
        iterations=iterations,
        residual_norm=residual_norm,
        normal_residual_norm=normal_residual_norm,
        operator_norm_estimate=operator_norm_estimate,
        condition_estimate=condition_estimate,
    )

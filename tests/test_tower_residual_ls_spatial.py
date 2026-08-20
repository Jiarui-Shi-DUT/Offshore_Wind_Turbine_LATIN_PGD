# -*- coding: utf-8 -*-
"""Unit tests for the matrix-free tower residual-LS spatial solver."""

import unittest

import numpy as np

from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_residual_ls_spatial import (
    _trapezoidal_weights,
    build_tower_residual_ls_linear_operator,
    solve_tower_residual_ls_spatial,
)
from latin.tower_state import MaterialPointLayout


def make_operator():
    layout = MaterialPointLayout(1, 2, 2)
    metric = MaterialPointMetric(
        np.array([1.0, 2.0, 1.5, 0.5])
    )
    H = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, -1.0],
        ],
        dtype=np.float64,
    )
    operator = TowerEquilibriumOperator(
        layout=layout,
        metric=metric,
        reference_modulus=np.array(
            [10.0, 20.0, 15.0, 25.0]
        ),
        compatibility_matrix=H,
        free_dofs=np.array([1, 3]),
        n_dof=4,
    )
    return metric, operator


def dense_stress_matrix(operator):
    nq = operator.n_material_points
    eye = np.eye(nq, dtype=np.float64)
    return np.column_stack(
        [
            operator.apply_spatial(
                eye[:, j]
            ).stress
            for j in range(nq)
        ]
    )


def explicit_stacked_ls(
    temporal_amplitude,
    temporal_rate,
    defect,
    time,
    H_sigma,
    metric,
    stress_matrix,
):
    lam = np.asarray(
        temporal_amplitude,
        dtype=np.float64,
    )
    rate = np.asarray(
        temporal_rate,
        dtype=np.float64,
    )
    t = np.asarray(time, dtype=np.float64)
    Hs = np.asarray(
        H_sigma,
        dtype=np.float64,
    )
    d = np.asarray(defect, dtype=np.float64)

    nt = t.size
    nq = metric.n_material_points
    weights = _trapezoidal_weights(t)
    sqrt_weight = np.sqrt(
        weights[:, None]
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

    C = np.vstack(blocks)
    rhs = (-sqrt_weight * d).reshape(-1)
    return C, rhs


class TestTowerResidualLSSpatial(unittest.TestCase):
    def setUp(self):
        self.metric, self.operator = make_operator()
        self.time = np.array(
            [0.0, 0.2, 0.55, 0.95, 1.4],
            dtype=np.float64,
        )
        self.amplitude = np.array(
            [0.0, 0.08, 0.21, 0.15, 0.29],
            dtype=np.float64,
        )
        self.rate = np.array(
            [0.17, 0.40, 0.3714285714285714, -0.15, 0.3111111111111111],
            dtype=np.float64,
        )
        self.H_sigma = (
            0.18
            + 0.01
            * np.arange(self.time.size)[:, None]
            + 0.006
            * np.arange(
                self.metric.n_material_points
            )[None, :]
        )

        rng = np.random.default_rng(20260820)
        self.defect = 0.03 * rng.standard_normal(
            (
                self.time.size,
                self.metric.n_material_points,
            )
        )

    def test_linear_operator_matches_explicit_dense_actions(self):
        linear_operator, rhs = (
            build_tower_residual_ls_linear_operator(
                self.amplitude,
                self.rate,
                self.defect,
                self.time,
                self.H_sigma,
                self.metric,
                self.operator,
            )
        )

        A_sigma = dense_stress_matrix(
            self.operator
        )
        C, dense_rhs = explicit_stacked_ls(
            self.amplitude,
            self.rate,
            self.defect,
            self.time,
            self.H_sigma,
            self.metric,
            A_sigma,
        )

        np.testing.assert_allclose(
            rhs,
            dense_rhs,
            rtol=0.0,
            atol=0.0,
        )

        rng = np.random.default_rng(7)
        p = rng.standard_normal(
            self.metric.n_material_points
        )
        y = rng.standard_normal(C.shape[0])

        np.testing.assert_allclose(
            linear_operator.matvec(p),
            C @ p,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            linear_operator.rmatvec(y),
            C.T @ y,
            rtol=1.0e-12,
            atol=1.0e-12,
        )

        lhs = float(
            np.dot(
                linear_operator.matvec(p),
                y,
            )
        )
        rhs_identity = float(
            np.dot(
                p,
                linear_operator.rmatvec(y),
            )
        )
        self.assertLess(
            abs(lhs - rhs_identity)
            / max(
                1.0,
                abs(lhs),
                abs(rhs_identity),
            ),
            1.0e-12,
        )

    def test_lsmr_solution_matches_explicit_dense_least_squares(self):
        A_sigma = dense_stress_matrix(
            self.operator
        )
        C, rhs = explicit_stacked_ls(
            self.amplitude,
            self.rate,
            self.defect,
            self.time,
            self.H_sigma,
            self.metric,
            A_sigma,
        )
        dense_raw = np.linalg.lstsq(
            C,
            rhs,
            rcond=1.0e-12,
        )[0]
        dense_p = (
            dense_raw
            / self.metric.norm(dense_raw)
        )
        dense_s = self.operator.apply_spatial(
            dense_p
        ).stress

        result = solve_tower_residual_ls_spatial(
            self.amplitude,
            self.rate,
            self.defect,
            self.time,
            self.H_sigma,
            self.metric,
            self.operator,
            atol=1.0e-12,
            btol=1.0e-12,
            conlim=1.0e12,
            max_iterations=2000,
        )

        self.assertTrue(
            result.converged,
            msg=(
                f"LSMR istop={result.istop}, "
                f"iterations={result.iterations}"
            ),
        )
        np.testing.assert_allclose(
            result.spatial_plastic_strain,
            dense_p,
            rtol=1.0e-9,
            atol=1.0e-10,
        )
        np.testing.assert_allclose(
            result.spatial_stress,
            dense_s,
            rtol=1.0e-9,
            atol=1.0e-10,
        )
        self.assertAlmostEqual(
            self.metric.norm(
                result.spatial_plastic_strain
            ),
            1.0,
            places=12,
        )

    def test_returned_stress_mode_is_equilibrated(self):
        result = solve_tower_residual_ls_spatial(
            self.amplitude,
            self.rate,
            self.defect,
            self.time,
            self.H_sigma,
            self.metric,
            self.operator,
        )
        residual = self.operator.equilibrium_residual(
            result.spatial_stress
        )
        self.assertLess(
            float(np.linalg.norm(residual)),
            1.0e-10,
        )

    def test_invalid_H_sigma_is_rejected(self):
        bad = self.H_sigma.copy()
        bad[2, 1] = 0.0
        with self.assertRaises(ValueError):
            solve_tower_residual_ls_spatial(
                self.amplitude,
                self.rate,
                self.defect,
                self.time,
                bad,
                self.metric,
                self.operator,
            )


if __name__ == "__main__":
    unittest.main()

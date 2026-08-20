# -*- coding: utf-8 -*-
"""Integration tests for residual-LS inside the tower enrichment transaction."""

import unittest

import numpy as np

from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_pgd_enrichment import enrich_tower_pgd_basis_once
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
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


def be_rate(amplitude, time):
    rate = np.zeros_like(amplitude)
    rate[0] = (
        (amplitude[1] - amplitude[0])
        / (time[1] - time[0])
    )
    rate[1:] = (
        np.diff(amplitude)
        / np.diff(time)
    )
    return rate


def two_mode_problem():
    metric, operator = make_operator()
    time = np.array(
        [0.0, 0.25, 0.55, 0.95, 1.35],
        dtype=np.float64,
    )

    p1 = np.array(
        [0.5, -0.2, 0.8, 0.1],
        dtype=np.float64,
    )
    p1 /= metric.norm(p1)

    q = np.array(
        [-0.3, 0.7, 0.1, 0.6],
        dtype=np.float64,
    )
    q -= metric.inner_product(p1, q) * p1
    p2 = q / metric.norm(q)

    s1 = operator.apply_spatial(p1).stress
    s2 = operator.apply_spatial(p2).stress

    a1 = np.array(
        [0.0, 0.08, 0.18, 0.12, 0.25],
        dtype=np.float64,
    )
    a2 = np.array(
        [0.0, -0.04, 0.06, 0.13, 0.05],
        dtype=np.float64,
    )
    r1 = be_rate(a1, time)
    r2 = be_rate(a2, time)

    H_sigma = (
        0.18
        + 0.008 * np.arange(time.size)[:, None]
        + 0.004 * np.arange(4)[None, :]
    )

    forcing = (
        np.outer(r1, p1)
        - H_sigma * np.outer(a1, s1)
        + np.outer(r2, p2)
        - H_sigma * np.outer(a2, s2)
    )

    basis = PGDBasisTower(
        4,
        time.size,
        (
            PGDModeTower(
                p1,
                s1,
                np.zeros(time.size),
                np.zeros(time.size),
            ),
        ),
    )
    trial_a = update_tower_pgd_time_functions(
        basis,
        time,
        forcing,
        H_sigma,
        metric,
        operator,
    )
    return (
        metric,
        operator,
        time,
        H_sigma,
        forcing,
        trial_a,
    )


class TestTowerPGDResidualLSEnrichment(unittest.TestCase):
    def test_default_strategy_matches_explicit_paper_galerkin(self):
        metric, operator = make_operator()
        time = np.array(
            [0.0, 0.2, 0.5, 0.9, 1.4],
            dtype=np.float64,
        )

        p = np.array(
            [0.6, -0.4, 0.9, -0.2],
            dtype=np.float64,
        )
        p /= metric.norm(p)
        s = operator.apply_spatial(p).stress
        amplitude = np.array(
            [0.0, 0.08, 0.20, 0.14, 0.30],
            dtype=np.float64,
        )
        rate = be_rate(amplitude, time)
        H_sigma = (
            0.20
            + 0.01 * np.arange(time.size)[:, None]
            + 0.005 * np.arange(4)[None, :]
        )
        forcing = (
            np.outer(rate, p)
            - H_sigma * np.outer(amplitude, s)
        )

        empty = PGDBasisTower(4, time.size)
        trial_a = update_tower_pgd_time_functions(
            empty,
            time,
            forcing,
            H_sigma,
            metric,
            operator,
        )

        common = dict(
            fixed_basis_result=trial_a,
            time=time,
            full_forcing=forcing,
            shifted_defect=trial_a.mechanical_residual,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-8,
            max_fixed_point_iterations=100,
        )

        default = enrich_tower_pgd_basis_once(
            **common
        )
        explicit = enrich_tower_pgd_basis_once(
            **common,
            spatial_strategy="paper_galerkin",
        )

        self.assertEqual(
            default.accepted,
            explicit.accepted,
        )
        self.assertEqual(
            default.failure_reason,
            explicit.failure_reason,
        )
        self.assertEqual(
            default.fixed_point_iterations,
            explicit.fixed_point_iterations,
        )
        np.testing.assert_allclose(
            default.fixed_point_history,
            explicit.fixed_point_history,
            rtol=0.0,
            atol=0.0,
        )

    def test_residual_ls_accepts_rank_growth(self):
        (
            metric,
            operator,
            time,
            H_sigma,
            forcing,
            trial_a,
        ) = two_mode_problem()

        before_rank = trial_a.basis.n_modes
        before_residual = np.linalg.norm(
            trial_a.mechanical_residual
        )

        result = enrich_tower_pgd_basis_once(
            trial_a,
            time,
            forcing,
            trial_a.mechanical_residual,
            H_sigma,
            metric,
            operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-8,
            max_fixed_point_iterations=200,
            spatial_strategy="residual_ls",
        )

        self.assertTrue(
            result.accepted,
            msg=result.failure_reason,
        )
        self.assertTrue(
            result.fixed_point_converged
        )
        self.assertEqual(
            result.n_modes,
            before_rank + 1,
        )
        self.assertGreater(
            result.residual_benefit,
            0.0,
        )
        self.assertLess(
            result.orthogonality_error,
            1.0e-8,
        )
        self.assertLess(
            np.linalg.norm(
                result
                .candidate_fixed_basis_result
                .mechanical_residual
            ),
            before_residual,
        )

    def test_residual_ls_rejection_preserves_trial_a(self):
        (
            metric,
            operator,
            time,
            H_sigma,
            forcing,
            trial_a,
        ) = two_mode_problem()

        before_basis = trial_a.basis.copy()
        before_residual = (
            trial_a.mechanical_residual.copy()
        )

        result = enrich_tower_pgd_basis_once(
            trial_a,
            time,
            forcing,
            trial_a.mechanical_residual,
            H_sigma,
            metric,
            operator,
            mode_significance_tolerance=1.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-8,
            max_fixed_point_iterations=200,
            spatial_strategy="residual_ls",
        )

        self.assertFalse(result.accepted)
        self.assertIsNone(
            result.candidate_fixed_basis_result
        )
        self.assertEqual(
            trial_a.basis.n_modes,
            before_basis.n_modes,
        )
        np.testing.assert_allclose(
            trial_a.basis.plastic_strain_correction(),
            before_basis.plastic_strain_correction(),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            trial_a.mechanical_residual,
            before_residual,
            rtol=0.0,
            atol=0.0,
        )

    def test_unknown_spatial_strategy_is_rejected(self):
        (
            metric,
            operator,
            time,
            H_sigma,
            forcing,
            trial_a,
        ) = two_mode_problem()

        with self.assertRaises(ValueError):
            enrich_tower_pgd_basis_once(
                trial_a,
                time,
                forcing,
                trial_a.mechanical_residual,
                H_sigma,
                metric,
                operator,
                mode_significance_tolerance=0.0,
                acceptance_tolerance=0.0,
                spatial_strategy="not_a_strategy",
            )


if __name__ == "__main__":
    unittest.main()

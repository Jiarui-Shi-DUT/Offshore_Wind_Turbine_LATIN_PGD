# -*- coding: utf-8 -*-
"""Regression tests for tower one-mode PGD enrichment."""

import unittest
import numpy as np

from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import MaterialPointMetric, TowerEquilibriumOperator
from latin.tower_pgd_enrichment import _post_fixed_point_transform, enrich_tower_pgd_basis_once
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_state import MaterialPointLayout


def make_operator():
    layout = MaterialPointLayout(1, 2, 2)
    metric = MaterialPointMetric(np.array([1.0, 2.0, 1.5, 0.5]))
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
        reference_modulus=np.array([10.0, 20.0, 15.0, 25.0]),
        compatibility_matrix=H,
        free_dofs=np.array([1, 3]),
        n_dof=4,
    )
    return metric, operator


def be_rate(amplitude, time):
    rate = np.zeros_like(amplitude)
    rate[0] = (amplitude[1] - amplitude[0]) / (time[1] - time[0])
    rate[1:] = np.diff(amplitude) / np.diff(time)
    return rate


class TestPostFixedPointTransform(unittest.TestCase):
    def test_exact_transform_preserves_fields(self):
        metric, operator = make_operator()
        time = np.array([0.0, 0.2, 0.5, 0.9, 1.4])

        p1 = np.array([0.4, -0.2, 0.8, 0.1])
        p1 /= metric.norm(p1)
        s1 = operator.apply_spatial(p1).stress
        a1 = np.array([0.0, 0.1, 0.2, 0.15, 0.3])
        r1 = be_rate(a1, time)
        basis = PGDBasisTower(4, time.size, (PGDModeTower(p1, s1, a1, r1),))

        q = np.array([0.2, 0.7, -0.3, 0.6])
        q -= metric.inner_product(p1, q) * p1
        q /= metric.norm(q)
        raw_p = 0.6 * p1 + 0.8 * q
        raw_s = operator.apply_spatial(raw_p).stress
        raw_a = np.array([0.0, -0.04, 0.08, 0.13, 0.05])
        raw_r = be_rate(raw_a, time)
        raw = PGDModeTower(raw_p, raw_s, raw_a, raw_r)

        out = _post_fixed_point_transform(
            basis, raw, time, metric, operator,
            1.0e-14, 1.0e-12, 0.0, 1.0e-10, 1.0e-12, 2,
        )
        transformed, coeff, c, gamma_sp, gamma_lambda, ortho, ep, er, es = out
        self.assertEqual(transformed.n_modes, 2)
        self.assertAlmostEqual(coeff[0], 0.6, places=12)
        self.assertAlmostEqual(c, 0.8, places=12)
        self.assertAlmostEqual(gamma_sp, 0.8, places=12)
        self.assertGreater(gamma_lambda, 0.0)
        self.assertLess(max(ortho, ep, er, es), 1.0e-11)


class TestTowerPGDEnrichment(unittest.TestCase):
    def test_empty_basis_exact_rank_one_forcing_is_accepted(self):
        metric, operator = make_operator()
        time = np.array([0.0, 0.2, 0.5, 0.9, 1.4])

        p = np.array([0.6, -0.4, 0.9, -0.2])
        p /= metric.norm(p)
        s = operator.apply_spatial(p).stress
        amp = np.array([0.0, 0.08, 0.20, 0.14, 0.30])
        rate = be_rate(amp, time)
        Hs = 0.20 + 0.01*np.arange(time.size)[:,None] + 0.005*np.arange(4)[None,:]
        forcing = np.outer(rate, p) - Hs*np.outer(amp, s)

        empty = PGDBasisTower(4, time.size)
        trial_a = update_tower_pgd_time_functions(
            empty, time, forcing, Hs, metric, operator
        )
        result = enrich_tower_pgd_basis_once(
            trial_a, time, forcing, trial_a.mechanical_residual, Hs, metric, operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-8,
            max_fixed_point_iterations=100,
        )
        self.assertTrue(result.accepted, msg=result.failure_reason)
        self.assertEqual(result.n_modes, 1)
        self.assertTrue(result.fixed_point_converged)
        self.assertGreater(result.residual_benefit, 0.99)

    def test_rejection_does_not_mutate_trial_a(self):
        metric, operator = make_operator()
        time = np.array([0.0, 0.3, 0.7, 1.0])
        p = np.array([0.5, -0.1, 0.7, 0.2])
        p /= metric.norm(p)
        s = operator.apply_spatial(p).stress
        amp = np.array([0.0, 0.1, 0.2, 0.15])
        rate = be_rate(amp, time)
        Hs = np.full((time.size, 4), 0.25)
        forcing = np.outer(rate, p) - Hs*np.outer(amp, s)
        empty = PGDBasisTower(4, time.size)
        trial_a = update_tower_pgd_time_functions(empty, time, forcing, Hs, metric, operator)
        before = trial_a.basis.copy()

        result = enrich_tower_pgd_basis_once(
            trial_a, time, forcing, trial_a.mechanical_residual, Hs, metric, operator,
            mode_significance_tolerance=1.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-7,
            max_fixed_point_iterations=100,
        )
        self.assertFalse(result.accepted)
        self.assertIsNone(result.candidate_fixed_basis_result)
        self.assertEqual(trial_a.basis.n_modes, before.n_modes)
        np.testing.assert_allclose(
            trial_a.basis.plastic_strain_correction(),
            before.plastic_strain_correction(),
            rtol=0.0, atol=0.0,
        )

    def test_existing_basis_accepts_one_new_orthogonal_component(self):
        metric, operator = make_operator()
        time = np.array([0.0, 0.25, 0.55, 0.95, 1.35])

        p1 = np.array([0.5, -0.2, 0.8, 0.1])
        p1 /= metric.norm(p1)
        q = np.array([-0.3, 0.7, 0.1, 0.6])
        q -= metric.inner_product(p1, q) * p1
        p2 = q / metric.norm(q)
        s1 = operator.apply_spatial(p1).stress
        s2 = operator.apply_spatial(p2).stress

        a1 = np.array([0.0, 0.08, 0.18, 0.12, 0.25])
        a2 = np.array([0.0, -0.04, 0.06, 0.13, 0.05])
        r1, r2 = be_rate(a1, time), be_rate(a2, time)
        Hs = 0.18 + 0.008*np.arange(time.size)[:,None] + 0.004*np.arange(4)[None,:]
        forcing = (
            np.outer(r1, p1) - Hs*np.outer(a1, s1)
            + np.outer(r2, p2) - Hs*np.outer(a2, s2)
        )

        basis = PGDBasisTower(
            4, time.size,
            (PGDModeTower(p1, s1, np.zeros(time.size), np.zeros(time.size)),)
        )
        trial_a = update_tower_pgd_time_functions(basis, time, forcing, Hs, metric, operator)
        before = np.linalg.norm(trial_a.mechanical_residual)
        result = enrich_tower_pgd_basis_once(
            trial_a, time, forcing, trial_a.mechanical_residual, Hs, metric, operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-7,
            max_fixed_point_iterations=120,
        )
        self.assertTrue(result.accepted, msg=result.failure_reason)
        self.assertEqual(result.n_modes, 2)
        self.assertGreater(result.residual_benefit, 0.0)
        self.assertLess(
            np.linalg.norm(result.candidate_fixed_basis_result.mechanical_residual),
            before,
        )

    def test_integer_controls_reject_bool_and_fractional_values(self):
        metric, operator = make_operator()
        time = np.array([0.0, 0.5, 1.0])
        H_sigma = np.full((time.size, 4), 0.25)
        forcing = np.ones((time.size, 4), dtype=np.float64)
        empty = PGDBasisTower(n_material_points=4, n_time=time.size)
        trial_a = update_tower_pgd_time_functions(
            basis=empty,
            time=time,
            forcing=forcing,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
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
        )
        with self.assertRaises(TypeError):
            enrich_tower_pgd_basis_once(
                **common,
                iteration_added=True,
            )
        with self.assertRaises(TypeError):
            enrich_tower_pgd_basis_once(
                **common,
                max_fixed_point_iterations=3.5,
            )
        with self.assertRaises(TypeError):
            enrich_tower_pgd_basis_once(
                **common,
                reorthogonalization_passes=False,
            )


if __name__ == "__main__":
    unittest.main()

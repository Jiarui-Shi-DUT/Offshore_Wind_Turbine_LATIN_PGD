# -*- coding: utf-8 -*-
"""Tests for asymmetric multi-cycle nonlinear tower response."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_multicycle_response import (
    critical_plastic_strain_drifts,
    cycle_displacement_drifts,
    normalized_cycle_displacement_drifts,
    ratcheting_metrics,
    run_asymmetric_multicycle_tower_analysis,
)
from fem.tower_loading import AsymmetricCyclicTopForceHistory
from material.viscoplastic_damage_1d import MaterialParameters


class TestAsymmetricMulticycleTowerResponse(unittest.TestCase):
    """Run a small two-cycle asymmetric finite-element analysis."""

    @classmethod
    def setUpClass(cls) -> None:
        configuration = TowerConfiguration(
            horizontal_force=1.0e6,
            n_elements=4,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        material = MaterialParameters()

        cls.result = run_asymmetric_multicycle_tower_analysis(
            configuration=configuration,
            material=material,
            maximum_force=1.0e6,
            force_ratio=-0.5,
            period=10.0,
            n_cycles=2,
            increments_per_cycle=8,
            similarity_tolerance=1.0e-3,
            max_iterations=40,
        )

    def test_loading_is_asymmetric_and_sign_reversing(self) -> None:
        loading = self.result.response.loading

        self.assertIsInstance(
            loading,
            AsymmetricCyclicTopForceHistory,
        )
        self.assertAlmostEqual(loading.force_ratio, -0.5, places=14)
        self.assertGreater(loading.mean_force, 0.0)
        self.assertGreater(loading.maximum_force, 0.0)
        self.assertLess(loading.minimum_force, 0.0)

    def test_exact_asymmetric_cycle_checkpoints(self) -> None:
        loading = self.result.response.loading
        quarter = loading.increments_per_cycle // 4
        indices = np.array(
            [
                0,
                quarter,
                2 * quarter,
                3 * quarter,
                4 * quarter,
                5 * quarter,
                6 * quarter,
                7 * quarter,
                8 * quarter,
            ],
            dtype=np.int64,
        )
        expected = np.array(
            [
                loading.mean_force,
                loading.maximum_force,
                loading.mean_force,
                loading.minimum_force,
                loading.mean_force,
                loading.maximum_force,
                loading.mean_force,
                loading.minimum_force,
                loading.mean_force,
            ],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            loading.forces[indices],
            expected,
            rtol=0.0,
            atol=1.0e-8,
        )

    def test_preload_is_not_counted_as_a_periodic_cycle(self) -> None:
        response = self.result.response
        loading = response.loading

        self.assertAlmostEqual(
            float(response.analysis_times[0]),
            0.25 * loading.period,
            places=14,
        )
        self.assertEqual(self.result.diagnostics.n_cycles, 2)
        self.assertEqual(
            response.top_displacements.shape,
            (loading.n_time_points,),
        )

    def test_cycle_displacement_drift_matches_same_force_endpoints(
        self,
    ) -> None:
        drifts = cycle_displacement_drifts(self.result)
        response = self.result.response

        self.assertEqual(drifts.shape, (2,))
        self.assertTrue(np.all(np.isfinite(drifts)))

        for index, cycle in enumerate(self.result.diagnostics.cycles):
            expected = (
                response.top_displacements[cycle.indices.end]
                - response.top_displacements[cycle.indices.start]
            )
            self.assertAlmostEqual(drifts[index], expected, places=14)

    def test_normalized_global_drift_is_finite_and_nonnegative(
        self,
    ) -> None:
        values = normalized_cycle_displacement_drifts(self.result)
        self.assertEqual(values.shape, (2,))
        self.assertTrue(np.all(np.isfinite(values)))
        self.assertTrue(np.all(values >= 0.0))

    def test_plastic_strain_drift_matches_cycle_diagnostics(
        self,
    ) -> None:
        values = critical_plastic_strain_drifts(self.result)
        np.testing.assert_allclose(
            values,
            self.result.diagnostics.critical_plastic_strain_increments,
            rtol=0.0,
            atol=0.0,
        )

    def test_ratcheting_metrics_have_one_value_per_cycle(self) -> None:
        global_drift, normalized_drift, plastic_drift = (
            ratcheting_metrics(self.result)
        )
        for values in (global_drift, normalized_drift, plastic_drift):
            self.assertEqual(values.shape, (2,))
            self.assertTrue(np.all(np.isfinite(values)))

    def test_damage_is_non_decreasing(self) -> None:
        response = self.result.response
        diagnostics = self.result.diagnostics

        self.assertTrue(
            np.all(np.diff(response.maximum_damages) >= -1.0e-14)
        )
        self.assertTrue(
            np.all(diagnostics.maximum_damage_increments >= -1.0e-14)
        )

    def test_solver_convergence_metrics(self) -> None:
        response = self.result.response

        self.assertTrue(np.all(response.iterations >= 1))
        self.assertLessEqual(int(np.max(response.iterations)), 40)
        self.assertLess(float(np.max(response.residual_norms)), 1.0e-1)

    def test_similarity_comparison_covers_adjacent_cycles(self) -> None:
        self.assertEqual(len(self.result.similarities), 1)
        similarity = self.result.similarities[0]
        self.assertEqual(similarity.previous_cycle, 1)
        self.assertEqual(similarity.current_cycle, 2)
        self.assertTrue(np.isfinite(similarity.maximum_indicator))


if __name__ == "__main__":
    unittest.main()

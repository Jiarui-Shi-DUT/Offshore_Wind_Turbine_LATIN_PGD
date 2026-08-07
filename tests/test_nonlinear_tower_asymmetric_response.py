# -*- coding: utf-8 -*-
"""Regression tests for asymmetric nonlinear tower cyclic response."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_response import (
    run_nonlinear_asymmetric_analysis,
)
from fem.tower_loading import (
    AsymmetricCyclicTopForceHistory,
    create_asymmetric_cyclic_top_force_history,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestNonlinearTowerAsymmetricResponse(unittest.TestCase):
    """Run a reduced asymmetric tower cycle with explicit preload."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.configuration = TowerConfiguration(
            horizontal_force=1.0e6,
            n_elements=4,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        cls.material = MaterialParameters()
        cls.loading = create_asymmetric_cyclic_top_force_history(
            maximum_force=1.0e6,
            force_ratio=-0.5,
            period=10.0,
            n_cycles=1,
            increments_per_cycle=8,
        )
        cls.response = run_nonlinear_asymmetric_analysis(
            configuration=cls.configuration,
            material=cls.material,
            loading=cls.loading,
            max_iterations=40,
        )

    def test_response_uses_asymmetric_loading(self) -> None:
        self.assertIsInstance(
            self.response.loading,
            AsymmetricCyclicTopForceHistory,
        )
        self.assertAlmostEqual(
            self.response.loading.force_ratio,
            -0.5,
            places=14,
        )

    def test_preload_ends_at_periodic_mean_state(self) -> None:
        expected_preload_duration = 0.25 * self.loading.period
        self.assertAlmostEqual(
            float(self.response.analysis_times[0]),
            expected_preload_duration,
            places=14,
        )
        self.assertAlmostEqual(
            float(self.loading.forces[0]),
            self.loading.mean_force,
            places=8,
        )
        self.assertGreater(
            abs(float(self.response.top_displacements[0])),
            1.0e-8,
        )

    def test_periodic_response_contains_only_loading_points(self) -> None:
        n_points = self.loading.n_time_points
        self.assertEqual(
            self.response.top_displacements.shape,
            (n_points,),
        )
        self.assertEqual(
            self.response.fiber_strains.shape[0],
            n_points,
        )
        self.assertAlmostEqual(
            float(self.response.analysis_times[-1]),
            0.25 * self.loading.period
            + self.loading.n_cycles * self.loading.period,
            places=14,
        )

    def test_force_history_is_sign_reversing_and_asymmetric(self) -> None:
        quarter = self.loading.increments_per_cycle // 4
        forces = self.loading.forces

        self.assertAlmostEqual(
            float(forces[quarter]),
            self.loading.maximum_force,
            places=8,
        )
        self.assertAlmostEqual(
            float(forces[3 * quarter]),
            self.loading.minimum_force,
            places=8,
        )
        self.assertGreater(
            self.loading.maximum_force,
            abs(self.loading.minimum_force),
        )
        self.assertGreater(
            self.loading.mean_force,
            0.0,
        )

    def test_complete_internal_states_are_finite(self) -> None:
        response = self.response
        self.assertEqual(
            response.fiber_states.shape,
            response.fiber_strains.shape + (4,),
        )
        self.assertTrue(
            np.all(np.isfinite(response.top_displacements))
        )
        self.assertTrue(
            np.all(np.isfinite(response.fiber_stresses))
        )
        self.assertTrue(
            np.all(np.isfinite(response.fiber_states))
        )

    def test_newton_solution_is_converged(self) -> None:
        response = self.response
        self.assertTrue(np.all(response.iterations >= 1))
        self.assertTrue(np.all(response.iterations <= 40))
        self.assertLess(
            float(np.max(response.residual_norms)),
            1.0e-1,
        )


if __name__ == "__main__":
    unittest.main()

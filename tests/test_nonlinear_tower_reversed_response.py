# -*- coding: utf-8 -*-
"""Regression tests for the fully reversed nonlinear tower example."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_reversed_response import (
    locate_critical_fiber,
    run_nonlinear_reversed_analysis,
)
from fem.tower_loading import create_reversed_top_force_history
from material.viscoplastic_damage_1d import MaterialParameters


class TestReversedCriticalFiberSelection(unittest.TestCase):
    """Test peak-history selection for the reversed example."""

    def test_peak_damage_is_primary_selection_criterion(self) -> None:
        states = np.zeros((4, 2, 2, 3, 4), dtype=np.float64)

        states[1, 1, 0, 2, 3] = 0.03
        states[-1, 0, 1, 1, 0] = 0.5

        location = locate_critical_fiber(
            fiber_states=states,
        )

        self.assertEqual(location, (1, 0, 2))

    def test_peak_plastic_strain_is_fallback_without_damage(self) -> None:
        states = np.zeros((4, 2, 2, 3, 4), dtype=np.float64)

        states[2, 0, 1, 1, 0] = -0.5

        location = locate_critical_fiber(
            fiber_states=states,
        )

        self.assertEqual(location, (0, 1, 1))


class TestNonlinearTowerReversedResponse(unittest.TestCase):
    """Run a reduced fully reversed tower cycle as a regression test."""

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
        cls.loading = create_reversed_top_force_history(
            force_amplitude=1.0e6,
            period=10.0,
            n_cycles=1,
            increments_per_cycle=8,
        )
        cls.response = run_nonlinear_reversed_analysis(
            configuration=cls.configuration,
            material=cls.material,
            loading=cls.loading,
            max_iterations=40,
        )

    def test_history_shapes_and_complete_states(self) -> None:
        response = self.response
        n_points = self.loading.n_time_points

        self.assertEqual(
            response.top_displacements.shape,
            (n_points,),
        )
        self.assertEqual(
            response.iterations.shape,
            (n_points,),
        )
        self.assertEqual(
            response.fiber_strains.shape,
            response.fiber_stresses.shape,
        )
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
        self.assertTrue(
            np.all(
                np.isfinite(
                    response.critical_fiber_yield_functions
                )
            )
        )

    def test_fixed_fiber_stress_reverses_between_peaks(self) -> None:
        quarter = self.loading.increments_per_cycle // 4
        positive_peak = quarter
        negative_peak = 3 * quarter
        stresses = self.response.critical_fiber_stresses

        self.assertLess(
            float(
                stresses[positive_peak]
                * stresses[negative_peak]
            ),
            0.0,
        )
        self.assertGreater(
            abs(float(stresses[positive_peak])),
            self.material.sigma_y,
        )
        self.assertGreater(
            abs(float(stresses[negative_peak])),
            self.material.sigma_y,
        )

    def test_plastic_flow_and_backstress_reverse(self) -> None:
        quarter = self.loading.increments_per_cycle // 4
        plastic_strains = (
            self.response.critical_fiber_plastic_strains
        )
        backstresses = (
            self.response.critical_fiber_backstresses
        )

        positive_half_increment = float(
            plastic_strains[quarter]
            - plastic_strains[0]
        )
        negative_half_increment = float(
            plastic_strains[3 * quarter]
            - plastic_strains[2 * quarter]
        )

        self.assertLess(
            positive_half_increment
            * negative_half_increment,
            0.0,
        )
        self.assertLess(
            float(
                backstresses[quarter]
                * backstresses[3 * quarter]
            ),
            0.0,
        )

    def test_cycle_develops_irreversible_state(self) -> None:
        response = self.response

        self.assertGreater(
            float(np.max(response.maximum_absolute_plastic_strains)),
            1.0e-8,
        )
        self.assertGreater(
            float(np.max(response.maximum_damages)),
            1.0e-6,
        )
        self.assertTrue(
            np.all(
                np.diff(response.fiber_damages, axis=0)
                >= -1.0e-14
            )
        )

    def test_zero_force_cycle_leaves_residual_displacement(self) -> None:
        response = self.response

        self.assertAlmostEqual(
            float(self.loading.forces[0]),
            0.0,
            places=8,
        )
        self.assertAlmostEqual(
            float(self.loading.forces[-1]),
            0.0,
            places=8,
        )
        self.assertGreater(
            abs(float(response.top_displacements[-1])),
            1.0e-3,
        )

    def test_critical_fiber_is_near_tower_base(self) -> None:
        response = self.response

        self.assertLess(
            response.critical_height,
            0.25 * self.configuration.height,
        )
        self.assertGreater(
            abs(response.critical_y_coordinate),
            0.9,
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

# -*- coding: utf-8 -*-
"""Regression tests for nonlinear tower diagnostic examples."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_load_probe import parse_force_levels
from examples.nonlinear_tower_pulsating_response import (
    locate_final_critical_fiber,
    run_nonlinear_pulsating_analysis,
)
from fem.tower_loading import create_pulsating_top_force_history
from material.viscoplastic_damage_1d import MaterialParameters


class TestNonlinearTowerLoadProbe(unittest.TestCase):
    """Test force-level parsing used by the monotonic probe."""

    def test_parse_force_levels_converts_meganewtons(self) -> None:
        forces = parse_force_levels("0.2, 0.4, 0.8")

        np.testing.assert_allclose(
            forces,
            (0.2e6, 0.4e6, 0.8e6),
            rtol=0.0,
            atol=0.0,
        )

    def test_parse_force_levels_rejects_invalid_sequences(self) -> None:
        invalid_inputs = (
            "",
            "0.2,",
            "0.4,0.2",
            "0.2,0.2",
            "-0.2,0.4",
            "not-a-number",
        )

        for value in invalid_inputs:
            with self.subTest(value=value):
                with self.assertRaises((TypeError, ValueError)):
                    parse_force_levels(value)


class TestCriticalFiberSelection(unittest.TestCase):
    """Test selection of one fixed critical fiber."""

    def test_final_damage_is_primary_selection_criterion(self) -> None:
        damages = np.zeros((3, 2, 2, 4), dtype=np.float64)
        plastic_strains = np.zeros_like(damages)

        damages[-1, 1, 0, 3] = 0.02
        plastic_strains[-1, 0, 1, 2] = 0.5

        location = locate_final_critical_fiber(
            fiber_damages=damages,
            fiber_plastic_strains=plastic_strains,
        )

        self.assertEqual(location, (1, 0, 3))

    def test_plastic_strain_is_fallback_when_damage_is_zero(self) -> None:
        damages = np.zeros((3, 2, 2, 4), dtype=np.float64)
        plastic_strains = np.zeros_like(damages)

        plastic_strains[-1, 0, 1, 2] = -0.5

        location = locate_final_critical_fiber(
            fiber_damages=damages,
            fiber_plastic_strains=plastic_strains,
        )

        self.assertEqual(location, (0, 1, 2))


class TestNonlinearTowerPulsatingResponse(unittest.TestCase):
    """Run a reduced nonlinear tower cycle as a regression test."""

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
        cls.loading = create_pulsating_top_force_history(
            maximum_force=1.0e6,
            force_ratio=0.1,
            period=10.0,
            n_cycles=1,
            increments_per_cycle=8,
        )
        cls.response = run_nonlinear_pulsating_analysis(
            configuration=cls.configuration,
            material=cls.material,
            loading=cls.loading,
            max_iterations=40,
        )

    def test_history_shapes_and_finite_values(self) -> None:
        response = self.response
        n_points = self.loading.n_time_points

        self.assertEqual(response.top_displacements.shape, (n_points,))
        self.assertEqual(response.iterations.shape, (n_points,))
        self.assertEqual(response.fiber_strains.shape[0], n_points)
        self.assertEqual(
            response.fiber_strains.shape,
            response.fiber_stresses.shape,
        )
        self.assertEqual(
            response.fiber_strains.shape,
            response.fiber_plastic_strains.shape,
        )
        self.assertEqual(
            response.fiber_strains.shape,
            response.fiber_damages.shape,
        )

        self.assertTrue(
            np.all(np.isfinite(response.top_displacements))
        )
        self.assertTrue(
            np.all(np.isfinite(response.fiber_stresses))
        )
        self.assertTrue(
            np.all(np.isfinite(response.fiber_plastic_strains))
        )
        self.assertTrue(
            np.all(np.isfinite(response.fiber_damages))
        )

    def test_cycle_develops_plasticity_and_damage(self) -> None:
        response = self.response

        self.assertGreater(
            float(np.max(response.maximum_absolute_plastic_strains)),
            1.0e-8,
        )
        self.assertGreater(
            float(np.max(response.maximum_damages)),
            1.0e-6,
        )
        self.assertGreater(
            response.maximum_damages[-1],
            response.maximum_damages[0],
        )

    def test_same_force_returns_with_residual_displacement(self) -> None:
        response = self.response

        self.assertAlmostEqual(
            float(self.loading.forces[0]),
            float(self.loading.forces[-1]),
            places=8,
        )
        self.assertGreater(
            float(
                response.top_displacements[-1]
                - response.top_displacements[0]
            ),
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

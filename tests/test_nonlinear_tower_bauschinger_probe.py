# -*- coding: utf-8 -*-
"""Regression tests for the tower Bauschinger-effect probe."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_bauschinger_probe import (
    directional_reverse_yield_limits,
    evaluate_bauschinger_metrics,
    find_yield_onset,
    interpolate_zero_crossing,
)
from examples.nonlinear_tower_reversed_response import (
    run_nonlinear_reversed_analysis,
)
from fem.tower_loading import create_reversed_top_force_history
from material.viscoplastic_damage_1d import MaterialParameters


class TestDirectionalReverseYieldLimits(unittest.TestCase):
    """Test the constitutive reverse-yield limit calculation."""

    def test_virgin_state_returns_initial_yield_stress(self) -> None:
        material = MaterialParameters()
        state = np.zeros(4, dtype=np.float64)

        (
            beta,
            isotropic_force,
            damage,
            limit_with_kinematic,
            limit_without_kinematic,
        ) = directional_reverse_yield_limits(
            material=material,
            state=state,
        )

        self.assertEqual(beta, 0.0)
        self.assertEqual(isotropic_force, 0.0)
        self.assertEqual(damage, 0.0)
        self.assertAlmostEqual(
            limit_with_kinematic,
            material.sigma_y,
            places=12,
        )
        self.assertAlmostEqual(
            limit_without_kinematic,
            material.sigma_y,
            places=12,
        )

    def test_negative_backstress_reduces_positive_reverse_limit(self) -> None:
        material = MaterialParameters()
        state = np.array(
            [
                0.0,
                -1.0e-4,
                1.0e-3,
                2.0e-3,
            ],
            dtype=np.float64,
        )

        (
            beta,
            _isotropic_force,
            _damage,
            limit_with_kinematic,
            limit_without_kinematic,
        ) = directional_reverse_yield_limits(
            material=material,
            state=state,
        )

        self.assertLess(beta, 0.0)
        self.assertLess(
            limit_with_kinematic,
            limit_without_kinematic,
        )

    def test_invalid_state_is_rejected(self) -> None:
        material = MaterialParameters()

        invalid_states = (
            np.zeros(3, dtype=np.float64),
            np.array(
                [0.0, 0.0, np.nan, 0.0],
                dtype=np.float64,
            ),
        )

        for state in invalid_states:
            with self.subTest(shape=state.shape):
                with self.assertRaises(ValueError):
                    directional_reverse_yield_limits(
                        material=material,
                        state=state,
                    )


class TestTowerBauschingerProbe(unittest.TestCase):
    """Run a reduced tower cycle and verify yield translation."""

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
            increments_per_cycle=40,
        )
        cls.response = run_nonlinear_reversed_analysis(
            configuration=cls.configuration,
            material=cls.material,
            loading=cls.loading,
            max_iterations=40,
        )
        cls.metrics = evaluate_bauschinger_metrics(
            cls.response
        )

    def test_initial_and_reverse_yield_crossings_are_found(self) -> None:
        quarter = self.loading.increments_per_cycle // 4
        half = 2 * quarter
        three_quarters = 3 * quarter

        initial = find_yield_onset(
            response=self.response,
            start_index=0,
            end_index=quarter,
        )
        reverse = find_yield_onset(
            response=self.response,
            start_index=half,
            end_index=three_quarters,
        )

        self.assertGreater(initial.normalized_time, 0.0)
        self.assertLess(initial.normalized_time, 0.25)
        self.assertGreater(reverse.normalized_time, 0.5)
        self.assertLess(reverse.normalized_time, 0.75)
        self.assertLess(initial.stress, 0.0)
        self.assertGreater(reverse.stress, 0.0)
        self.assertEqual(initial.yield_function, 0.0)
        self.assertEqual(reverse.yield_function, 0.0)

    def test_initial_yield_matches_virgin_yield_stress(self) -> None:
        self.assertAlmostEqual(
            self.metrics.observed_initial_yield_magnitude,
            self.material.sigma_y,
            delta=1.0e-2,
        )

    def test_reverse_yield_is_lower_than_initial_yield(self) -> None:
        self.assertGreater(
            self.metrics.observed_change_from_initial,
            0.0,
        )
        self.assertLess(
            self.metrics.observed_reverse_yield_magnitude,
            self.metrics.observed_initial_yield_magnitude,
        )

    def test_backstress_causes_positive_kinematic_reduction(self) -> None:
        self.assertLess(
            self.metrics.reversal_backstress,
            0.0,
        )
        self.assertGreater(
            self.metrics.kinematic_reverse_yield_reduction,
            0.0,
        )
        self.assertGreater(
            self.metrics.kinematic_reduction_ratio,
            0.0,
        )
        self.assertLess(
            self.metrics.kinematic_reduction_ratio,
            0.1,
        )

    def test_observed_reverse_onset_matches_translated_limit(self) -> None:
        self.assertAlmostEqual(
            self.metrics.observed_reverse_yield_magnitude,
            self.metrics.reverse_limit_with_kinematic,
            delta=1.0e-2,
        )

    def test_counterfactual_limit_exceeds_translated_limit(self) -> None:
        self.assertGreater(
            self.metrics.reverse_limit_without_kinematic,
            self.metrics.reverse_limit_with_kinematic,
        )
        self.assertAlmostEqual(
            self.metrics.kinematic_reverse_yield_reduction,
            (
                self.metrics.reverse_limit_without_kinematic
                - self.metrics.reverse_limit_with_kinematic
            ),
            places=12,
        )

    def test_interpolation_rejects_non_bracketing_indices(self) -> None:
        with self.assertRaises(ValueError):
            interpolate_zero_crossing(
                response=self.response,
                lower_index=0,
                upper_index=2,
            )

        with self.assertRaises(ValueError):
            interpolate_zero_crossing(
                response=self.response,
                lower_index=0,
                upper_index=1,
            )


if __name__ == "__main__":
    unittest.main()

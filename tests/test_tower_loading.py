# -*- coding: utf-8 -*-
"""Regression tests for offshore-wind-turbine tower loading histories."""

import unittest

import numpy as np

from fem.tower_loading import (
    create_pulsating_top_force_history,
    create_reversed_top_force_history,
    evaluate_pulsating_top_force,
    evaluate_reversed_top_force,
)


class TestTowerLoading(unittest.TestCase):
    """Verify the positive-mean sinusoidal tower-top force definition."""

    def setUp(self) -> None:
        self.maximum_force = 1.0e6
        self.force_ratio = 0.1
        self.period = 2.0

    def test_force_parameters_for_ratio_point_one(self) -> None:
        """R_F = 0.1 must give the prescribed mean and amplitude."""
        history = create_pulsating_top_force_history(
            maximum_force=self.maximum_force,
            force_ratio=self.force_ratio,
            period=self.period,
            n_cycles=1,
            increments_per_cycle=100,
        )

        self.assertAlmostEqual(
            history.minimum_force,
            0.1 * self.maximum_force,
            places=8,
        )
        self.assertAlmostEqual(
            history.mean_force,
            0.55 * self.maximum_force,
            places=8,
        )
        self.assertAlmostEqual(
            history.force_amplitude,
            0.45 * self.maximum_force,
            places=8,
        )
        self.assertAlmostEqual(
            history.angular_frequency,
            np.pi,
            places=14,
        )
        self.assertAlmostEqual(
            history.time_increment,
            self.period / 100.0,
            places=14,
        )

    def test_quarter_cycle_values_are_exact(self) -> None:
        """The discrete history must contain mean, maximum, and minimum."""
        history = create_pulsating_top_force_history(
            maximum_force=self.maximum_force,
            force_ratio=self.force_ratio,
            period=self.period,
            n_cycles=1,
            increments_per_cycle=100,
        )

        quarter = history.increments_per_cycle // 4
        expected = np.array(
            [
                history.mean_force,
                history.maximum_force,
                history.mean_force,
                history.minimum_force,
                history.mean_force,
            ],
            dtype=np.float64,
        )
        actual = history.forces[
            np.array(
                [
                    0,
                    quarter,
                    2 * quarter,
                    3 * quarter,
                    4 * quarter,
                ],
                dtype=np.int64,
            )
        ]

        np.testing.assert_allclose(
            actual,
            expected,
            rtol=0.0,
            atol=1.0e-9,
        )
        self.assertAlmostEqual(
            float(np.max(history.forces)),
            history.maximum_force,
            places=8,
        )
        self.assertAlmostEqual(
            float(np.min(history.forces)),
            history.minimum_force,
            places=8,
        )

    def test_multiple_cycles_have_correct_size_and_periodicity(self) -> None:
        """A multi-cycle history must repeat and include the final endpoint."""
        n_cycles = 3
        increments_per_cycle = 40
        history = create_pulsating_top_force_history(
            maximum_force=self.maximum_force,
            force_ratio=self.force_ratio,
            period=self.period,
            n_cycles=n_cycles,
            increments_per_cycle=increments_per_cycle,
        )

        self.assertEqual(
            history.n_time_points,
            n_cycles * increments_per_cycle + 1,
        )
        self.assertEqual(
            history.times.shape,
            history.forces.shape,
        )
        self.assertAlmostEqual(
            history.times[0],
            0.0,
            places=14,
        )
        self.assertAlmostEqual(
            history.times[-1],
            n_cycles * self.period,
            places=14,
        )

        first_cycle = history.forces[: increments_per_cycle + 1]
        second_cycle = history.forces[
            increments_per_cycle:
            2 * increments_per_cycle + 1
        ]
        np.testing.assert_allclose(
            first_cycle,
            second_cycle,
            rtol=0.0,
            atol=1.0e-9,
        )

    def test_direct_evaluation_matches_history(self) -> None:
        """Direct force evaluation and generated history must be identical."""
        history = create_pulsating_top_force_history(
            maximum_force=self.maximum_force,
            force_ratio=self.force_ratio,
            period=self.period,
            n_cycles=2,
            increments_per_cycle=40,
        )
        evaluated = evaluate_pulsating_top_force(
            time=history.times,
            maximum_force=self.maximum_force,
            force_ratio=self.force_ratio,
            period=self.period,
        )

        np.testing.assert_allclose(
            evaluated,
            history.forces,
            rtol=0.0,
            atol=0.0,
        )

    def test_unit_force_ratio_produces_constant_force(self) -> None:
        """R_F = 1 must reduce the sinusoid to a constant force."""
        history = create_pulsating_top_force_history(
            maximum_force=self.maximum_force,
            force_ratio=1.0,
            period=self.period,
            n_cycles=2,
            increments_per_cycle=20,
        )

        np.testing.assert_allclose(
            history.forces,
            self.maximum_force,
            rtol=0.0,
            atol=0.0,
        )
        self.assertAlmostEqual(
            history.force_amplitude,
            0.0,
            places=14,
        )

    def test_invalid_loading_parameters_are_rejected(self) -> None:
        """Physically or numerically invalid inputs must raise errors."""
        with self.assertRaises(ValueError):
            create_pulsating_top_force_history(
                maximum_force=0.0,
            )

        with self.assertRaises(ValueError):
            create_pulsating_top_force_history(
                maximum_force=self.maximum_force,
                force_ratio=-0.1,
            )

        with self.assertRaises(ValueError):
            create_pulsating_top_force_history(
                maximum_force=self.maximum_force,
                force_ratio=1.1,
            )

        with self.assertRaises(ValueError):
            create_pulsating_top_force_history(
                maximum_force=self.maximum_force,
                period=0.0,
            )

        with self.assertRaises(ValueError):
            create_pulsating_top_force_history(
                maximum_force=self.maximum_force,
                n_cycles=0,
            )

        with self.assertRaises(ValueError):
            create_pulsating_top_force_history(
                maximum_force=self.maximum_force,
                increments_per_cycle=10,
            )

        with self.assertRaises(TypeError):
            create_pulsating_top_force_history(
                maximum_force=self.maximum_force,
                n_cycles=1.5,
            )


class TestReversedTowerLoading(unittest.TestCase):
    """Verify the zero-mean fully reversed tower-top force."""

    def setUp(self) -> None:
        self.force_amplitude = 1.0e6
        self.period = 2.0

    def test_reversed_force_parameters(self) -> None:
        """The history must have equal positive and negative peaks."""
        history = create_reversed_top_force_history(
            force_amplitude=self.force_amplitude,
            period=self.period,
            n_cycles=1,
            increments_per_cycle=100,
        )

        self.assertAlmostEqual(
            history.maximum_force,
            self.force_amplitude,
            places=8,
        )
        self.assertAlmostEqual(
            history.minimum_force,
            -self.force_amplitude,
            places=8,
        )
        self.assertAlmostEqual(
            history.mean_force,
            0.0,
            places=14,
        )
        self.assertAlmostEqual(
            history.angular_frequency,
            np.pi,
            places=14,
        )
        self.assertAlmostEqual(
            history.time_increment,
            self.period / 100.0,
            places=14,
        )

    def test_quarter_cycle_values_are_exact(self) -> None:
        """The discrete cycle must contain both force extrema."""
        history = create_reversed_top_force_history(
            force_amplitude=self.force_amplitude,
            period=self.period,
            n_cycles=1,
            increments_per_cycle=100,
        )

        quarter = history.increments_per_cycle // 4
        indices = np.array(
            [
                0,
                quarter,
                2 * quarter,
                3 * quarter,
                4 * quarter,
            ],
            dtype=np.int64,
        )
        expected = np.array(
            [
                0.0,
                self.force_amplitude,
                0.0,
                -self.force_amplitude,
                0.0,
            ],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            history.forces[indices],
            expected,
            rtol=0.0,
            atol=1.0e-9,
        )

    def test_multiple_cycles_repeat_and_include_endpoint(self) -> None:
        """The generated history must be periodic over complete cycles."""
        n_cycles = 3
        increments_per_cycle = 40
        history = create_reversed_top_force_history(
            force_amplitude=self.force_amplitude,
            period=self.period,
            n_cycles=n_cycles,
            increments_per_cycle=increments_per_cycle,
        )

        self.assertEqual(
            history.n_time_points,
            n_cycles * increments_per_cycle + 1,
        )
        self.assertAlmostEqual(
            history.times[-1],
            n_cycles * self.period,
            places=14,
        )

        first_cycle = history.forces[: increments_per_cycle + 1]
        second_cycle = history.forces[
            increments_per_cycle:
            2 * increments_per_cycle + 1
        ]
        np.testing.assert_allclose(
            first_cycle,
            second_cycle,
            rtol=32.0 * np.finfo(np.float64).eps,
            atol=1.0e-9,
        )

    def test_direct_evaluation_matches_history(self) -> None:
        """Direct evaluation and generated history must be identical."""
        history = create_reversed_top_force_history(
            force_amplitude=self.force_amplitude,
            period=self.period,
            n_cycles=2,
            increments_per_cycle=40,
        )
        evaluated = evaluate_reversed_top_force(
            time=history.times,
            force_amplitude=self.force_amplitude,
            period=self.period,
        )

        np.testing.assert_allclose(
            evaluated,
            history.forces,
            rtol=0.0,
            atol=0.0,
        )

    def test_invalid_reversed_loading_parameters_are_rejected(self) -> None:
        """Invalid amplitudes and discretisations must raise errors."""
        with self.assertRaises(ValueError):
            create_reversed_top_force_history(
                force_amplitude=0.0,
            )

        with self.assertRaises(ValueError):
            create_reversed_top_force_history(
                force_amplitude=-1.0,
            )

        with self.assertRaises(ValueError):
            create_reversed_top_force_history(
                force_amplitude=self.force_amplitude,
                period=0.0,
            )

        with self.assertRaises(ValueError):
            create_reversed_top_force_history(
                force_amplitude=self.force_amplitude,
                n_cycles=0,
            )

        with self.assertRaises(ValueError):
            create_reversed_top_force_history(
                force_amplitude=self.force_amplitude,
                increments_per_cycle=10,
            )

        with self.assertRaises(TypeError):
            create_reversed_top_force_history(
                force_amplitude=self.force_amplitude,
                n_cycles=1.5,
            )


if __name__ == "__main__":
    unittest.main()

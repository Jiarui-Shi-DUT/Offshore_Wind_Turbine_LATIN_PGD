# -*- coding: utf-8 -*-
"""Focused policy tests for OPT-6 residual-LS spatial Aitken acceleration."""

import unittest

import numpy as np

from latin.tower_equilibrium_operator import MaterialPointMetric
from latin.tower_pgd_enrichment import (
    _AITKEN_ACTIVATION_THRESHOLD,
    _AITKEN_MONOTONE_WINDOW,
    _AITKEN_OMEGA_MAX,
    _aitken_relaxation_factor,
    _aitken_tail_ready,
    _align_spatial_sign,
    _residual_ls_fixed_point_converged,
)


class TestTowerPGDOPT6AitkenPolicy(unittest.TestCase):
    def setUp(self) -> None:
        self.metric = MaterialPointMetric(
            np.array([1.0, 2.0, 0.5], dtype=np.float64)
        )

    def test_tail_gate_requires_small_strictly_decreasing_window(self) -> None:
        threshold = _AITKEN_ACTIVATION_THRESHOLD
        window = _AITKEN_MONOTONE_WINDOW

        self.assertFalse(
            _aitken_tail_ready([0.8 * threshold] * (window - 1))
        )

        not_small = np.linspace(
            2.0 * threshold,
            1.1 * threshold,
            window,
        ).tolist()
        self.assertFalse(_aitken_tail_ready(not_small))

        self.assertFalse(
            _aitken_tail_ready(
                [8.0e-3, 6.0e-3, 6.5e-3, 4.0e-3, 3.0e-3]
            )
        )
        self.assertTrue(
            _aitken_tail_ready(
                [9.0e-3, 7.0e-3, 5.0e-3, 3.0e-3, 2.0e-3]
            )
        )

    def test_vector_aitken_clips_to_calibrated_acceleration_range(self) -> None:
        previous = np.array(
            [1.0, -0.5, 0.25],
            dtype=np.float64,
        )
        current = 0.8 * previous

        omega = _aitken_relaxation_factor(
            1.0,
            previous,
            current,
            self.metric,
        )
        self.assertAlmostEqual(omega, _AITKEN_OMEGA_MAX)

        damping_candidate = _aitken_relaxation_factor(
            1.0,
            previous,
            -previous,
            self.metric,
        )
        self.assertEqual(damping_candidate, 1.0)

    def test_sign_alignment_removes_only_pgd_plus_minus_gauge(self) -> None:
        reference = np.array(
            [0.5, -0.2, 0.9],
            dtype=np.float64,
        )
        candidate = -reference
        aligned = _align_spatial_sign(
            reference,
            candidate,
            self.metric,
        )
        np.testing.assert_allclose(
            aligned,
            reference,
            rtol=0.0,
            atol=0.0,
        )

    def test_dual_gate_rejects_false_accelerated_convergence(self) -> None:
        tolerance = 1.0e-5
        self.assertFalse(
            _residual_ls_fixed_point_converged(
                1.0e-6,
                2.0e-5,
                tolerance,
            )
        )
        self.assertFalse(
            _residual_ls_fixed_point_converged(
                2.0e-5,
                1.0e-6,
                tolerance,
            )
        )
        self.assertTrue(
            _residual_ls_fixed_point_converged(
                1.0e-6,
                2.0e-6,
                tolerance,
            )
        )


if __name__ == "__main__":
    unittest.main()

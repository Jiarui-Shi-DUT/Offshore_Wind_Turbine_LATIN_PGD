# -*- coding: utf-8 -*-
"""Regression tests for the elastic tower pulsating-response example."""

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.elastic_tower_pulsating_response import (
    build_scaled_history,
    verify_direct_checkpoint_solutions,
)
from fem.tower_loading import create_pulsating_top_force_history


class TestElasticTowerPulsatingResponse(unittest.TestCase):
    """Verify global and fiber responses under positive-mean loading."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.configuration = TowerConfiguration(
            horizontal_force=1.0e6,
            n_elements=4,
            n_gauss=4,
            n_circumferential=8,
            n_radial=2,
        )
        cls.loading = create_pulsating_top_force_history(
            maximum_force=1.0e6,
            force_ratio=0.1,
            period=1.0,
            n_cycles=1,
            increments_per_cycle=20,
        )
        cls.history = build_scaled_history(
            configuration=cls.configuration,
            loading=cls.loading,
        )
        cls.check = verify_direct_checkpoint_solutions(
            configuration=cls.configuration,
            history=cls.history,
        )

    def test_response_histories_match_loading_shape(self) -> None:
        """Every response history must use the loading time grid."""
        expected_shape = self.loading.times.shape

        self.assertEqual(
            self.history.top_displacements.shape,
            expected_shape,
        )
        self.assertEqual(
            self.history.top_rotations.shape,
            expected_shape,
        )
        self.assertEqual(
            self.history.critical_curvatures.shape,
            expected_shape,
        )
        self.assertEqual(
            self.history.critical_minimum_stresses.shape,
            expected_shape,
        )
        self.assertEqual(
            self.history.critical_maximum_stresses.shape,
            expected_shape,
        )

        self.assertTrue(
            np.all(np.isfinite(self.history.top_displacements))
        )
        self.assertTrue(
            np.all(np.isfinite(self.history.top_rotations))
        )
        self.assertTrue(
            np.all(np.isfinite(self.history.critical_curvatures))
        )
        self.assertTrue(
            np.all(
                np.isfinite(
                    self.history.critical_minimum_stresses
                )
            )
        )
        self.assertTrue(
            np.all(
                np.isfinite(
                    self.history.critical_maximum_stresses
                )
            )
        )

    def test_minimum_to_maximum_response_ratio_is_point_one(
        self,
    ) -> None:
        """Elastic displacement, curvature, and stress must follow R_F."""
        quarter = self.loading.increments_per_cycle // 4
        maximum_index = quarter
        minimum_index = 3 * quarter
        expected_ratio = self.loading.force_ratio

        response_pairs = (
            self.history.top_displacements,
            self.history.top_rotations,
            self.history.critical_curvatures,
            self.history.critical_minimum_stresses,
            self.history.critical_maximum_stresses,
        )

        for response in response_pairs:
            actual_ratio = (
                response[minimum_index]
                / response[maximum_index]
            )
            self.assertAlmostEqual(
                float(actual_ratio),
                expected_ratio,
                places=12,
            )

    def test_extreme_fiber_signs_do_not_reverse(self) -> None:
        """
        Each section side must retain its stress sign under positive loading.

        The plotted minimum and maximum stresses refer to opposite sides of
        the annular section. Their opposite signs are caused by bending and
        do not represent stress reversal at one material point.
        """
        self.assertTrue(
            np.all(
                self.history.critical_minimum_stresses < 0.0
            )
        )
        self.assertTrue(
            np.all(
                self.history.critical_maximum_stresses > 0.0
            )
        )

    def test_critical_gauss_point_is_inside_tower(self) -> None:
        """The reported discrete critical point must be a valid Gauss point."""
        self.assertGreaterEqual(
            self.history.critical_element_index,
            0,
        )
        self.assertLess(
            self.history.critical_element_index,
            self.configuration.n_elements,
        )
        self.assertGreaterEqual(
            self.history.critical_gauss_index,
            0,
        )
        self.assertLess(
            self.history.critical_gauss_index,
            self.configuration.n_gauss,
        )
        self.assertGreater(
            self.history.critical_height,
            0.0,
        )
        self.assertLess(
            self.history.critical_height,
            self.configuration.height,
        )

    def test_scaled_history_matches_direct_finite_element_solves(
        self,
    ) -> None:
        """Five independent FE solutions must confirm linear scaling."""
        tolerance = 1.0e-9

        self.assertLess(
            self.check.top_displacement_error,
            tolerance,
        )
        self.assertLess(
            self.check.top_rotation_error,
            tolerance,
        )
        self.assertLess(
            self.check.curvature_error,
            tolerance,
        )
        self.assertLess(
            self.check.minimum_stress_error,
            tolerance,
        )
        self.assertLess(
            self.check.maximum_stress_error,
            tolerance,
        )


if __name__ == "__main__":
    unittest.main()

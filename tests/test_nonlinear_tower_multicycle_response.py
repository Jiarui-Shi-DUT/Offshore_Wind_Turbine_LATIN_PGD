# -*- coding: utf-8 -*-
"""Fast tests for the real multi-cycle nonlinear tower example."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_multicycle_response import (
    compare_consecutive_cycles,
    normalized_drift,
    relative_change,
    run_multicycle_tower_analysis,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestSimilarityUtilities(unittest.TestCase):
    """Test dimensionless cycle-comparison utilities."""

    def test_relative_change(self) -> None:
        self.assertAlmostEqual(
            relative_change(
                current=12.0,
                previous=10.0,
            ),
            2.0 / 12.0,
            places=15,
        )
        self.assertEqual(
            relative_change(
                current=0.0,
                previous=0.0,
            ),
            0.0,
        )

    def test_normalized_drift(self) -> None:
        self.assertAlmostEqual(
            normalized_drift(
                current=0.3,
                previous=0.1,
                cycle_scale=2.0,
            ),
            0.1,
            places=15,
        )

    def test_invalid_similarity_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            relative_change(
                current=np.nan,
                previous=1.0,
            )
        with self.assertRaises(ValueError):
            relative_change(
                current=1.0,
                previous=1.0,
                scale_floor=0.0,
            )
        with self.assertRaises(ValueError):
            normalized_drift(
                current=1.0,
                previous=0.0,
                cycle_scale=-1.0,
            )


class TestRealMulticycleTowerResponse(unittest.TestCase):
    """Run one small two-cycle finite-element analysis."""

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

        cls.result = run_multicycle_tower_analysis(
            configuration=configuration,
            material=material,
            force_amplitude=1.0e6,
            period=10.0,
            n_cycles=2,
            increments_per_cycle=8,
            similarity_tolerance=1.0e-3,
            max_iterations=40,
        )

    def test_history_dimensions(self) -> None:
        response = self.result.response
        loading = response.loading

        self.assertEqual(loading.n_cycles, 2)
        self.assertEqual(
            loading.n_time_points,
            17,
        )
        self.assertEqual(
            response.top_displacements.shape,
            (17,),
        )
        self.assertEqual(
            response.iterations.shape,
            (17,),
        )
        self.assertEqual(
            response.fiber_states.shape[0],
            17,
        )

    def test_exact_reversed_loading_checkpoints(self) -> None:
        forces = self.result.response.loading.forces

        expected = np.array(
            [
                0.0,
                1.0e6,
                0.0,
                -1.0e6,
                0.0,
                1.0e6,
                0.0,
                -1.0e6,
                0.0,
            ],
            dtype=np.float64,
        )
        checkpoint_indices = np.array(
            [0, 2, 4, 6, 8, 10, 12, 14, 16],
            dtype=np.int64,
        )

        np.testing.assert_allclose(
            forces[checkpoint_indices],
            expected,
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_cycle_diagnostics_cover_both_cycles(self) -> None:
        diagnostics = self.result.diagnostics

        self.assertEqual(
            diagnostics.n_cycles,
            2,
        )
        np.testing.assert_array_equal(
            diagnostics.cycle_numbers,
            np.array([1, 2], dtype=np.int64),
        )
        self.assertEqual(
            len(self.result.similarities),
            1,
        )

    def test_stress_reverses_in_each_cycle(self) -> None:
        for cycle in self.result.diagnostics.cycles:
            with self.subTest(
                cycle=cycle.cycle_number
            ):
                self.assertLess(
                    cycle.critical_stress_at_positive_peak
                    * cycle.critical_stress_at_negative_peak,
                    0.0,
                )
                self.assertGreater(
                    cycle.critical_stress_range,
                    0.0,
                )

    def test_damage_is_non_decreasing(self) -> None:
        response = self.result.response
        diagnostics = self.result.diagnostics

        self.assertTrue(
            np.all(
                np.diff(response.maximum_damages)
                >= -1.0e-14
            )
        )
        self.assertTrue(
            np.all(
                diagnostics.maximum_damage_increments
                >= -1.0e-14
            )
        )
        self.assertGreaterEqual(
            diagnostics.maximum_damage_ends[-1],
            diagnostics.maximum_damage_ends[0],
        )

    def test_cycle_metrics_are_finite(self) -> None:
        diagnostics = self.result.diagnostics

        arrays = (
            diagnostics.residual_displacements,
            diagnostics.displacement_ranges,
            diagnostics.critical_stress_ranges,
            diagnostics.critical_plastic_strain_ends,
            diagnostics.critical_backstress_ends,
            diagnostics.critical_r_bar_ends,
            diagnostics.maximum_damage_ends,
            diagnostics.external_work_magnitudes,
        )
        for values in arrays:
            self.assertTrue(
                np.all(np.isfinite(values))
            )

        self.assertTrue(
            np.all(
                diagnostics.displacement_ranges > 0.0
            )
        )
        self.assertTrue(
            np.all(
                diagnostics.critical_stress_ranges > 0.0
            )
        )
        self.assertTrue(
            np.all(
                diagnostics.external_work_magnitudes > 0.0
            )
        )

    def test_solver_convergence_metrics(self) -> None:
        response = self.result.response

        self.assertTrue(
            np.all(response.iterations >= 1)
        )
        self.assertLessEqual(
            int(np.max(response.iterations)),
            40,
        )
        self.assertTrue(
            np.all(response.residual_norms >= 0.0)
        )
        self.assertLess(
            float(np.max(response.residual_norms)),
            1.0e-1,
        )

    def test_similarity_metrics_are_consistent(self) -> None:
        similarity = self.result.similarities[0]

        components = np.array(
            [
                similarity.displacement_range_change,
                similarity.stress_range_change,
                similarity.external_work_change,
                similarity.normalized_residual_drift,
                similarity.normalized_plastic_strain_drift,
            ],
            dtype=np.float64,
        )

        self.assertTrue(
            np.all(np.isfinite(components))
        )
        self.assertTrue(
            np.all(components >= 0.0)
        )
        self.assertAlmostEqual(
            similarity.maximum_indicator,
            float(np.max(components)),
            places=15,
        )
        self.assertEqual(
            similarity.within_tolerance,
            similarity.maximum_indicator
            <= similarity.tolerance,
        )

    def test_similarity_flag_changes_with_selected_tolerance(
        self,
    ) -> None:
        previous = self.result.diagnostics.cycles[0]
        current = self.result.diagnostics.cycles[1]

        loose = compare_consecutive_cycles(
            previous=previous,
            current=current,
            tolerance=1.0,
        )
        strict_tolerance = max(
            loose.maximum_indicator * 0.5,
            1.0e-15,
        )
        strict = compare_consecutive_cycles(
            previous=previous,
            current=current,
            tolerance=strict_tolerance,
        )

        self.assertTrue(loose.within_tolerance)
        if loose.maximum_indicator > 0.0:
            self.assertFalse(
                strict.within_tolerance
            )

    def test_final_similarity_property_matches_last_pair(
        self,
    ) -> None:
        self.assertEqual(
            self.result.final_cycle_is_within_tolerance,
            self.result.similarities[-1].within_tolerance,
        )


if __name__ == "__main__":
    unittest.main()

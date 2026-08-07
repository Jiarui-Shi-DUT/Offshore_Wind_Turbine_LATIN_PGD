# -*- coding: utf-8 -*-
"""Fast tests for the tower damage-mechanism separation probe."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_damage_mechanism_probe import (
    compare_cycle_diagnostics,
    relative_difference,
    response_with_common_critical_location,
    run_damage_mechanism_comparison,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestRelativeDifference(unittest.TestCase):
    """Test signed relative differences used in paired comparisons."""

    def test_positive_and_negative_differences(self) -> None:
        self.assertAlmostEqual(
            relative_difference(
                coupled=12.0,
                reference=10.0,
            ),
            0.2,
            places=15,
        )
        self.assertAlmostEqual(
            relative_difference(
                coupled=8.0,
                reference=10.0,
            ),
            -0.2,
            places=15,
        )

    def test_zero_reference_is_protected(self) -> None:
        value = relative_difference(
            coupled=1.0e-15,
            reference=0.0,
            scale_floor=1.0e-14,
        )
        self.assertAlmostEqual(
            value,
            0.1,
            places=15,
        )

    def test_invalid_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            relative_difference(
                coupled=np.nan,
                reference=1.0,
            )
        with self.assertRaises(ValueError):
            relative_difference(
                coupled=1.0,
                reference=1.0,
                scale_floor=0.0,
            )


class TestDamageMechanismComparison(unittest.TestCase):
    """Run one small paired nonlinear finite-element study."""

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

        cls.comparison = run_damage_mechanism_comparison(
            configuration=configuration,
            material=material,
            force_amplitude=1.0e6,
            period=10.0,
            n_cycles=2,
            increments_per_cycle=8,
            similarity_tolerance=1.0e-3,
            max_iterations=40,
        )

    def test_paired_loading_is_identical(self) -> None:
        coupled = self.comparison.coupled.response.loading
        disabled = (
            self.comparison.damage_disabled.response.loading
        )

        self.assertEqual(
            coupled.force_amplitude,
            disabled.force_amplitude,
        )
        self.assertEqual(
            coupled.period,
            disabled.period,
        )
        self.assertEqual(
            coupled.n_cycles,
            disabled.n_cycles,
        )
        self.assertEqual(
            coupled.increments_per_cycle,
            disabled.increments_per_cycle,
        )
        np.testing.assert_allclose(
            coupled.times,
            disabled.times,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            coupled.forces,
            disabled.forces,
            rtol=0.0,
            atol=0.0,
        )

    def test_damage_switch_changes_only_damage_parameter(self) -> None:
        coupled_material = (
            self.comparison.coupled.response.material
        )
        disabled_material = (
            self.comparison.damage_disabled.response.material
        )

        self.assertGreater(
            coupled_material.k_damage,
            0.0,
        )
        self.assertEqual(
            disabled_material.k_damage,
            0.0,
        )

        coupled_values = coupled_material.__dict__.copy()
        disabled_values = disabled_material.__dict__.copy()
        coupled_values.pop("k_damage")
        disabled_values.pop("k_damage")

        self.assertEqual(
            coupled_values,
            disabled_values,
        )

    def test_common_physical_fiber_is_used(self) -> None:
        coupled_response = (
            self.comparison.coupled.response
        )
        disabled_response = (
            self.comparison.damage_disabled.response
        )

        self.assertEqual(
            coupled_response.critical_location,
            disabled_response.critical_location,
        )
        self.assertAlmostEqual(
            coupled_response.critical_height,
            disabled_response.critical_height,
            places=15,
        )
        self.assertAlmostEqual(
            coupled_response.critical_y_coordinate,
            disabled_response.critical_y_coordinate,
            places=15,
        )

    def test_damage_disabled_case_remains_undamaged(self) -> None:
        disabled_response = (
            self.comparison.damage_disabled.response
        )

        np.testing.assert_allclose(
            disabled_response.maximum_damages,
            np.zeros_like(
                disabled_response.maximum_damages
            ),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            self.comparison.disabled_damage_ends,
            np.zeros(
                self.comparison.n_cycles,
                dtype=np.float64,
            ),
            rtol=0.0,
            atol=0.0,
        )

    def test_coupled_case_accumulates_damage(self) -> None:
        coupled_damage_ends = (
            self.comparison.coupled_damage_ends
        )

        self.assertTrue(
            np.all(coupled_damage_ends > 0.0)
        )
        self.assertTrue(
            np.all(
                np.diff(coupled_damage_ends)
                >= -1.0e-14
            )
        )
        self.assertGreater(
            coupled_damage_ends[-1],
            coupled_damage_ends[0],
        )

    def test_cycle_comparisons_are_complete_and_ordered(
        self,
    ) -> None:
        self.assertEqual(
            self.comparison.n_cycles,
            2,
        )
        np.testing.assert_array_equal(
            self.comparison.cycle_numbers,
            np.array([1, 2], dtype=np.int64),
        )

        for expected_number, cycle in enumerate(
            self.comparison.cycles,
            start=1,
        ):
            self.assertEqual(
                cycle.cycle_number,
                expected_number,
            )

    def test_relative_differences_are_consistent(self) -> None:
        for cycle in self.comparison.cycles:
            with self.subTest(
                cycle=cycle.cycle_number
            ):
                expected_displacement = relative_difference(
                    cycle.coupled_displacement_range,
                    cycle.undamaged_displacement_range,
                )
                expected_stress = relative_difference(
                    cycle.coupled_stress_range,
                    cycle.undamaged_stress_range,
                )
                expected_path = relative_difference(
                    cycle.coupled_plastic_strain_path_length,
                    cycle.undamaged_plastic_strain_path_length,
                )
                expected_cumulative_path = relative_difference(
                    cycle.coupled_cumulative_plastic_strain_path_end,
                    cycle.undamaged_cumulative_plastic_strain_path_end,
                )
                expected_work = relative_difference(
                    cycle.coupled_external_work,
                    cycle.undamaged_external_work,
                )

                self.assertAlmostEqual(
                    cycle
                    .displacement_range_relative_difference,
                    expected_displacement,
                    places=15,
                )
                self.assertAlmostEqual(
                    cycle.stress_range_relative_difference,
                    expected_stress,
                    places=15,
                )
                self.assertAlmostEqual(
                    cycle
                    .plastic_strain_path_length_relative_difference,
                    expected_path,
                    places=15,
                )
                self.assertAlmostEqual(
                    cycle
                    .cumulative_plastic_strain_path_relative_difference,
                    expected_cumulative_path,
                    places=15,
                )
                self.assertAlmostEqual(
                    cycle.external_work_relative_difference,
                    expected_work,
                    places=15,
                )

    def test_plastic_path_metrics_match_cycle_diagnostics(
        self,
    ) -> None:
        coupled_cycles = self.comparison.coupled.diagnostics.cycles
        disabled_cycles = (
            self.comparison.damage_disabled.diagnostics.cycles
        )

        for comparison_cycle, coupled, disabled in zip(
            self.comparison.cycles,
            coupled_cycles,
            disabled_cycles,
        ):
            with self.subTest(
                cycle=comparison_cycle.cycle_number
            ):
                self.assertAlmostEqual(
                    comparison_cycle
                    .coupled_plastic_strain_path_length,
                    coupled.critical_plastic_strain_path_length,
                    places=15,
                )
                self.assertAlmostEqual(
                    comparison_cycle
                    .undamaged_plastic_strain_path_length,
                    disabled.critical_plastic_strain_path_length,
                    places=15,
                )
                self.assertAlmostEqual(
                    comparison_cycle
                    .coupled_cumulative_plastic_strain_path_end,
                    coupled
                    .critical_cumulative_plastic_strain_path_at_end,
                    places=15,
                )
                self.assertAlmostEqual(
                    comparison_cycle
                    .undamaged_cumulative_plastic_strain_path_end,
                    disabled
                    .critical_cumulative_plastic_strain_path_at_end,
                    places=15,
                )
                self.assertAlmostEqual(
                    comparison_cycle
                    .coupled_plastic_net_to_path_ratio,
                    coupled.critical_plastic_net_to_path_ratio,
                    places=15,
                )
                self.assertAlmostEqual(
                    comparison_cycle
                    .undamaged_plastic_net_to_path_ratio,
                    disabled.critical_plastic_net_to_path_ratio,
                    places=15,
                )
                self.assertGreaterEqual(
                    comparison_cycle
                    .coupled_plastic_strain_path_length,
                    0.0,
                )
                self.assertGreaterEqual(
                    comparison_cycle
                    .undamaged_plastic_strain_path_length,
                    0.0,
                )
                self.assertGreaterEqual(
                    comparison_cycle
                    .coupled_plastic_net_to_path_ratio,
                    0.0,
                )
                self.assertLessEqual(
                    comparison_cycle
                    .coupled_plastic_net_to_path_ratio,
                    1.0 + 1.0e-12,
                )
                self.assertGreaterEqual(
                    comparison_cycle
                    .undamaged_plastic_net_to_path_ratio,
                    0.0,
                )
                self.assertLessEqual(
                    comparison_cycle
                    .undamaged_plastic_net_to_path_ratio,
                    1.0 + 1.0e-12,
                )

    def test_response_quantities_are_finite(self) -> None:
        arrays = (
            self.comparison
            .displacement_range_relative_differences,
            self.comparison
            .stress_range_relative_differences,
            self.comparison
            .plastic_strain_path_length_relative_differences,
            self.comparison
            .cumulative_plastic_strain_path_relative_differences,
            self.comparison
            .plastic_net_to_path_ratio_differences,
            self.comparison
            .external_work_relative_differences,
            self.comparison.coupled_damage_ends,
            self.comparison.disabled_damage_ends,
        )
        for values in arrays:
            self.assertTrue(
                np.all(np.isfinite(values))
            )

        for cycle in self.comparison.cycles:
            scalar_values = (
                cycle.coupled_residual_displacement,
                cycle.undamaged_residual_displacement,
                cycle.residual_displacement_difference,
                cycle.coupled_plastic_strain_end,
                cycle.undamaged_plastic_strain_end,
                cycle.plastic_strain_end_difference,
                cycle.coupled_plastic_strain_path_length,
                cycle.undamaged_plastic_strain_path_length,
                cycle
                .plastic_strain_path_length_relative_difference,
                cycle
                .coupled_cumulative_plastic_strain_path_end,
                cycle
                .undamaged_cumulative_plastic_strain_path_end,
                cycle
                .cumulative_plastic_strain_path_relative_difference,
                cycle.coupled_plastic_net_to_path_ratio,
                cycle.undamaged_plastic_net_to_path_ratio,
                cycle.plastic_net_to_path_ratio_difference,
            )
            self.assertTrue(
                np.all(
                    np.isfinite(
                        np.asarray(
                            scalar_values,
                            dtype=np.float64,
                        )
                    )
                )
            )

    def test_both_solvers_converge(self) -> None:
        responses = (
            self.comparison.coupled.response,
            self.comparison.damage_disabled.response,
        )

        for response in responses:
            with self.subTest(
                k_damage=response.material.k_damage
            ):
                self.assertTrue(
                    np.all(response.iterations >= 1)
                )
                self.assertLessEqual(
                    int(np.max(response.iterations)),
                    40,
                )
                self.assertTrue(
                    np.all(
                        response.residual_norms >= 0.0
                    )
                )
                self.assertLess(
                    float(
                        np.max(
                            response.residual_norms
                        )
                    ),
                    1.0e-1,
                )

    def test_reanchoring_preserves_global_histories(self) -> None:
        disabled_response = (
            self.comparison.damage_disabled.response
        )
        coupled_response = (
            self.comparison.coupled.response
        )

        reanchored = response_with_common_critical_location(
            response=disabled_response,
            reference=coupled_response,
        )

        np.testing.assert_allclose(
            reanchored.top_displacements,
            disabled_response.top_displacements,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            reanchored.fiber_states,
            disabled_response.fiber_states,
            rtol=0.0,
            atol=0.0,
        )
        self.assertEqual(
            reanchored.critical_location,
            coupled_response.critical_location,
        )

    def test_cycle_comparison_rebuild_is_reproducible(
        self,
    ) -> None:
        rebuilt = compare_cycle_diagnostics(
            coupled=self.comparison.coupled.diagnostics,
            damage_disabled=(
                self.comparison
                .damage_disabled
                .diagnostics
            ),
        )

        self.assertEqual(
            len(rebuilt),
            len(self.comparison.cycles),
        )
        for actual, expected in zip(
            rebuilt,
            self.comparison.cycles,
        ):
            self.assertEqual(
                actual,
                expected,
            )


if __name__ == "__main__":
    unittest.main()

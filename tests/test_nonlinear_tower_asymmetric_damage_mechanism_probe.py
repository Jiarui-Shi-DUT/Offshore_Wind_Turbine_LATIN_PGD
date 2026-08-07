# -*- coding: utf-8 -*-
"""Tests for asymmetric damage-mechanism separation."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_damage_mechanism_probe import (
    run_asymmetric_damage_mechanism_comparison,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestAsymmetricDamageMechanismComparison(unittest.TestCase):
    """Run a reduced two-cycle paired asymmetric analysis."""

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

        cls.comparison = (
            run_asymmetric_damage_mechanism_comparison(
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
        )

    def test_paired_loading_is_identical(self) -> None:
        comparison = self.comparison
        coupled_loading = (
            comparison.mechanism.coupled.response.loading
        )
        disabled_loading = (
            comparison.mechanism
            .damage_disabled.response.loading
        )

        self.assertAlmostEqual(
            coupled_loading.force_ratio,
            -0.5,
            places=14,
        )
        self.assertEqual(
            coupled_loading.maximum_force,
            disabled_loading.maximum_force,
        )
        self.assertEqual(
            coupled_loading.force_ratio,
            disabled_loading.force_ratio,
        )
        np.testing.assert_allclose(
            coupled_loading.forces,
            disabled_loading.forces,
            rtol=0.0,
            atol=0.0,
        )

    def test_damage_switch_is_isolated(self) -> None:
        comparison = self.comparison
        coupled_material = (
            comparison.mechanism.coupled.response.material
        )
        disabled_material = (
            comparison.mechanism
            .damage_disabled.response.material
        )

        self.assertGreater(coupled_material.k_damage, 0.0)
        self.assertEqual(disabled_material.k_damage, 0.0)

    def test_common_critical_fiber_is_used(self) -> None:
        comparison = self.comparison
        coupled = comparison.mechanism.coupled.response
        disabled = (
            comparison.mechanism
            .damage_disabled.response
        )

        self.assertEqual(
            coupled.critical_location,
            disabled.critical_location,
        )
        self.assertEqual(
            coupled.critical_height,
            disabled.critical_height,
        )
        self.assertEqual(
            coupled.critical_y_coordinate,
            disabled.critical_y_coordinate,
        )

    def test_ratcheting_histories_have_one_value_per_cycle(
        self,
    ) -> None:
        comparison = self.comparison

        histories = (
            comparison.coupled_cycle_displacement_drifts,
            comparison.damage_disabled_cycle_displacement_drifts,
            comparison.cycle_displacement_drift_differences,
            comparison.coupled_plastic_strain_drifts,
            comparison.damage_disabled_plastic_strain_drifts,
            comparison.plastic_strain_drift_differences,
        )
        for values in histories:
            self.assertEqual(values.shape, (2,))
            self.assertTrue(np.all(np.isfinite(values)))

    def test_drift_differences_are_exact_pairwise_differences(
        self,
    ) -> None:
        comparison = self.comparison

        np.testing.assert_allclose(
            comparison.cycle_displacement_drift_differences,
            (
                comparison.coupled_cycle_displacement_drifts
                - comparison
                .damage_disabled_cycle_displacement_drifts
            ),
            rtol=1.0e-12,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            comparison.plastic_strain_drift_differences,
            (
                comparison.coupled_plastic_strain_drifts
                - comparison
                .damage_disabled_plastic_strain_drifts
            ),
            rtol=1.0e-12,
            atol=1.0e-15,
        )

    def test_damage_disabled_case_remains_zero_damage(self) -> None:
        disabled = (
            self.comparison.mechanism.damage_disabled
        )

        np.testing.assert_allclose(
            disabled.diagnostics.maximum_damage_ends,
            np.zeros(2, dtype=np.float64),
            rtol=0.0,
            atol=1.0e-15,
        )

    def test_coupled_damage_is_non_decreasing(self) -> None:
        coupled = self.comparison.mechanism.coupled

        self.assertTrue(
            np.all(
                np.diff(
                    coupled.diagnostics.maximum_damage_ends
                )
                >= -1.0e-14
            )
        )

    def test_base_mechanism_comparison_covers_all_cycles(
        self,
    ) -> None:
        comparison = self.comparison

        self.assertEqual(comparison.n_cycles, 2)
        np.testing.assert_array_equal(
            comparison.cycle_numbers,
            np.array([1, 2], dtype=np.int64),
        )
        self.assertEqual(
            comparison.mechanism.n_cycles,
            2,
        )

    def test_both_nonlinear_solves_converge(self) -> None:
        comparison = self.comparison

        for result in (
            comparison.mechanism.coupled,
            comparison.mechanism.damage_disabled,
        ):
            response = result.response
            self.assertTrue(
                np.all(response.iterations >= 1)
            )
            self.assertLessEqual(
                int(np.max(response.iterations)),
                40,
            )
            self.assertLess(
                float(np.max(response.residual_norms)),
                1.0e-1,
            )


if __name__ == "__main__":
    unittest.main()

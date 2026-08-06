# -*- coding: utf-8 -*-
"""Unit tests for per-cycle nonlinear tower diagnostics."""

from __future__ import annotations

import unittest

import numpy as np

from examples.nonlinear_tower_multicycle_diagnostics import (
    cycle_indices,
    extract_cycle_diagnostics,
    extract_multicycle_diagnostics,
    plastic_strain_path_length,
    signed_force_displacement_work,
)
from examples.nonlinear_tower_reversed_response import (
    NonlinearReversedResponse,
)
from fem.tower_loading import create_reversed_top_force_history
from material.viscoplastic_damage_1d import MaterialParameters


def create_synthetic_multicycle_response(
) -> NonlinearReversedResponse:
    """Create a deterministic two-cycle response without FE solving."""
    loading = create_reversed_top_force_history(
        force_amplitude=100.0,
        period=10.0,
        n_cycles=2,
        increments_per_cycle=4,
    )
    material = MaterialParameters()

    top_displacements = np.array(
        [
            0.0,
            1.0,
            0.2,
            -0.8,
            0.1,
            0.9,
            0.15,
            -0.7,
            0.2,
        ],
        dtype=np.float64,
    )
    top_rotations = 0.01 * top_displacements
    base_horizontal_reactions = -loading.forces
    base_moment_reactions = (
        -10.0 * loading.forces
    )

    iterations = np.array(
        [1, 2, 3, 2, 1, 2, 4, 3, 2],
        dtype=np.int64,
    )
    residual_norms = np.array(
        [
            0.0,
            1.0e-3,
            2.0e-3,
            4.0e-3,
            1.0e-3,
            2.0e-3,
            5.0e-3,
            3.0e-3,
            1.0e-3,
        ],
        dtype=np.float64,
    )

    fiber_strains = np.array(
        [
            0.0,
            -8.0e-4,
            -1.0e-4,
            7.0e-4,
            1.0e-4,
            -7.0e-4,
            -5.0e-5,
            6.5e-4,
            1.5e-4,
        ],
        dtype=np.float64,
    ).reshape((-1, 1, 1, 1))

    fiber_stresses = np.array(
        [
            0.0,
            -80.0,
            5.0,
            82.0,
            -3.0,
            -78.0,
            4.0,
            79.0,
            -2.0,
        ],
        dtype=np.float64,
    ).reshape((-1, 1, 1, 1))

    plastic_strains = np.array(
        [
            0.0,
            -1.0e-2,
            -1.5e-2,
            5.0e-3,
            8.0e-3,
            -5.0e-3,
            -1.0e-2,
            9.0e-3,
            1.2e-2,
        ],
        dtype=np.float64,
    )
    alphas = np.array(
        [
            0.0,
            -1.0e-3,
            -2.0e-3,
            5.0e-4,
            1.0e-3,
            -8.0e-4,
            -1.5e-3,
            7.0e-4,
            1.2e-3,
        ],
        dtype=np.float64,
    )
    r_bars = np.array(
        [
            0.0,
            1.0e-3,
            2.0e-3,
            3.0e-3,
            4.0e-3,
            5.0e-3,
            6.0e-3,
            7.0e-3,
            8.0e-3,
        ],
        dtype=np.float64,
    )
    damages = np.array(
        [
            0.00,
            0.01,
            0.02,
            0.03,
            0.04,
            0.05,
            0.06,
            0.07,
            0.08,
        ],
        dtype=np.float64,
    )

    fiber_states = np.column_stack(
        (
            plastic_strains,
            alphas,
            r_bars,
            damages,
        )
    ).reshape((-1, 1, 1, 1, 4))

    return NonlinearReversedResponse(
        loading=loading,
        material=material,
        analysis_times=(
            loading.time_increment + loading.times
        ),
        top_displacements=top_displacements,
        top_rotations=top_rotations,
        base_horizontal_reactions=(
            base_horizontal_reactions
        ),
        base_moment_reactions=base_moment_reactions,
        iterations=iterations,
        residual_norms=residual_norms,
        fiber_strains=fiber_strains,
        fiber_stresses=fiber_stresses,
        fiber_states=fiber_states,
        critical_location=(0, 0, 0),
        critical_height=1.0,
        critical_y_coordinate=-2.0,
    )


class TestCycleIndices(unittest.TestCase):
    """Test exact indexing of complete reversed cycles."""

    def setUp(self) -> None:
        self.loading = create_reversed_top_force_history(
            force_amplitude=100.0,
            period=10.0,
            n_cycles=2,
            increments_per_cycle=4,
        )

    def test_first_cycle_indices(self) -> None:
        indices = cycle_indices(
            loading=self.loading,
            cycle_number=1,
        )

        self.assertEqual(indices.start, 0)
        self.assertEqual(indices.positive_peak, 1)
        self.assertEqual(indices.first_zero, 2)
        self.assertEqual(indices.negative_peak, 3)
        self.assertEqual(indices.end, 4)
        self.assertEqual(indices.slice, slice(0, 5))

    def test_second_cycle_starts_at_previous_cycle_end(self) -> None:
        first = cycle_indices(
            loading=self.loading,
            cycle_number=1,
        )
        second = cycle_indices(
            loading=self.loading,
            cycle_number=2,
        )

        self.assertEqual(second.start, first.end)
        self.assertEqual(second.positive_peak, 5)
        self.assertEqual(second.first_zero, 6)
        self.assertEqual(second.negative_peak, 7)
        self.assertEqual(second.end, 8)

    def test_invalid_cycle_numbers_are_rejected(self) -> None:
        invalid_values = (0, 3, -1, 1.5, True)

        for value in invalid_values:
            with self.subTest(value=value):
                with self.assertRaises(
                    (TypeError, ValueError)
                ):
                    cycle_indices(
                        loading=self.loading,
                        cycle_number=value,
                    )


class TestSignedForceDisplacementWork(unittest.TestCase):
    """Test the discrete force-displacement line integral."""

    def test_known_path_work(self) -> None:
        forces = np.array(
            [0.0, 100.0, 0.0, -100.0, 0.0],
            dtype=np.float64,
        )
        displacements = np.array(
            [0.0, 1.0, 0.2, -0.8, 0.1],
            dtype=np.float64,
        )

        work = signed_force_displacement_work(
            forces=forces,
            displacements=displacements,
        )

        self.assertAlmostEqual(work, 15.0, places=12)

    def test_zero_force_path_has_zero_work(self) -> None:
        work = signed_force_displacement_work(
            forces=np.zeros(4, dtype=np.float64),
            displacements=np.array(
                [0.0, 1.0, -1.0, 0.5],
                dtype=np.float64,
            ),
        )

        self.assertEqual(work, 0.0)

    def test_invalid_paths_are_rejected(self) -> None:
        invalid_cases = (
            (
                np.zeros((2, 2), dtype=np.float64),
                np.zeros(4, dtype=np.float64),
            ),
            (
                np.zeros(3, dtype=np.float64),
                np.zeros(4, dtype=np.float64),
            ),
            (
                np.zeros(1, dtype=np.float64),
                np.zeros(1, dtype=np.float64),
            ),
            (
                np.array([0.0, np.nan]),
                np.zeros(2, dtype=np.float64),
            ),
        )

        for forces, displacements in invalid_cases:
            with self.subTest(
                force_shape=forces.shape,
                displacement_shape=displacements.shape,
            ):
                with self.assertRaises(ValueError):
                    signed_force_displacement_work(
                        forces=forces,
                        displacements=displacements,
                    )



class TestPlasticStrainPathLength(unittest.TestCase):
    """Test accumulated absolute plastic-strain path length."""

    def test_reversing_path_exceeds_net_increment(self) -> None:
        plastic_strains = np.array(
            [0.0, -0.01, -0.015, 0.005, 0.008],
            dtype=np.float64,
        )

        path_length = plastic_strain_path_length(
            plastic_strains
        )
        net_increment = float(
            plastic_strains[-1] - plastic_strains[0]
        )

        self.assertAlmostEqual(
            path_length,
            0.038,
            places=15,
        )
        self.assertGreater(
            path_length,
            abs(net_increment),
        )

    def test_monotonic_path_equals_absolute_net_increment(
        self,
    ) -> None:
        plastic_strains = np.array(
            [0.0, 0.01, 0.03, 0.05],
            dtype=np.float64,
        )

        path_length = plastic_strain_path_length(
            plastic_strains
        )

        self.assertAlmostEqual(
            path_length,
            0.05,
            places=15,
        )
        self.assertAlmostEqual(
            path_length,
            abs(
                float(
                    plastic_strains[-1]
                    - plastic_strains[0]
                )
            ),
            places=15,
        )

    def test_invalid_plastic_histories_are_rejected(
        self,
    ) -> None:
        invalid_histories = (
            np.zeros((2, 2), dtype=np.float64),
            np.zeros(1, dtype=np.float64),
            np.array([0.0, np.nan], dtype=np.float64),
        )

        for values in invalid_histories:
            with self.subTest(shape=values.shape):
                with self.assertRaises(ValueError):
                    plastic_strain_path_length(values)


class TestCycleDiagnosticsExtraction(unittest.TestCase):
    """Test scalar extraction from a synthetic multi-cycle history."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.response = create_synthetic_multicycle_response()

    def test_first_cycle_global_diagnostics(self) -> None:
        cycle = extract_cycle_diagnostics(
            response=self.response,
            cycle_number=1,
        )

        self.assertEqual(cycle.cycle_number, 1)
        self.assertAlmostEqual(
            cycle.displacement_at_positive_peak,
            1.0,
        )
        self.assertAlmostEqual(
            cycle.displacement_at_negative_peak,
            -0.8,
        )
        self.assertAlmostEqual(
            cycle.maximum_displacement,
            1.0,
        )
        self.assertAlmostEqual(
            cycle.minimum_displacement,
            -0.8,
        )
        self.assertAlmostEqual(
            cycle.residual_displacement,
            0.1,
        )
        self.assertAlmostEqual(
            cycle.displacement_range,
            1.8,
        )
        self.assertAlmostEqual(
            cycle.signed_external_work,
            15.0,
        )
        self.assertAlmostEqual(
            cycle.external_work_magnitude,
            15.0,
        )
        self.assertEqual(
            cycle.maximum_newton_iterations,
            3,
        )
        self.assertAlmostEqual(
            cycle.maximum_residual_norm,
            4.0e-3,
        )

    def test_first_cycle_fixed_fiber_diagnostics(self) -> None:
        cycle = extract_cycle_diagnostics(
            response=self.response,
            cycle_number=1,
        )

        self.assertAlmostEqual(
            cycle.critical_stress_at_positive_peak,
            -80.0,
        )
        self.assertAlmostEqual(
            cycle.critical_stress_at_negative_peak,
            82.0,
        )
        self.assertAlmostEqual(
            cycle.maximum_critical_stress,
            82.0,
        )
        self.assertAlmostEqual(
            cycle.minimum_critical_stress,
            -80.0,
        )
        self.assertAlmostEqual(
            cycle.critical_stress_range,
            162.0,
        )

        self.assertAlmostEqual(
            cycle.critical_plastic_strain_at_end,
            8.0e-3,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_strain_increment,
            8.0e-3,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_strain_range,
            2.3e-2,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_strain_path_length,
            3.8e-2,
        )
        self.assertAlmostEqual(
            cycle
            .critical_cumulative_plastic_strain_path_at_end,
            3.8e-2,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_net_to_path_ratio,
            8.0 / 38.0,
            places=15,
        )

        self.assertAlmostEqual(
            cycle.critical_backstress_at_end,
            5.5,
        )
        self.assertAlmostEqual(
            cycle.minimum_critical_backstress,
            -11.0,
        )
        self.assertAlmostEqual(
            cycle.maximum_critical_backstress,
            5.5,
        )
        self.assertAlmostEqual(
            cycle.critical_backstress_range,
            16.5,
        )

        self.assertAlmostEqual(
            cycle.critical_r_bar_at_end,
            4.0e-3,
        )
        self.assertAlmostEqual(
            cycle.critical_r_bar_increment,
            4.0e-3,
        )

    def test_first_cycle_damage_diagnostics(self) -> None:
        cycle = extract_cycle_diagnostics(
            response=self.response,
            cycle_number=1,
        )

        self.assertAlmostEqual(
            cycle.maximum_damage_at_start,
            0.0,
        )
        self.assertAlmostEqual(
            cycle.maximum_damage_at_end,
            0.04,
        )
        self.assertAlmostEqual(
            cycle.maximum_damage_increment,
            0.04,
        )
        self.assertAlmostEqual(
            cycle.critical_damage_at_start,
            0.0,
        )
        self.assertAlmostEqual(
            cycle.critical_damage_at_end,
            0.04,
        )
        self.assertAlmostEqual(
            cycle.critical_damage_increment,
            0.04,
        )

    def test_second_cycle_uses_shared_boundary_state(self) -> None:
        cycle = extract_cycle_diagnostics(
            response=self.response,
            cycle_number=2,
        )

        self.assertEqual(cycle.indices.start, 4)
        self.assertAlmostEqual(
            cycle.maximum_damage_at_start,
            0.04,
        )
        self.assertAlmostEqual(
            cycle.maximum_damage_at_end,
            0.08,
        )
        self.assertAlmostEqual(
            cycle.maximum_damage_increment,
            0.04,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_strain_increment,
            4.0e-3,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_strain_path_length,
            4.0e-2,
        )
        self.assertAlmostEqual(
            cycle
            .critical_cumulative_plastic_strain_path_at_end,
            7.8e-2,
        )
        self.assertAlmostEqual(
            cycle.critical_plastic_net_to_path_ratio,
            0.1,
            places=15,
        )
        self.assertAlmostEqual(
            cycle.residual_displacement,
            0.2,
        )
        self.assertAlmostEqual(
            cycle.signed_external_work,
            0.0,
            places=12,
        )


class TestMulticycleDiagnostics(unittest.TestCase):
    """Test ordered array views across all cycles."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.response = create_synthetic_multicycle_response()
        cls.diagnostics = extract_multicycle_diagnostics(
            cls.response
        )

    def test_cycle_count_and_numbering(self) -> None:
        self.assertEqual(self.diagnostics.n_cycles, 2)
        np.testing.assert_array_equal(
            self.diagnostics.cycle_numbers,
            np.array([1, 2], dtype=np.int64),
        )

    def test_array_histories_follow_cycle_order(self) -> None:
        np.testing.assert_allclose(
            self.diagnostics.residual_displacements,
            np.array([0.1, 0.2]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            self.diagnostics.maximum_damage_ends,
            np.array([0.04, 0.08]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            self.diagnostics.maximum_damage_increments,
            np.array([0.04, 0.04]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            self.diagnostics.external_work_magnitudes,
            np.array([15.0, 0.0]),
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            self.diagnostics
            .critical_plastic_strain_path_lengths,
            np.array([0.038, 0.040]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            self.diagnostics
            .critical_cumulative_plastic_strain_path_ends,
            np.array([0.038, 0.078]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            self.diagnostics
            .critical_plastic_net_to_path_ratios,
            np.array([8.0 / 38.0, 0.1]),
            rtol=0.0,
            atol=1.0e-15,
        )

    def test_extracted_arrays_are_defensive(self) -> None:
        values = self.diagnostics.residual_displacements
        values[0] = 999.0

        self.assertAlmostEqual(
            self.diagnostics.residual_displacements[0],
            0.1,
        )


if __name__ == "__main__":
    unittest.main()

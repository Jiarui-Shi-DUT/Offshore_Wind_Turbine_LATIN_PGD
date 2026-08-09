# -*- coding: utf-8 -*-
"""Regression tests for tower cycle-phase SVD diagnostics."""

from __future__ import annotations

import unittest

import numpy as np

from examples.nonlinear_tower_low_rank_diagnostics import (
    analyze_field_low_rank,
    analyze_tensor_low_rank,
    analyze_tower_low_rank,
    compute_svd_spectrum,
    cycle_increment_field,
    flatten_cycle_phase_field,
    mode_unfolding,
)
from examples.nonlinear_tower_snapshot_tensor import (
    TowerCyclePhaseSnapshots,
)


class TestTowerLowRankHelpers(unittest.TestCase):
    """Verify tensor preparation and mode unfolding."""

    def test_flatten_preserves_cycle_phase_and_flattens_space(self) -> None:
        values = np.arange(
            2 * 3 * 2 * 4,
            dtype=np.float64,
        ).reshape((2, 3, 2, 4))

        tensor = flatten_cycle_phase_field(values)

        self.assertEqual(tensor.shape, (2, 3, 8))
        np.testing.assert_array_equal(
            tensor[1, 2],
            values[1, 2].reshape(-1),
        )

    def test_cycle_increment_uses_each_cycle_start(self) -> None:
        values = np.array(
            [
                [[10.0, 20.0], [12.0, 19.0], [13.0, 25.0]],
                [[30.0, 40.0], [31.0, 44.0], [35.0, 41.0]],
            ],
            dtype=np.float64,
        )

        increment = cycle_increment_field(values)

        np.testing.assert_allclose(
            increment[:, 0],
            np.zeros((2, 2), dtype=np.float64),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            increment[0, -1],
            np.array([3.0, 5.0]),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            increment[1, -1],
            np.array([5.0, 1.0]),
            rtol=0.0,
            atol=0.0,
        )

    def test_mode_unfolding_shapes_are_exact(self) -> None:
        tensor = np.arange(
            3 * 4 * 5,
            dtype=np.float64,
        ).reshape((3, 4, 5))

        self.assertEqual(
            mode_unfolding(tensor, "cycle").shape,
            (3, 20),
        )
        self.assertEqual(
            mode_unfolding(tensor, "phase").shape,
            (4, 15),
        )
        self.assertEqual(
            mode_unfolding(tensor, "space").shape,
            (5, 12),
        )

        np.testing.assert_array_equal(
            mode_unfolding(tensor, "cycle")[2],
            tensor[2].reshape(-1),
        )
        np.testing.assert_array_equal(
            mode_unfolding(tensor, "phase")[1],
            tensor[:, 1, :].reshape(-1),
        )
        np.testing.assert_array_equal(
            mode_unfolding(tensor, "space")[3],
            tensor[:, :, 3].reshape(-1),
        )

    def test_invalid_mode_is_rejected(self) -> None:
        tensor = np.zeros((2, 3, 4), dtype=np.float64)

        with self.assertRaises(ValueError):
            mode_unfolding(tensor, "time")

    def test_nonfinite_field_is_rejected(self) -> None:
        values = np.zeros((2, 3, 4), dtype=np.float64)
        values[0, 0, 0] = np.nan

        with self.assertRaises(ValueError):
            flatten_cycle_phase_field(values)


class TestSvdSpectrum(unittest.TestCase):
    """Verify singular-value energy and rank diagnostics."""

    def test_rank_one_tensor_is_rank_one_in_all_modes(self) -> None:
        cycle = np.array([1.0, 2.0, 4.0], dtype=np.float64)
        phase = np.array(
            [1.0, -1.0, 2.0, 0.5],
            dtype=np.float64,
        )
        space = np.array(
            [2.0, 3.0, -1.0, 4.0, 5.0],
            dtype=np.float64,
        )

        tensor = (
            cycle[:, np.newaxis, np.newaxis]
            * phase[np.newaxis, :, np.newaxis]
            * space[np.newaxis, np.newaxis, :]
        )

        diagnostics = analyze_tensor_low_rank(tensor)

        for mode in ("cycle", "phase", "space"):
            spectrum = diagnostics.mode(mode)
            self.assertEqual(
                spectrum.rank_for_energy(0.999999999999),
                1,
            )
            self.assertAlmostEqual(
                spectrum.error_after_rank(1),
                0.0,
                places=12,
            )

    def test_known_singular_values_give_correct_energy_rank(self) -> None:
        matrix = np.diag(
            np.array([4.0, 3.0], dtype=np.float64)
        )
        spectrum = compute_svd_spectrum(matrix)

        np.testing.assert_allclose(
            spectrum.singular_values,
            np.array([4.0, 3.0]),
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            spectrum.energy_fractions,
            np.array([16.0 / 25.0, 9.0 / 25.0]),
            rtol=0.0,
            atol=1.0e-14,
        )
        self.assertEqual(spectrum.rank_for_energy(0.60), 1)
        self.assertEqual(spectrum.rank_for_energy(0.90), 2)
        self.assertAlmostEqual(
            spectrum.error_after_rank(1),
            3.0 / 5.0,
            places=14,
        )
        self.assertAlmostEqual(
            spectrum.error_after_rank(2),
            0.0,
            places=14,
        )

    def test_zero_matrix_has_exact_rank_zero_reconstruction(self) -> None:
        spectrum = compute_svd_spectrum(
            np.zeros((3, 5), dtype=np.float64)
        )

        self.assertTrue(spectrum.is_zero)
        self.assertEqual(spectrum.rank_for_energy(0.999), 0)
        self.assertEqual(spectrum.error_after_rank(0), 0.0)
        self.assertEqual(spectrum.error_after_rank(1), 0.0)

    def test_rank_energy_validation(self) -> None:
        spectrum = compute_svd_spectrum(
            np.eye(2, dtype=np.float64)
        )

        with self.assertRaises(ValueError):
            spectrum.rank_for_energy(0.0)
        with self.assertRaises(ValueError):
            spectrum.rank_for_energy(1.01)
        with self.assertRaises(ValueError):
            spectrum.error_after_rank(-1)
        with self.assertRaises(ValueError):
            spectrum.error_after_rank(3)


class TestTowerLowRankDiagnostics(unittest.TestCase):
    """Verify target-field extraction without running the nonlinear solver."""

    @staticmethod
    def create_snapshots() -> TowerCyclePhaseSnapshots:
        n_cycles = 3
        n_phase = 4
        n_dof = 6
        n_elements = 2
        n_gauss = 1
        n_fibers = 2

        cycle_numbers = np.arange(
            1,
            n_cycles + 1,
            dtype=np.int64,
        )
        phase_times = np.linspace(
            0.0,
            10.0,
            n_phase,
            dtype=np.float64,
        )
        phase_fractions = phase_times / 10.0
        phase_forces = np.array(
            [0.25, 1.0, -0.5, 0.25],
            dtype=np.float64,
        )
        analysis_times = (
            np.arange(n_cycles, dtype=np.float64)[:, np.newaxis]
            * 10.0
            + phase_times[np.newaxis, :]
        )

        cycle_factor = np.array(
            [1.0, 1.2, 1.5],
            dtype=np.float64,
        )[:, np.newaxis, np.newaxis]
        phase_factor = np.array(
            [0.5, 1.0, -0.4, 0.55],
            dtype=np.float64,
        )[np.newaxis, :, np.newaxis]

        displacement_space = np.linspace(
            0.0,
            1.0,
            n_dof,
            dtype=np.float64,
        )[np.newaxis, np.newaxis, :]
        nodal_displacements = (
            cycle_factor
            * phase_factor
            * displacement_space
        )

        fiber_space = np.array(
            [1.0, 2.0, 3.0, 4.0],
            dtype=np.float64,
        ).reshape(
            (1, 1, n_elements, n_gauss, n_fibers)
        )
        fiber_base = (
            cycle_factor[..., np.newaxis, np.newaxis]
            * phase_factor[..., np.newaxis, np.newaxis]
            * fiber_space
        )
        fiber_strains = 1.0e-3 * fiber_base
        fiber_stresses = 100.0 * fiber_base

        fiber_states = np.zeros(
            fiber_strains.shape + (4,),
            dtype=np.float64,
        )
        fiber_states[..., 0] = 2.0e-4 * fiber_base
        fiber_states[..., 1] = 1.0e-4 * fiber_base
        fiber_states[..., 2] = 1.0e-3 * np.abs(fiber_base)
        fiber_states[..., 3] = (
            1.0e-2 * np.abs(fiber_base)
        )

        return TowerCyclePhaseSnapshots(
            cycle_numbers=cycle_numbers,
            phase_times=phase_times,
            phase_fractions=phase_fractions,
            phase_forces=phase_forces,
            analysis_times=analysis_times,
            nodal_displacements=nodal_displacements,
            fiber_strains=fiber_strains,
            fiber_stresses=fiber_stresses,
            fiber_states=fiber_states,
        )

    def test_target_fields_have_expected_flattened_shapes(self) -> None:
        snapshots = self.create_snapshots()
        diagnostics = analyze_tower_low_rank(snapshots)

        self.assertEqual(
            diagnostics.nodal_displacements.raw.tensor_shape,
            (3, 4, 6),
        )
        for field in (
            diagnostics.fiber_stresses,
            diagnostics.fiber_plastic_strains,
            diagnostics.fiber_damages,
        ):
            self.assertEqual(
                field.raw.tensor_shape,
                (3, 4, 4),
            )
            self.assertEqual(
                field.cycle_increment.tensor_shape,
                (3, 4, 4),
            )

    def test_target_field_dictionary_uses_physical_names(self) -> None:
        diagnostics = analyze_tower_low_rank(
            self.create_snapshots()
        )

        self.assertEqual(
            set(diagnostics.as_dict().keys()),
            {"u", "sigma", "eps_p", "D"},
        )

    def test_cycle_increment_diagnostics_remove_start_baseline(self) -> None:
        snapshots = self.create_snapshots()

        direct = analyze_field_low_rank(
            snapshots.nodal_displacements
        )
        tower = analyze_tower_low_rank(snapshots)

        self.assertEqual(
            direct.cycle_increment.tensor_shape,
            tower.nodal_displacements.cycle_increment.tensor_shape,
        )
        self.assertEqual(
            tower.nodal_displacements.cycle_increment
            .cycle_mode.rank_for_energy(0.999999),
            1,
        )

    def test_wrong_snapshot_type_is_rejected(self) -> None:
        with self.assertRaises(TypeError):
            analyze_tower_low_rank(
                np.zeros((2, 3, 4), dtype=np.float64)
            )


if __name__ == "__main__":
    unittest.main()

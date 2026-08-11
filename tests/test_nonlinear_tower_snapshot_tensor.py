# -*- coding: utf-8 -*-
"""Regression tests for slow-cycle x fast-phase tower snapshots."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_response import (
    run_nonlinear_asymmetric_analysis,
)
from examples.nonlinear_tower_snapshot_tensor import (
    TowerCyclePhaseSnapshots,
    build_tower_cycle_phase_snapshots,
    load_tower_cycle_phase_snapshots,
    save_tower_cycle_phase_snapshots,
    select_tower_cycle_range,
)
from fem.tower_loading import (
    create_asymmetric_cyclic_top_force_history,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestNonlinearTowerSnapshotTensor(unittest.TestCase):
    """Verify exact reorganization and persistence of tower snapshots."""

    @classmethod
    def setUpClass(cls) -> None:
        cls.configuration = TowerConfiguration(
            horizontal_force=1.0e6,
            n_elements=2,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        cls.material = MaterialParameters()
        cls.loading = create_asymmetric_cyclic_top_force_history(
            maximum_force=1.0e6,
            force_ratio=-0.5,
            period=10.0,
            n_cycles=2,
            increments_per_cycle=8,
        )
        cls.response = run_nonlinear_asymmetric_analysis(
            configuration=cls.configuration,
            material=cls.material,
            loading=cls.loading,
            max_iterations=40,
        )
        cls.snapshots = build_tower_cycle_phase_snapshots(
            cls.response
        )

    def test_cycle_phase_grid_shape(self) -> None:
        snapshots = self.snapshots

        self.assertEqual(snapshots.n_cycles, 2)
        self.assertEqual(snapshots.n_phase_points, 9)
        self.assertEqual(
            snapshots.analysis_times.shape,
            (2, 9),
        )
        self.assertEqual(
            snapshots.nodal_displacements.shape,
            (2, 9, 9),
        )
        self.assertEqual(
            snapshots.fiber_strains.shape[:2],
            (2, 9),
        )
        self.assertEqual(
            snapshots.fiber_stresses.shape,
            snapshots.fiber_strains.shape,
        )
        self.assertEqual(
            snapshots.fiber_states.shape,
            snapshots.fiber_strains.shape + (4,),
        )

    def test_cycle_numbers_and_fast_phase_are_exact(self) -> None:
        snapshots = self.snapshots

        np.testing.assert_array_equal(
            snapshots.cycle_numbers,
            np.array([1, 2], dtype=np.int64),
        )
        np.testing.assert_allclose(
            snapshots.phase_times,
            np.linspace(0.0, 10.0, 9),
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            snapshots.phase_fractions,
            np.linspace(0.0, 1.0, 9),
            rtol=0.0,
            atol=1.0e-14,
        )

    def test_each_cycle_matches_original_time_history(self) -> None:
        snapshots = self.snapshots
        response = self.response

        np.testing.assert_allclose(
            snapshots.analysis_times[0],
            response.analysis_times[0:9],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.analysis_times[1],
            response.analysis_times[8:17],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.nodal_displacements[0],
            response.nodal_displacements[0:9],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.nodal_displacements[1],
            response.nodal_displacements[8:17],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_strains[0],
            response.fiber_strains[0:9],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_strains[1],
            response.fiber_strains[8:17],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_stresses[0],
            response.fiber_stresses[0:9],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_stresses[1],
            response.fiber_stresses[8:17],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_states[0],
            response.fiber_states[0:9],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_states[1],
            response.fiber_states[8:17],
            rtol=0.0,
            atol=0.0,
        )

    def test_shared_cycle_boundary_is_duplicated_exactly(self) -> None:
        snapshots = self.snapshots

        self.assertEqual(
            float(snapshots.analysis_times[0, -1]),
            float(snapshots.analysis_times[1, 0]),
        )
        np.testing.assert_allclose(
            snapshots.nodal_displacements[0, -1],
            snapshots.nodal_displacements[1, 0],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_strains[0, -1],
            snapshots.fiber_strains[1, 0],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_stresses[0, -1],
            snapshots.fiber_stresses[1, 0],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_states[0, -1],
            snapshots.fiber_states[1, 0],
            rtol=0.0,
            atol=0.0,
        )

    def test_phase_force_contains_exact_asymmetric_checkpoints(self) -> None:
        snapshots = self.snapshots
        loading = self.loading
        quarter = loading.increments_per_cycle // 4

        np.testing.assert_allclose(
            snapshots.phase_forces,
            loading.forces[:9],
            rtol=0.0,
            atol=0.0,
        )
        self.assertAlmostEqual(
            float(snapshots.phase_forces[0]),
            loading.mean_force,
            places=8,
        )
        self.assertAlmostEqual(
            float(snapshots.phase_forces[quarter]),
            loading.maximum_force,
            places=8,
        )
        self.assertAlmostEqual(
            float(snapshots.phase_forces[2 * quarter]),
            loading.mean_force,
            places=8,
        )
        self.assertAlmostEqual(
            float(snapshots.phase_forces[3 * quarter]),
            loading.minimum_force,
            places=8,
        )
        self.assertAlmostEqual(
            float(snapshots.phase_forces[4 * quarter]),
            loading.mean_force,
            places=8,
        )

    def test_material_state_channels_are_preserved(self) -> None:
        snapshots = self.snapshots

        np.testing.assert_allclose(
            snapshots.fiber_plastic_strains,
            snapshots.fiber_states[..., 0],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_alphas,
            snapshots.fiber_states[..., 1],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_r_bars,
            snapshots.fiber_states[..., 2],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            snapshots.fiber_damages,
            snapshots.fiber_states[..., 3],
            rtol=0.0,
            atol=0.0,
        )


    def test_global_cycle_numbers_may_start_above_one(self) -> None:
        """Stage snapshots must retain their original global cycle numbers."""
        snapshots = self.snapshots

        shifted = TowerCyclePhaseSnapshots(
            cycle_numbers=np.array(
                [21, 22],
                dtype=np.int64,
            ),
            phase_times=snapshots.phase_times,
            phase_fractions=snapshots.phase_fractions,
            phase_forces=snapshots.phase_forces,
            analysis_times=snapshots.analysis_times,
            nodal_displacements=snapshots.nodal_displacements,
            fiber_strains=snapshots.fiber_strains,
            fiber_stresses=snapshots.fiber_stresses,
            fiber_states=snapshots.fiber_states,
        )

        np.testing.assert_array_equal(
            shifted.cycle_numbers,
            np.array(
                [21, 22],
                dtype=np.int64,
            ),
        )
        self.assertEqual(
            shifted.n_cycles,
            2,
        )

    def test_global_cycle_numbers_must_be_positive_and_consecutive(self) -> None:
        """Non-positive, gapped, or reversed global numbering is invalid."""
        snapshots = self.snapshots

        invalid_cycle_numbers = (
            np.array([0, 1], dtype=np.int64),
            np.array([21, 23], dtype=np.int64),
            np.array([22, 21], dtype=np.int64),
        )

        for cycle_numbers in invalid_cycle_numbers:
            with self.subTest(
                cycle_numbers=cycle_numbers.tolist()
            ):
                with self.assertRaises(ValueError):
                    TowerCyclePhaseSnapshots(
                        cycle_numbers=cycle_numbers,
                        phase_times=snapshots.phase_times,
                        phase_fractions=snapshots.phase_fractions,
                        phase_forces=snapshots.phase_forces,
                        analysis_times=snapshots.analysis_times,
                        nodal_displacements=(
                            snapshots.nodal_displacements
                        ),
                        fiber_strains=snapshots.fiber_strains,
                        fiber_stresses=snapshots.fiber_stresses,
                        fiber_states=snapshots.fiber_states,
                    )

    def test_select_cycle_range_preserves_global_numbering_and_data(self) -> None:
        """Stage selection must not renumber cycles locally."""
        snapshots = self.snapshots

        shifted = TowerCyclePhaseSnapshots(
            cycle_numbers=np.array(
                [21, 22],
                dtype=np.int64,
            ),
            phase_times=snapshots.phase_times,
            phase_fractions=snapshots.phase_fractions,
            phase_forces=snapshots.phase_forces,
            analysis_times=snapshots.analysis_times,
            nodal_displacements=snapshots.nodal_displacements,
            fiber_strains=snapshots.fiber_strains,
            fiber_stresses=snapshots.fiber_stresses,
            fiber_states=snapshots.fiber_states,
        )

        selected = select_tower_cycle_range(
            snapshots=shifted,
            first_cycle=22,
            last_cycle=22,
        )

        np.testing.assert_array_equal(
            selected.cycle_numbers,
            np.array(
                [22],
                dtype=np.int64,
            ),
        )
        self.assertEqual(
            selected.n_cycles,
            1,
        )
        np.testing.assert_array_equal(
            selected.phase_times,
            shifted.phase_times,
        )
        np.testing.assert_array_equal(
            selected.phase_fractions,
            shifted.phase_fractions,
        )
        np.testing.assert_array_equal(
            selected.phase_forces,
            shifted.phase_forces,
        )
        np.testing.assert_array_equal(
            selected.analysis_times,
            shifted.analysis_times[1:2],
        )
        np.testing.assert_array_equal(
            selected.nodal_displacements,
            shifted.nodal_displacements[1:2],
        )
        np.testing.assert_array_equal(
            selected.fiber_strains,
            shifted.fiber_strains[1:2],
        )
        np.testing.assert_array_equal(
            selected.fiber_stresses,
            shifted.fiber_stresses[1:2],
        )
        np.testing.assert_array_equal(
            selected.fiber_states,
            shifted.fiber_states[1:2],
        )

    def test_select_cycle_range_rejects_invalid_bounds(self) -> None:
        snapshots = self.snapshots

        invalid_ranges = (
            (0, 1),
            (2, 1),
            (1, 3),
        )

        for first_cycle, last_cycle in invalid_ranges:
            with self.subTest(
                first_cycle=first_cycle,
                last_cycle=last_cycle,
            ):
                with self.assertRaises(ValueError):
                    select_tower_cycle_range(
                        snapshots=snapshots,
                        first_cycle=first_cycle,
                        last_cycle=last_cycle,
                    )

        with self.assertRaises(TypeError):
            select_tower_cycle_range(
                snapshots=snapshots,
                first_cycle=True,
                last_cycle=2,
            )

    def test_npz_round_trip_preserves_global_cycle_numbers(self) -> None:
        """Persistence must retain non-one-based global cycle numbering."""
        snapshots = self.snapshots

        shifted = TowerCyclePhaseSnapshots(
            cycle_numbers=np.array(
                [47, 48],
                dtype=np.int64,
            ),
            phase_times=snapshots.phase_times,
            phase_fractions=snapshots.phase_fractions,
            phase_forces=snapshots.phase_forces,
            analysis_times=snapshots.analysis_times,
            nodal_displacements=snapshots.nodal_displacements,
            fiber_strains=snapshots.fiber_strains,
            fiber_stresses=snapshots.fiber_stresses,
            fiber_states=snapshots.fiber_states,
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            file_path = (
                Path(temporary_directory)
                / "global_cycle_snapshots.npz"
            )

            save_tower_cycle_phase_snapshots(
                snapshots=shifted,
                file_path=file_path,
            )
            loaded = load_tower_cycle_phase_snapshots(
                file_path=file_path,
            )

            np.testing.assert_array_equal(
                loaded.cycle_numbers,
                np.array(
                    [47, 48],
                    dtype=np.int64,
                ),
            )
            np.testing.assert_array_equal(
                loaded.fiber_states,
                shifted.fiber_states,
            )

    def test_npz_round_trip_preserves_all_snapshot_arrays_exactly(self) -> None:
        """Saving and loading must reproduce the snapshot object exactly."""
        snapshots = self.snapshots

        with tempfile.TemporaryDirectory() as temporary_directory:
            file_path = (
                Path(temporary_directory)
                / "tower_cycle_phase_snapshots.npz"
            )

            saved_path = save_tower_cycle_phase_snapshots(
                snapshots=snapshots,
                file_path=file_path,
            )
            loaded = load_tower_cycle_phase_snapshots(
                file_path=saved_path,
            )

            self.assertEqual(saved_path, file_path)
            self.assertTrue(saved_path.is_file())
            self.assertIsInstance(
                loaded,
                TowerCyclePhaseSnapshots,
            )

            np.testing.assert_array_equal(
                loaded.cycle_numbers,
                snapshots.cycle_numbers,
            )
            np.testing.assert_array_equal(
                loaded.phase_times,
                snapshots.phase_times,
            )
            np.testing.assert_array_equal(
                loaded.phase_fractions,
                snapshots.phase_fractions,
            )
            np.testing.assert_array_equal(
                loaded.phase_forces,
                snapshots.phase_forces,
            )
            np.testing.assert_array_equal(
                loaded.analysis_times,
                snapshots.analysis_times,
            )
            np.testing.assert_array_equal(
                loaded.nodal_displacements,
                snapshots.nodal_displacements,
            )
            np.testing.assert_array_equal(
                loaded.fiber_strains,
                snapshots.fiber_strains,
            )
            np.testing.assert_array_equal(
                loaded.fiber_stresses,
                snapshots.fiber_stresses,
            )
            np.testing.assert_array_equal(
                loaded.fiber_states,
                snapshots.fiber_states,
            )

    def test_npz_contains_only_version_and_primitive_snapshot_arrays(self) -> None:
        """The frozen FOM file must not duplicate derived diagnostics."""
        snapshots = self.snapshots

        expected_keys = {
            "snapshot_format_version",
            "cycle_numbers",
            "phase_times",
            "phase_fractions",
            "phase_forces",
            "analysis_times",
            "nodal_displacements",
            "fiber_strains",
            "fiber_stresses",
            "fiber_states",
        }

        with tempfile.TemporaryDirectory() as temporary_directory:
            file_path = (
                Path(temporary_directory)
                / "tower_cycle_phase_snapshots.npz"
            )

            save_tower_cycle_phase_snapshots(
                snapshots=snapshots,
                file_path=file_path,
            )

            with np.load(
                file_path,
                allow_pickle=False,
            ) as archive:
                self.assertEqual(
                    set(archive.files),
                    expected_keys,
                )
                self.assertEqual(
                    int(
                        np.asarray(
                            archive["snapshot_format_version"]
                        ).reshape(-1)[0]
                    ),
                    1,
                )


if __name__ == "__main__":
    unittest.main()

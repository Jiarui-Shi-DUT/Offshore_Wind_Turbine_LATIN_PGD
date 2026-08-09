# -*- coding: utf-8 -*-
"""Regression tests for slow-cycle x fast-phase tower snapshots."""

from __future__ import annotations

import unittest

import numpy as np

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_response import (
    run_nonlinear_asymmetric_analysis,
)
from examples.nonlinear_tower_snapshot_tensor import (
    build_tower_cycle_phase_snapshots,
)
from fem.tower_loading import (
    create_asymmetric_cyclic_top_force_history,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestNonlinearTowerSnapshotTensor(unittest.TestCase):
    """Verify exact reorganization of FOM histories into (n, tau) form."""

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


if __name__ == "__main__":
    unittest.main()

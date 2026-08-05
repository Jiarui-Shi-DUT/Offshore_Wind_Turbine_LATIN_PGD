# -*- coding: utf-8 -*-
"""Regression tests for sequential nonlinear tower load histories."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.tower_loading import (
    create_pulsating_top_force_history,
)
from fem.tower_system_2d import (
    assemble_elastic_tower_stiffness,
    cantilever_base_fixed_dofs,
    solve_linear_static_system,
    top_horizontal_load_vector,
)
from fem.viscoplastic_tower_history_2d import (
    NonlinearPulsatingTowerHistory,
    solve_pulsating_tower_history,
)
from fem.viscoplastic_tower_system_2d import (
    NonlinearTowerConvergenceError,
    ViscoplasticDamageTowerSystem2D,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestViscoplasticTowerHistory2D(unittest.TestCase):
    """Verify sequential loading, storage, continuation, and rollback."""

    def setUp(self) -> None:
        self.height = 2.0
        self.mesh = create_uniform_vertical_tower_mesh(
            height=self.height,
            n_elements=2,
        )
        self.geometry = LinearTaperedTowerGeometry(
            height=self.height,
            base_outer_diameter=1.0,
            top_outer_diameter=1.0,
            base_thickness=0.1,
            top_thickness=0.1,
        )
        self.material = MaterialParameters(
            sigma_y=1.0e9,
            k_damage=0.0,
        )

    def create_system(self) -> ViscoplasticDamageTowerSystem2D:
        """Create a small elastic-limit nonlinear tower."""
        return ViscoplasticDamageTowerSystem2D(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            material=self.material,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )

    def create_loading(self):
        """Create one coarsely sampled pulsating cycle."""
        return create_pulsating_top_force_history(
            maximum_force=1.0e5,
            force_ratio=0.1,
            period=2.0,
            n_cycles=1,
            increments_per_cycle=4,
        )

    def test_history_shapes_times_and_preload(self) -> None:
        """Stored arrays must match the loading and preload definition."""
        system = self.create_system()
        loading = self.create_loading()

        history = solve_pulsating_tower_history(
            system=system,
            loading=loading,
        )

        self.assertIsInstance(
            history,
            NonlinearPulsatingTowerHistory,
        )
        self.assertEqual(
            history.displacements.shape,
            (loading.n_time_points, self.mesh.n_dof),
        )
        self.assertEqual(
            history.reactions.shape,
            history.displacements.shape,
        )
        self.assertEqual(
            history.iterations.shape,
            loading.times.shape,
        )
        self.assertEqual(
            history.maximum_damages.shape,
            loading.times.shape,
        )

        expected_preload = loading.time_increment
        expected_analysis_times = (
            expected_preload + loading.times
        )
        self.assertAlmostEqual(
            history.preload_duration,
            expected_preload,
        )
        np.testing.assert_allclose(
            history.analysis_times,
            expected_analysis_times,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            history.physical_times,
            loading.times,
            rtol=0.0,
            atol=0.0,
        )
        self.assertAlmostEqual(
            system.committed_time,
            expected_analysis_times[-1],
        )
        self.assertFalse(system.has_uncommitted_trial)

    def test_elastic_history_matches_independent_linear_solutions(
        self,
    ) -> None:
        """Every sequential nonlinear point must recover the elastic FE result."""
        system = self.create_system()
        loading = self.create_loading()

        history = solve_pulsating_tower_history(
            system=system,
            loading=loading,
            relative_residual_tolerance=1.0e-9,
            absolute_residual_tolerance=1.0e-6,
        )

        elastic = assemble_elastic_tower_stiffness(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            elastic_modulus=self.material.E * 1.0e6,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        fixed = cantilever_base_fixed_dofs(self.mesh)

        expected_displacements = []
        expected_reactions = []

        for force in loading.forces:
            load_vector = top_horizontal_load_vector(
                mesh=self.mesh,
                horizontal_force=float(force),
            )
            solution = solve_linear_static_system(
                global_stiffness=elastic.global_stiffness,
                load_vector=load_vector,
                fixed_dofs=fixed,
            )
            expected_displacements.append(
                solution.displacements
            )
            expected_reactions.append(solution.reactions)

        np.testing.assert_allclose(
            history.displacements,
            np.stack(expected_displacements, axis=0),
            rtol=1.0e-7,
            atol=1.0e-11,
        )
        np.testing.assert_allclose(
            history.reactions,
            np.stack(expected_reactions, axis=0),
            rtol=1.0e-7,
            atol=1.0e-5,
        )
        np.testing.assert_allclose(
            history.maximum_damages,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        self.assertLessEqual(
            int(np.max(history.iterations)),
            3,
        )

    def test_convenience_response_histories(self) -> None:
        """Convenience properties must select the correct global DOFs."""
        system = self.create_system()
        history = solve_pulsating_tower_history(
            system=system,
            loading=self.create_loading(),
        )

        np.testing.assert_allclose(
            history.top_horizontal_displacements,
            history.displacements[:, -3],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            history.top_rotations,
            history.displacements[:, -1],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            history.base_horizontal_reactions,
            history.reactions[:, 0],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            history.base_bending_reactions,
            history.reactions[:, 2],
            rtol=0.0,
            atol=0.0,
        )
        self.assertEqual(history.final_maximum_damage, 0.0)

    def test_custom_preload_duration_is_used(self) -> None:
        """A user-specified preload duration must shift analysis time only."""
        system = self.create_system()
        loading = self.create_loading()
        preload_duration = 0.125

        history = solve_pulsating_tower_history(
            system=system,
            loading=loading,
            preload_duration=preload_duration,
        )

        self.assertAlmostEqual(
            history.preload_duration,
            preload_duration,
        )
        np.testing.assert_allclose(
            history.analysis_times,
            preload_duration + loading.times,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            history.physical_times,
            loading.times,
            rtol=0.0,
            atol=0.0,
        )

    def test_history_can_continue_from_committed_state(self) -> None:
        """A second history must start after the first committed endpoint."""
        system = self.create_system()
        loading = self.create_loading()

        first = solve_pulsating_tower_history(
            system=system,
            loading=loading,
        )
        first_end_time = float(first.analysis_times[-1])

        second = solve_pulsating_tower_history(
            system=system,
            loading=loading,
        )

        expected_second_times = (
            first_end_time
            + loading.time_increment
            + loading.times
        )
        np.testing.assert_allclose(
            second.analysis_times,
            expected_second_times,
            rtol=0.0,
            atol=0.0,
        )
        self.assertGreater(
            second.analysis_times[0],
            first.analysis_times[-1],
        )
        self.assertAlmostEqual(
            system.committed_time,
            second.analysis_times[-1],
        )

    def test_failed_history_step_preserves_last_commit(self) -> None:
        """Failure of the first history step must leave the initial state intact."""
        system = self.create_system()
        loading = self.create_loading()

        with self.assertRaises(
            NonlinearTowerConvergenceError
        ):
            solve_pulsating_tower_history(
                system=system,
                loading=loading,
                max_iterations=1,
            )

        self.assertEqual(system.committed_time, 0.0)
        self.assertFalse(system.has_uncommitted_trial)
        np.testing.assert_allclose(
            system.committed_displacements,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            system.committed_response.internal_force,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

    def test_invalid_history_arguments_are_rejected(self) -> None:
        """Invalid preload, iteration count, and object types must fail."""
        system = self.create_system()
        loading = self.create_loading()

        with self.assertRaises(TypeError):
            solve_pulsating_tower_history(
                system=object(),
                loading=loading,
            )

        with self.assertRaises(TypeError):
            solve_pulsating_tower_history(
                system=system,
                loading=object(),
            )

        with self.assertRaises(ValueError):
            solve_pulsating_tower_history(
                system=system,
                loading=loading,
                preload_duration=0.0,
            )

        with self.assertRaises(ValueError):
            solve_pulsating_tower_history(
                system=system,
                loading=loading,
                max_iterations=0,
            )


if __name__ == "__main__":
    unittest.main()

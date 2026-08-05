# -*- coding: utf-8 -*-
"""Regression tests for the nonlinear viscoplastic tower system."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.tower_system_2d import (
    assemble_elastic_tower_stiffness,
    cantilever_base_fixed_dofs,
    solve_linear_static_system,
    top_horizontal_load_vector,
)
from fem.viscoplastic_tower_system_2d import (
    NonlinearTowerConvergenceError,
    ViscoplasticDamageTowerSystem2D,
    solve_nonlinear_tower_load_step,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestViscoplasticDamageTowerSystem2D(unittest.TestCase):
    """Verify global assembly, Newton solution, and state management."""

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
        self.elastic_material = MaterialParameters(
            sigma_y=1.0e9,
            k_damage=0.0,
        )

    def create_system(
        self,
        material: MaterialParameters = None,
    ) -> ViscoplasticDamageTowerSystem2D:
        """Create a small two-element tower system."""
        if material is None:
            material = self.elastic_material

        return ViscoplasticDamageTowerSystem2D(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            material=material,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )

    @staticmethod
    def stacked_trial_states(
        system: ViscoplasticDamageTowerSystem2D,
    ) -> np.ndarray:
        """Return all element, Gauss-point, and fiber trial states."""
        return np.stack(
            [
                np.stack(
                    [
                        section.trial_states
                        for section in element.sections
                    ],
                    axis=0,
                )
                for element in system.elements
            ],
            axis=0,
        )

    @staticmethod
    def stacked_committed_states(
        system: ViscoplasticDamageTowerSystem2D,
    ) -> np.ndarray:
        """Return all committed material states in the tower."""
        return np.stack(
            [
                np.stack(
                    [
                        section.committed_states
                        for section in element.sections
                    ],
                    axis=0,
                )
                for element in system.elements
            ],
            axis=0,
        )

    def test_initial_global_system_state(self) -> None:
        """The tower must start with zero response and valid DOF maps."""
        system = self.create_system()
        response = system.committed_response

        self.assertEqual(len(system.elements), self.mesh.n_elements)
        self.assertEqual(system.element_dofs.shape, (2, 6))
        self.assertEqual(
            response.displacements.shape,
            (self.mesh.n_dof,),
        )
        self.assertEqual(
            response.internal_force.shape,
            (self.mesh.n_dof,),
        )
        self.assertEqual(len(response.element_responses), 2)
        self.assertFalse(system.has_uncommitted_trial)
        self.assertEqual(system.committed_time, 0.0)

        np.testing.assert_allclose(
            response.displacements,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            response.internal_force,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            self.stacked_committed_states(system),
            0.0,
            rtol=0.0,
            atol=0.0,
        )

        expected_dofs = np.array(
            [
                [0, 1, 2, 3, 4, 5],
                [3, 4, 5, 6, 7, 8],
            ],
            dtype=np.int64,
        )
        np.testing.assert_array_equal(
            system.element_dofs,
            expected_dofs,
        )

    def test_elastic_global_tangent_matches_existing_assembly(
        self,
    ) -> None:
        """The nonlinear assembly must recover the elastic global matrix."""
        system = self.create_system()

        response = system.set_trial_displacements(
            time=0.1,
            displacements=np.zeros(self.mesh.n_dof),
            compute_tangent=True,
        )
        elastic = assemble_elastic_tower_stiffness(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            elastic_modulus=(
                self.elastic_material.E * 1.0e6
            ),
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )

        self.assertIsNotNone(response.tangent)
        np.testing.assert_allclose(
            response.tangent,
            elastic.global_stiffness,
            rtol=1.0e-7,
            atol=1.0e-2,
        )
        np.testing.assert_allclose(
            response.tangent,
            response.tangent.T,
            rtol=0.0,
            atol=1.0e-8,
        )

        system.revert_to_last_commit()

    def test_elastic_newton_solution_matches_linear_solver(
        self,
    ) -> None:
        """One nonlinear load step must reproduce the elastic solution."""
        system = self.create_system()
        load_vector = top_horizontal_load_vector(
            mesh=self.mesh,
            horizontal_force=1.0e5,
        )
        fixed_dofs = cantilever_base_fixed_dofs(self.mesh)

        nonlinear = solve_nonlinear_tower_load_step(
            system=system,
            time=0.1,
            load_vector=load_vector,
            fixed_dofs=fixed_dofs,
            max_iterations=5,
            relative_residual_tolerance=1.0e-9,
            absolute_residual_tolerance=1.0e-6,
        )

        elastic = assemble_elastic_tower_stiffness(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            elastic_modulus=(
                self.elastic_material.E * 1.0e6
            ),
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        linear = solve_linear_static_system(
            global_stiffness=elastic.global_stiffness,
            load_vector=load_vector,
            fixed_dofs=fixed_dofs,
        )

        np.testing.assert_allclose(
            nonlinear.displacements,
            linear.displacements,
            rtol=1.0e-7,
            atol=1.0e-11,
        )
        np.testing.assert_allclose(
            nonlinear.reactions,
            linear.reactions,
            rtol=1.0e-7,
            atol=1.0e-5,
        )
        self.assertLessEqual(nonlinear.iterations, 3)
        self.assertLess(nonlinear.residual_norm, 1.0e-3)
        self.assertEqual(system.committed_time, 0.1)
        self.assertFalse(system.has_uncommitted_trial)

    def test_repeated_global_trial_does_not_accumulate_damage(
        self,
    ) -> None:
        """Repeated Newton trials must restart from the last commit."""
        system = self.create_system(MaterialParameters())
        displacements = np.zeros(self.mesh.n_dof)
        top_vertical_dof = 3 * (self.mesh.n_nodes - 1) + 1
        displacements[top_vertical_dof] = (
            1.2e-3 * self.height
        )

        first = system.set_trial_displacements(
            time=0.1,
            displacements=displacements,
            compute_tangent=False,
        )
        first_states = self.stacked_trial_states(system)

        second = system.set_trial_displacements(
            time=0.1,
            displacements=displacements,
            compute_tangent=False,
        )
        second_states = self.stacked_trial_states(system)

        np.testing.assert_allclose(
            second_states,
            first_states,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            second.internal_force,
            first.internal_force,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            self.stacked_committed_states(system),
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        self.assertGreater(first.maximum_damage, 0.0)

    def test_global_commit_rollback_and_restart(self) -> None:
        """Global state operations must propagate through all elements."""
        system = self.create_system(MaterialParameters())
        first_displacements = np.zeros(self.mesh.n_dof)
        top_vertical_dof = 3 * (self.mesh.n_nodes - 1) + 1
        first_displacements[top_vertical_dof] = (
            1.2e-3 * self.height
        )

        system.set_trial_displacements(
            time=0.1,
            displacements=first_displacements,
            compute_tangent=False,
        )
        committed = system.commit_state()
        committed_states = self.stacked_committed_states(
            system
        )

        self.assertGreater(committed.maximum_damage, 0.0)
        self.assertFalse(system.has_uncommitted_trial)

        second_displacements = np.zeros(self.mesh.n_dof)
        second_displacements[top_vertical_dof] = (
            -1.2e-3 * self.height
        )
        system.set_trial_displacements(
            time=0.2,
            displacements=second_displacements,
            compute_tangent=False,
        )
        self.assertTrue(system.has_uncommitted_trial)

        system.revert_to_last_commit()

        self.assertFalse(system.has_uncommitted_trial)
        np.testing.assert_allclose(
            system.trial_displacements,
            first_displacements,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            self.stacked_trial_states(system),
            committed_states,
            rtol=0.0,
            atol=0.0,
        )

        restarted = system.revert_to_start()

        self.assertEqual(system.committed_time, 0.0)
        self.assertFalse(system.has_uncommitted_trial)
        np.testing.assert_allclose(
            restarted.displacements,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            self.stacked_committed_states(system),
            0.0,
            rtol=0.0,
            atol=0.0,
        )

    def test_failed_newton_step_rolls_back_all_states(self) -> None:
        """A non-converged step must not alter committed tower history."""
        system = self.create_system()
        load_vector = top_horizontal_load_vector(
            mesh=self.mesh,
            horizontal_force=1.0e5,
        )

        with self.assertRaises(
            NonlinearTowerConvergenceError
        ):
            solve_nonlinear_tower_load_step(
                system=system,
                time=0.1,
                load_vector=load_vector,
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
            self.stacked_committed_states(system),
            0.0,
            rtol=0.0,
            atol=0.0,
        )

    def test_invalid_global_operations_are_rejected(self) -> None:
        """Invalid vector sizes, times, and empty commits must fail."""
        system = self.create_system()

        with self.assertRaises(RuntimeError):
            system.commit_state()

        with self.assertRaises(ValueError):
            system.set_trial_displacements(
                time=0.1,
                displacements=np.zeros(
                    self.mesh.n_dof - 1
                ),
            )

        with self.assertRaises(ValueError):
            solve_nonlinear_tower_load_step(
                system=system,
                time=0.0,
                load_vector=np.zeros(self.mesh.n_dof),
            )

        with self.assertRaises(ValueError):
            solve_nonlinear_tower_load_step(
                system=system,
                time=0.1,
                load_vector=np.zeros(
                    self.mesh.n_dof - 1
                ),
            )


if __name__ == "__main__":
    unittest.main()

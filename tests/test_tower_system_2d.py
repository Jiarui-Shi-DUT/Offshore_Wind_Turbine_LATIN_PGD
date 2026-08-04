# -*- coding: utf-8 -*-
"""Regression tests for global assembly and static tower solution."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.fiber_section import create_annular_fiber_section
from fem.tower_system_2d import (
    assemble_elastic_tower_stiffness,
    cantilever_base_fixed_dofs,
    element_dof_indices,
    free_dofs_from_fixed,
    solve_elastic_tower_top_force,
    solve_linear_static_system,
    top_horizontal_load_vector,
)


class TestTowerSystem2D(unittest.TestCase):
    """Verify assembly, boundary conditions, loading, and static solution."""

    def setUp(self) -> None:
        self.elastic_modulus = 210.0e9
        self.force = 1.0e6

    @staticmethod
    def constant_tower(
        height: float,
        outer_diameter: float = 6.0,
        thickness: float = 0.027,
    ) -> LinearTaperedTowerGeometry:
        """Return a constant annular tower represented by the taper class."""
        return LinearTaperedTowerGeometry(
            height=height,
            base_outer_diameter=outer_diameter,
            top_outer_diameter=outer_diameter,
            base_thickness=thickness,
            top_thickness=thickness,
        )

    def test_element_dof_indices(self) -> None:
        dofs = element_dof_indices([2, 5])

        np.testing.assert_array_equal(
            dofs,
            np.array([6, 7, 8, 15, 16, 17], dtype=np.int64),
        )

        with self.assertRaises(ValueError):
            element_dof_indices([1, 1])

        with self.assertRaises(ValueError):
            element_dof_indices([-1, 2])

    def test_boundary_conditions_and_top_load(self) -> None:
        mesh = create_uniform_vertical_tower_mesh(
            height=10.0,
            n_elements=4,
        )

        fixed = cantilever_base_fixed_dofs(mesh)
        free = free_dofs_from_fixed(mesh.n_dof, fixed)
        loads = top_horizontal_load_vector(
            mesh=mesh,
            horizontal_force=self.force,
        )

        np.testing.assert_array_equal(
            fixed,
            np.array([0, 1, 2], dtype=np.int64),
        )
        self.assertEqual(free.size, 12)
        self.assertEqual(loads.shape, (15,))
        self.assertEqual(
            int(np.count_nonzero(loads)),
            1,
        )
        self.assertAlmostEqual(
            loads[12],
            self.force,
            places=10,
        )
        self.assertAlmostEqual(
            float(np.sum(loads)),
            self.force,
            places=10,
        )

    def test_global_assembly_shape_and_symmetry(self) -> None:
        mesh = create_uniform_vertical_tower_mesh(
            height=12.0,
            n_elements=3,
        )
        tower = self.constant_tower(height=12.0)

        assembly = assemble_elastic_tower_stiffness(
            mesh=mesh,
            tower_geometry=tower,
            elastic_modulus=self.elastic_modulus,
        )

        self.assertEqual(
            assembly.global_stiffness.shape,
            (12, 12),
        )
        self.assertEqual(
            assembly.element_dofs.shape,
            (3, 6),
        )
        self.assertEqual(
            len(assembly.element_results),
            3,
        )
        np.testing.assert_allclose(
            assembly.global_stiffness,
            assembly.global_stiffness.T,
            rtol=0.0,
            atol=1.0e-7,
        )
        self.assertTrue(
            np.all(np.isfinite(assembly.global_stiffness))
        )

    def test_single_element_constant_section_cantilever(self) -> None:
        height = 10.0
        outer_diameter = 6.0
        thickness = 0.027

        mesh = create_uniform_vertical_tower_mesh(
            height=height,
            n_elements=1,
        )
        tower = self.constant_tower(
            height=height,
            outer_diameter=outer_diameter,
            thickness=thickness,
        )

        _, solution = solve_elastic_tower_top_force(
            mesh=mesh,
            tower_geometry=tower,
            elastic_modulus=self.elastic_modulus,
            horizontal_force=self.force,
        )

        section = create_annular_fiber_section(
            outer_diameter=outer_diameter,
            thickness=thickness,
            n_circumferential=32,
            n_radial=2,
        )
        bending_rigidity = (
            self.elastic_modulus * section.second_moment_x
        )

        expected_tip_displacement = (
            self.force * height ** 3
            / (3.0 * bending_rigidity)
        )
        expected_tip_rotation = (
            -self.force * height ** 2
            / (2.0 * bending_rigidity)
        )

        np.testing.assert_allclose(
            solution.displacements[:3],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        self.assertAlmostEqual(
            solution.displacements[3],
            expected_tip_displacement,
            delta=1.0e-11
            * max(1.0, abs(expected_tip_displacement)),
        )
        self.assertAlmostEqual(
            solution.displacements[4],
            0.0,
            delta=1.0e-14,
        )
        self.assertAlmostEqual(
            solution.displacements[5],
            expected_tip_rotation,
            delta=1.0e-11
            * max(1.0, abs(expected_tip_rotation)),
        )

        expected_base_reactions = np.array(
            [
                -self.force,
                0.0,
                self.force * height,
            ],
            dtype=np.float64,
        )
        np.testing.assert_allclose(
            solution.reactions[:3],
            expected_base_reactions,
            rtol=1.0e-11,
            atol=1.0e-5,
        )
        np.testing.assert_allclose(
            solution.reactions[3:],
            0.0,
            rtol=0.0,
            atol=1.0e-5,
        )

    def test_multi_element_constant_section_matches_cantilever_solution(
        self,
    ) -> None:
        height = 10.0
        outer_diameter = 6.0
        thickness = 0.027

        mesh = create_uniform_vertical_tower_mesh(
            height=height,
            n_elements=8,
        )
        tower = self.constant_tower(
            height=height,
            outer_diameter=outer_diameter,
            thickness=thickness,
        )

        _, solution = solve_elastic_tower_top_force(
            mesh=mesh,
            tower_geometry=tower,
            elastic_modulus=self.elastic_modulus,
            horizontal_force=self.force,
        )

        section = create_annular_fiber_section(
            outer_diameter=outer_diameter,
            thickness=thickness,
            n_circumferential=32,
            n_radial=2,
        )
        bending_rigidity = (
            self.elastic_modulus * section.second_moment_x
        )
        expected_tip_displacement = (
            self.force * height ** 3
            / (3.0 * bending_rigidity)
        )
        expected_tip_rotation = (
            -self.force * height ** 2
            / (2.0 * bending_rigidity)
        )

        np.testing.assert_allclose(
            solution.displacements[-3],
            expected_tip_displacement,
            rtol=1.0e-9,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            solution.displacements[-2],
            0.0,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            solution.displacements[-1],
            expected_tip_rotation,
            rtol=1.0e-9,
            atol=1.0e-12,
        )

    def test_tapered_nrel_tower_solution_is_in_equilibrium(self) -> None:
        height = 87.6
        mesh = create_uniform_vertical_tower_mesh(
            height=height,
            n_elements=40,
        )
        tower = LinearTaperedTowerGeometry(
            height=height,
            base_outer_diameter=6.0,
            top_outer_diameter=3.87,
            base_thickness=0.027,
            top_thickness=0.019,
        )

        assembly, solution = solve_elastic_tower_top_force(
            mesh=mesh,
            tower_geometry=tower,
            elastic_modulus=self.elastic_modulus,
            horizontal_force=self.force,
        )

        self.assertEqual(
            assembly.global_stiffness.shape,
            (123, 123),
        )
        self.assertTrue(
            np.all(np.isfinite(solution.displacements))
        )
        self.assertTrue(
            np.all(np.isfinite(solution.reactions))
        )
        self.assertGreater(
            solution.displacements[-3],
            0.0,
        )

        expected_base_reactions = np.array(
            [
                -self.force,
                0.0,
                self.force * height,
            ],
            dtype=np.float64,
        )
        np.testing.assert_allclose(
            solution.reactions[:3],
            expected_base_reactions,
            rtol=1.0e-8,
            atol=1.0e-2,
        )
        np.testing.assert_allclose(
            solution.reactions[solution.free_dofs],
            0.0,
            rtol=0.0,
            atol=1.0e-2,
        )

    def test_direct_linear_solver(self) -> None:
        stiffness = np.array(
            [
                [2.0, -1.0],
                [-1.0, 1.0],
            ],
            dtype=np.float64,
        )
        loads = np.array([0.0, 3.0], dtype=np.float64)

        solution = solve_linear_static_system(
            global_stiffness=stiffness,
            load_vector=loads,
            fixed_dofs=[0],
        )

        np.testing.assert_allclose(
            solution.displacements,
            np.array([0.0, 3.0]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            solution.reactions,
            np.array([-3.0, 0.0]),
            rtol=0.0,
            atol=1.0e-15,
        )

    def test_invalid_system_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            free_dofs_from_fixed(
                n_dof=6,
                fixed_dofs=[0, 0],
            )

        with self.assertRaises(ValueError):
            solve_linear_static_system(
                global_stiffness=np.eye(3),
                load_vector=np.zeros(2),
                fixed_dofs=[0],
            )

        with self.assertRaises(np.linalg.LinAlgError):
            solve_linear_static_system(
                global_stiffness=np.zeros((2, 2)),
                load_vector=np.ones(2),
                fixed_dofs=[],
            )


if __name__ == "__main__":
    unittest.main()

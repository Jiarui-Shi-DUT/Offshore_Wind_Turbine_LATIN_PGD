# -*- coding: utf-8 -*-
"""Regression tests for the two-dimensional Euler-Bernoulli beam module."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    compute_elastic_beam_element_stiffness,
    create_uniform_vertical_tower_mesh,
    euler_bernoulli_strain_displacement,
    gauss_legendre_rule,
    planar_frame_transformation,
)
from fem.fiber_section import create_annular_fiber_section


class TestBeamColumn2D(unittest.TestCase):
    """Verify tower geometry, mesh, kinematics, and element stiffness."""

    def setUp(self) -> None:
        self.elastic_modulus = 210.0e9
        self.tower = LinearTaperedTowerGeometry(
            height=87.6,
            base_outer_diameter=6.0,
            top_outer_diameter=3.87,
            base_thickness=0.027,
            top_thickness=0.019,
        )

    def test_linear_tower_geometry(self) -> None:
        self.assertAlmostEqual(
            self.tower.outer_diameter_at(0.0),
            6.0,
            places=15,
        )
        self.assertAlmostEqual(
            self.tower.outer_diameter_at(87.6),
            3.87,
            places=15,
        )
        self.assertAlmostEqual(
            self.tower.outer_diameter_at(43.8),
            0.5 * (6.0 + 3.87),
            places=15,
        )

        self.assertAlmostEqual(
            self.tower.thickness_at(0.0),
            0.027,
            places=15,
        )
        self.assertAlmostEqual(
            self.tower.thickness_at(87.6),
            0.019,
            places=15,
        )
        self.assertAlmostEqual(
            self.tower.thickness_at(43.8),
            0.5 * (0.027 + 0.019),
            places=15,
        )

        self.assertGreater(
            self.tower.exact_area_at(0.0),
            self.tower.exact_area_at(87.6),
        )
        self.assertGreater(
            self.tower.exact_second_moment_at(0.0),
            self.tower.exact_second_moment_at(87.6),
        )

    def test_uniform_vertical_tower_mesh(self) -> None:
        mesh = create_uniform_vertical_tower_mesh(
            height=87.6,
            n_elements=40,
        )

        self.assertEqual(mesh.n_nodes, 41)
        self.assertEqual(mesh.n_elements, 40)
        self.assertEqual(mesh.n_dof, 123)
        self.assertEqual(mesh.coordinates.shape, (41, 2))
        self.assertEqual(mesh.connectivity.shape, (40, 2))

        np.testing.assert_allclose(
            mesh.coordinates[:, 0],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            mesh.coordinates[:, 1],
            np.linspace(0.0, 87.6, 41),
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            mesh.element_lengths,
            np.full(40, 87.6 / 40.0),
            rtol=1.0e-14,
            atol=1.0e-14,
        )

    def test_four_point_gauss_legendre_rule(self) -> None:
        coordinates, weights = gauss_legendre_rule(4)

        self.assertEqual(coordinates.shape, (4,))
        self.assertEqual(weights.shape, (4,))
        np.testing.assert_allclose(
            coordinates,
            -coordinates[::-1],
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            weights,
            weights[::-1],
            rtol=0.0,
            atol=1.0e-15,
        )
        self.assertAlmostEqual(
            float(np.sum(weights)),
            2.0,
            places=15,
        )

        for degree in range(8):
            numerical = float(
                np.dot(weights, coordinates ** degree)
            )
            if degree % 2 == 0:
                exact = 2.0 / (degree + 1)
            else:
                exact = 0.0

            self.assertAlmostEqual(
                numerical,
                exact,
                places=14,
            )

    def test_vertical_element_transformation(self) -> None:
        coordinates = np.array(
            [
                [0.0, 0.0],
                [0.0, 2.0],
            ],
            dtype=np.float64,
        )

        transformation = planar_frame_transformation(
            coordinates
        )

        expected_rotation = np.array(
            [
                [0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
            ],
            dtype=np.float64,
        )
        expected = np.zeros((6, 6), dtype=np.float64)
        expected[:3, :3] = expected_rotation
        expected[3:, 3:] = expected_rotation

        np.testing.assert_allclose(
            transformation,
            expected,
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            transformation @ transformation.T,
            np.eye(6),
            rtol=0.0,
            atol=1.0e-15,
        )

    def test_strain_displacement_rejects_rigid_modes(self) -> None:
        length = 3.2
        axial_translation = np.array(
            [1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
            dtype=np.float64,
        )
        transverse_translation = np.array(
            [0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
            dtype=np.float64,
        )
        rigid_rotation = np.array(
            [0.0, 0.0, 1.0, 0.0, length, 1.0],
            dtype=np.float64,
        )

        for natural_coordinate in (-1.0, -0.3, 0.0, 0.6, 1.0):
            strain_displacement = (
                euler_bernoulli_strain_displacement(
                    length=length,
                    natural_coordinate=natural_coordinate,
                )
            )

            np.testing.assert_allclose(
                strain_displacement @ axial_translation,
                0.0,
                rtol=0.0,
                atol=1.0e-15,
            )
            np.testing.assert_allclose(
                strain_displacement @ transverse_translation,
                0.0,
                rtol=0.0,
                atol=1.0e-15,
            )
            np.testing.assert_allclose(
                strain_displacement @ rigid_rotation,
                0.0,
                rtol=0.0,
                atol=1.0e-15,
            )

    def test_constant_section_stiffness_matches_closed_form(self) -> None:
        length = 4.0
        diameter = 6.0
        thickness = 0.027

        constant_tower = LinearTaperedTowerGeometry(
            height=length,
            base_outer_diameter=diameter,
            top_outer_diameter=diameter,
            base_thickness=thickness,
            top_thickness=thickness,
        )
        coordinates = np.array(
            [
                [0.0, 0.0],
                [length, 0.0],
            ],
            dtype=np.float64,
        )

        result = compute_elastic_beam_element_stiffness(
            node_coordinates=coordinates,
            tower_axis_start=0.0,
            tower_axis_end=length,
            tower_geometry=constant_tower,
            elastic_modulus=self.elastic_modulus,
            n_gauss=4,
            n_circumferential=32,
            n_radial=2,
        )

        section = create_annular_fiber_section(
            outer_diameter=diameter,
            thickness=thickness,
            n_circumferential=32,
            n_radial=2,
        )
        axial_rigidity = self.elastic_modulus * section.area
        bending_rigidity = (
            self.elastic_modulus * section.second_moment_x
        )

        expected = np.zeros((6, 6), dtype=np.float64)

        axial = axial_rigidity / length
        expected[0, 0] = axial
        expected[0, 3] = -axial
        expected[3, 0] = -axial
        expected[3, 3] = axial

        coefficient = bending_rigidity / length ** 3
        expected[np.ix_([1, 2, 4, 5], [1, 2, 4, 5])] = (
            coefficient
            * np.array(
                [
                    [12.0, 6.0 * length, -12.0, 6.0 * length],
                    [
                        6.0 * length,
                        4.0 * length ** 2,
                        -6.0 * length,
                        2.0 * length ** 2,
                    ],
                    [
                        -12.0,
                        -6.0 * length,
                        12.0,
                        -6.0 * length,
                    ],
                    [
                        6.0 * length,
                        2.0 * length ** 2,
                        -6.0 * length,
                        4.0 * length ** 2,
                    ],
                ],
                dtype=np.float64,
            )
        )

        np.testing.assert_allclose(
            result.local_stiffness,
            expected,
            rtol=1.0e-12,
            atol=1.0e-4,
        )
        np.testing.assert_allclose(
            result.global_stiffness,
            expected,
            rtol=1.0e-12,
            atol=1.0e-4,
        )

    def test_tapered_element_uses_gauss_point_geometry(self) -> None:
        mesh = create_uniform_vertical_tower_mesh(
            height=self.tower.height,
            n_elements=40,
        )
        first_element = mesh.connectivity[0]
        node_coordinates = mesh.coordinates[first_element]

        result = compute_elastic_beam_element_stiffness(
            node_coordinates=node_coordinates,
            tower_axis_start=0.0,
            tower_axis_end=mesh.tower_axis_coordinates[1],
            tower_geometry=self.tower,
            elastic_modulus=self.elastic_modulus,
            n_gauss=4,
            n_circumferential=32,
            n_radial=2,
        )

        self.assertEqual(result.gauss_heights.shape, (4,))
        self.assertTrue(
            np.all(np.diff(result.gauss_heights) > 0.0)
        )
        self.assertTrue(
            np.all(np.diff(result.outer_diameters) < 0.0)
        )
        self.assertTrue(
            np.all(np.diff(result.thicknesses) < 0.0)
        )
        self.assertTrue(
            np.all(np.diff(result.section_areas) < 0.0)
        )
        self.assertTrue(
            np.all(
                np.diff(result.section_second_moments) < 0.0
            )
        )

        expected_areas = np.array(
            [
                self.tower.exact_area_at(height)
                for height in result.gauss_heights
            ],
            dtype=np.float64,
        )
        np.testing.assert_allclose(
            result.section_areas,
            expected_areas,
            rtol=0.0,
            atol=1.0e-14,
        )

        np.testing.assert_allclose(
            result.local_stiffness,
            result.local_stiffness.T,
            rtol=0.0,
            atol=1.0e-8,
        )
        np.testing.assert_allclose(
            result.global_stiffness,
            result.global_stiffness.T,
            rtol=0.0,
            atol=1.0e-8,
        )

    def test_invalid_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            LinearTaperedTowerGeometry(
                height=0.0,
                base_outer_diameter=6.0,
                top_outer_diameter=3.87,
                base_thickness=0.027,
                top_thickness=0.019,
            )

        with self.assertRaises(ValueError):
            create_uniform_vertical_tower_mesh(
                height=87.6,
                n_elements=0,
            )

        with self.assertRaises(ValueError):
            self.tower.outer_diameter_at(-1.0)

        with self.assertRaises(ValueError):
            euler_bernoulli_strain_displacement(
                length=2.0,
                natural_coordinate=1.1,
            )


if __name__ == "__main__":
    unittest.main()

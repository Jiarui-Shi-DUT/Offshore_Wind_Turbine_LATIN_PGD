# -*- coding: utf-8 -*-
"""Regression tests for elastic tower response recovery."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.fiber_section import create_annular_fiber_section
from fem.tower_response_2d import recover_elastic_tower_response
from fem.tower_system_2d import solve_elastic_tower_top_force


class TestTowerResponse2D(unittest.TestCase):
    """Verify displacement-to-section and fiber-response recovery."""

    def setUp(self) -> None:
        self.height = 10.0
        self.n_elements = 4
        self.outer_diameter = 6.0
        self.thickness = 0.027
        self.elastic_modulus = 210.0e9
        self.force = 1.0e6
        self.n_gauss = 4
        self.n_circumferential = 32
        self.n_radial = 2

        self.mesh = create_uniform_vertical_tower_mesh(
            height=self.height,
            n_elements=self.n_elements,
        )
        self.tower = LinearTaperedTowerGeometry(
            height=self.height,
            base_outer_diameter=self.outer_diameter,
            top_outer_diameter=self.outer_diameter,
            base_thickness=self.thickness,
            top_thickness=self.thickness,
        )
        self.assembly, self.solution = (
            solve_elastic_tower_top_force(
                mesh=self.mesh,
                tower_geometry=self.tower,
                elastic_modulus=self.elastic_modulus,
                horizontal_force=self.force,
                n_gauss=self.n_gauss,
                n_circumferential=self.n_circumferential,
                n_radial=self.n_radial,
            )
        )
        self.response = recover_elastic_tower_response(
            assembly=self.assembly,
            solution=self.solution,
            tower_geometry=self.tower,
            elastic_modulus=self.elastic_modulus,
            n_circumferential=self.n_circumferential,
            n_radial=self.n_radial,
        )

    def test_response_array_shapes(self) -> None:
        """Recovered arrays must follow the element/Gauss/fiber layout."""
        n_fibers = (
            self.n_circumferential * self.n_radial
        )

        self.assertEqual(
            self.response.local_displacements.shape,
            (self.n_elements, 6),
        )
        self.assertEqual(
            self.response.gauss_heights.shape,
            (self.n_elements, self.n_gauss),
        )
        self.assertEqual(
            self.response.curvatures.shape,
            (self.n_elements, self.n_gauss),
        )
        self.assertEqual(
            self.response.bending_moments.shape,
            (self.n_elements, self.n_gauss),
        )
        self.assertEqual(
            self.response.fiber_strains.shape,
            (
                self.n_elements,
                self.n_gauss,
                n_fibers,
            ),
        )
        self.assertEqual(
            self.response.fiber_stresses.shape,
            (
                self.n_elements,
                self.n_gauss,
                n_fibers,
            ),
        )
        self.assertEqual(
            self.response.n_elements,
            self.n_elements,
        )
        self.assertEqual(
            self.response.n_gauss,
            self.n_gauss,
        )
        self.assertEqual(
            self.response.n_fibers,
            n_fibers,
        )

    def test_gauss_heights_are_ascending_and_inside_tower(self) -> None:
        """Gauss-point heights must preserve the physical tower order."""
        heights = self.response.flattened_gauss_heights

        self.assertTrue(np.all(np.diff(heights) > 0.0))
        self.assertGreater(heights[0], 0.0)
        self.assertLess(heights[-1], self.height)

    def test_pure_bending_has_negligible_axial_response(self) -> None:
        """A horizontal tip load must not create axial section force."""
        np.testing.assert_allclose(
            self.response.axial_strains,
            0.0,
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            self.response.axial_forces,
            0.0,
            rtol=0.0,
            atol=1.0e-3,
        )

    def test_recovered_moment_matches_cantilever_solution(self) -> None:
        """
        Internal moment must match the analytical tip-force distribution.

        For the local sign convention used by the vertical beam,

            M(z) = -F * (H - z).
        """
        expected_moments = (
            -self.force
            * (
                self.height
                - self.response.gauss_heights
            )
        )

        np.testing.assert_allclose(
            self.response.bending_moments,
            expected_moments,
            rtol=1.0e-10,
            atol=1.0e-3,
        )

    def test_curvature_and_fiber_stress_are_consistent(self) -> None:
        """Recovered moment, curvature, strain, and stress must agree."""
        section = create_annular_fiber_section(
            outer_diameter=self.outer_diameter,
            thickness=self.thickness,
            n_circumferential=self.n_circumferential,
            n_radial=self.n_radial,
        )
        expected_moments = (
            self.elastic_modulus
            * section.second_moment_x
            * self.response.curvatures
        )

        np.testing.assert_allclose(
            self.response.bending_moments,
            expected_moments,
            rtol=1.0e-12,
            atol=1.0e-4,
        )
        np.testing.assert_allclose(
            self.response.fiber_stresses,
            self.elastic_modulus
            * self.response.fiber_strains,
            rtol=1.0e-14,
            atol=1.0e-6,
        )

    def test_opposite_section_sides_have_opposite_stress(self) -> None:
        """Pure bending must create tensile and compressive fiber stress."""
        minimum = self.response.minimum_fiber_stresses
        maximum = self.response.maximum_fiber_stresses

        self.assertTrue(np.all(minimum < 0.0))
        self.assertTrue(np.all(maximum > 0.0))
        np.testing.assert_allclose(
            maximum,
            -minimum,
            rtol=1.0e-12,
            atol=1.0e-5,
        )


if __name__ == "__main__":
    unittest.main()

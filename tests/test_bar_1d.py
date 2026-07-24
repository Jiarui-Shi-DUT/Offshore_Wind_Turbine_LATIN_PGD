# -*- coding: utf-8 -*-
"""Regression tests for the one-dimensional bar finite-element module."""

import unittest

import numpy as np

from fem.bar_1d import (
    create_uniform_bar_mesh,
    solve_displacement_controlled_bar,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestBar1D(unittest.TestCase):
    """Verify mesh construction and the elastic constant-strain patch test."""

    def test_uniform_mesh(self) -> None:
        mesh = create_uniform_bar_mesh(
            length=1000.0,
            n_elements=90,
        )

        self.assertEqual(mesh.n_nodes, 91)
        self.assertEqual(mesh.n_elements, 90)
        np.testing.assert_allclose(
            mesh.element_lengths,
            np.full(90, 1000.0 / 90.0),
            rtol=0.0,
            atol=1.0e-12,
        )

    def test_elastic_patch(self) -> None:
        length = 1000.0
        area = 100.0
        imposed_displacement = 0.1
        n_elements = 10

        mesh = create_uniform_bar_mesh(
            length=length,
            n_elements=n_elements,
        )

        elastic_material = MaterialParameters(
            sigma_y=1.0e12,
        )

        response = solve_displacement_controlled_bar(
            mesh=mesh,
            area=area,
            materials=[elastic_material] * n_elements,
            time=np.array([0.0, 1.0], dtype=np.float64),
            right_displacement=np.array(
                [0.0, imposed_displacement],
                dtype=np.float64,
            ),
        )

        expected_strain = imposed_displacement / length
        expected_stress = elastic_material.E * expected_strain
        expected_reaction = expected_stress * area
        expected_displacement = (
            mesh.coordinates * imposed_displacement / length
        )

        np.testing.assert_allclose(
            response.displacement[-1],
            expected_displacement,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            response.strain[-1],
            expected_strain,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            response.stress[-1],
            expected_stress,
            rtol=0.0,
            atol=1.0e-10,
        )
        self.assertAlmostEqual(
            response.reaction_right[-1],
            expected_reaction,
            places=9,
        )
        self.assertAlmostEqual(
            response.reaction_left[-1],
            -expected_reaction,
            places=9,
        )
        np.testing.assert_allclose(
            response.state[-1],
            0.0,
            rtol=0.0,
            atol=1.0e-14,
        )


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""Regression tests for the one-dimensional LATIN equilibrium operator."""

import unittest

import numpy as np

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
)
from fem.bar_1d import create_uniform_bar_mesh
from latin.equilibrium_operator import apply_equilibrium_operator


class TestLatinEquilibriumOperator(unittest.TestCase):
    """Verify compatibility, equilibrium and constitutive consistency."""

    def test_full_space_time_projection(self) -> None:
        configuration = BenchmarkConfiguration()

        mesh = create_uniform_bar_mesh(
            length=configuration.length,
            n_elements=configuration.n_elements,
        )
        materials = create_three_material_distribution(
            configuration.n_elements
        )
        time = create_time_grid(
            total_time=configuration.total_time,
            time_step=configuration.time_step,
        )

        element_centres = 0.5 * (
            mesh.coordinates[:-1] + mesh.coordinates[1:]
        )
        source_strain = (
            1.0e-3
            * np.sin(
                2.0 * np.pi * time / configuration.period
            )[:, np.newaxis]
            * (
                element_centres / configuration.length
            )[np.newaxis, :]
        )

        projection = apply_equilibrium_operator(
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            source_strain=source_strain,
        )

        self.assertEqual(
            projection.compatible_strain.shape,
            (2001, 90),
        )
        self.assertEqual(
            projection.stress.shape,
            (2001, 90),
        )
        self.assertEqual(
            projection.displacement.shape,
            (2001, 91),
        )

        np.testing.assert_allclose(
            projection.displacement[:, 0],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            projection.displacement[:, -1],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            np.ptp(projection.stress, axis=1),
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            projection.reaction_left
            + projection.reaction_right,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

        elastic_moduli = np.array(
            [material.E for material in materials],
            dtype=np.float64,
        )
        constitutive_stress = (
            elastic_moduli[np.newaxis, :]
            * (
                projection.compatible_strain
                - source_strain
            )
        )

        np.testing.assert_allclose(
            projection.stress,
            constitutive_stress,
            rtol=0.0,
            atol=2.0e-14,
        )

        extension = (
            projection.compatible_strain
            @ mesh.element_lengths
        )
        np.testing.assert_allclose(
            extension,
            0.0,
            rtol=0.0,
            atol=1.0e-12,
        )


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""Regression tests for the LATIN elastic initialization."""

import unittest

import numpy as np

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
    prescribed_displacement,
)
from fem.bar_1d import create_uniform_bar_mesh
from latin.initialization import compute_elastic_initialization


class TestLatinElasticInitialization(unittest.TestCase):
    """Verify the 90-element full-time elastic initialization."""

    def test_full_three_material_bar_initialization(self) -> None:
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
        right_displacement = prescribed_displacement(
            time=time,
            amplitude=configuration.displacement_amplitude,
            period=configuration.period,
        )

        initialization = compute_elastic_initialization(
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            time=time,
            right_displacement=right_displacement,
        )

        self.assertEqual(
            initialization.displacement.shape,
            (2001, 91),
        )
        self.assertEqual(
            initialization.strain.shape,
            (2001, 90),
        )
        self.assertEqual(
            initialization.stress.shape,
            (2001, 90),
        )
        self.assertEqual(
            initialization.state.shape,
            (2001, 90, 4),
        )

        np.testing.assert_allclose(
            initialization.displacement[:, 0],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            initialization.displacement[:, -1],
            right_displacement,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            np.ptp(initialization.stress, axis=1),
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            initialization.state,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            initialization.reaction_left
            + initialization.reaction_right,
            0.0,
            rtol=0.0,
            atol=0.0,
        )


if __name__ == "__main__":
    unittest.main()

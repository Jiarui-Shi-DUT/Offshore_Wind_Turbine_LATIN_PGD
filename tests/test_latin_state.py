# -*- coding: utf-8 -*-
"""Regression tests for the LATIN space-time state container."""

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
from latin.state import LatinState


class TestLatinState(unittest.TestCase):
    """Verify construction of the initial LATIN space-time state."""

    def test_from_full_elastic_initialization(self) -> None:
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

        elastic = compute_elastic_initialization(
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            time=time,
            right_displacement=right_displacement,
        )
        state = LatinState.from_elastic_initialization(
            initialization=elastic,
            materials=materials,
        )

        self.assertEqual(state.field_shape, (2001, 90))
        self.assertEqual(state.n_time, 2001)
        self.assertEqual(state.n_elements, 90)

        np.testing.assert_allclose(
            state.stress,
            elastic.stress,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            state.elastic_strain,
            elastic.strain,
            rtol=0.0,
            atol=0.0,
        )

        zero_fields = (
            state.plastic_strain_rate,
            state.alpha_rate,
            state.r_bar_rate,
            state.damage_rate,
            state.plastic_strain,
            state.alpha,
            state.r_bar,
            state.damage,
            state.beta,
            state.R_bar,
        )
        for field in zero_fields:
            np.testing.assert_allclose(
                field,
                0.0,
                rtol=0.0,
                atol=0.0,
            )

        self.assertGreaterEqual(
            float(np.min(state.energy_release_rate)),
            0.0,
        )
        self.assertAlmostEqual(
            float(np.max(state.energy_release_rate)),
            0.09648,
            places=12,
        )

    def test_copy_is_independent(self) -> None:
        time = np.array([0.0, 0.1], dtype=np.float64)
        state = LatinState.zeros(time=time, n_elements=2)
        copied = state.copy()

        copied.stress[1, 0] = 10.0
        copied.damage[1, 1] = 0.2

        self.assertEqual(float(state.stress[1, 0]), 0.0)
        self.assertEqual(float(state.damage[1, 1]), 0.0)


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""Regression tests for LATIN descent search directions."""

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
from latin.local_stage import solve_local_stage
from latin.search_directions import (
    compute_descent_search_directions,
)
from latin.state import LatinState


class TestLatinSearchDirections(unittest.TestCase):
    """Verify the regularized descent operators on the full benchmark."""

    def test_full_three_material_bar_search_directions(self) -> None:
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
        global_state = LatinState.from_elastic_initialization(
            initialization=elastic,
            materials=materials,
        )
        local_state = solve_local_stage(
            global_state=global_state,
            materials=materials,
        )

        directions = compute_descent_search_directions(
            local_state=local_state,
            materials=materials,
            regularization=0.15,
        )

        self.assertEqual(directions.field_shape, (2001, 90))
        self.assertAlmostEqual(
            directions.regularization,
            0.15,
            places=15,
        )

        for field in (
            directions.H_sigma,
            directions.H_beta,
            directions.H_R_bar,
        ):
            self.assertTrue(np.all(np.isfinite(field)))
            self.assertTrue(np.all(field > 0.0))

        np.testing.assert_allclose(
            directions.b_damage,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

        self.assertAlmostEqual(
            float(np.min(directions.H_sigma)),
            1.1194029850746267e-06,
            places=18,
        )
        self.assertAlmostEqual(
            float(np.max(directions.H_sigma)),
            5.315952958867752e-04,
            places=16,
        )
        self.assertAlmostEqual(
            float(np.min(directions.H_beta)),
            2.7272727272727273e-05,
            places=18,
        )
        self.assertAlmostEqual(
            float(np.max(directions.H_beta)),
            7.943011909254297e-04,
            places=16,
        )
        self.assertAlmostEqual(
            float(np.min(directions.H_R_bar)),
            5.0e-03,
            places=15,
        )
        self.assertAlmostEqual(
            float(np.max(directions.H_R_bar)),
            5.758336731217285e-03,
            places=15,
        )


if __name__ == "__main__":
    unittest.main()

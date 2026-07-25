# -*- coding: utf-8 -*-
"""Regression test for the full-order LATIN outer iteration driver."""

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
from latin.solver import solve_full_order_latin
from latin.state import LatinState


class TestFullOrderLatinSolver(unittest.TestCase):
    """Verify iteration, relaxation history and convergence stopping."""

    def test_two_iteration_convergence_control(self) -> None:
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
        initial_state = LatinState.from_elastic_initialization(
            initialization=elastic,
            materials=materials,
        )

        result = solve_full_order_latin(
            initial_state=initial_state,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            tolerance=0.15,
            max_iterations=2,
            relaxation=0.8,
        )

        self.assertTrue(result.converged)
        self.assertEqual(result.iterations, 2)
        self.assertEqual(result.indicator_history.shape, (2,))
        self.assertTrue(np.all(np.isfinite(result.indicator_history)))
        self.assertTrue(np.all(np.isfinite(result.state.stress)))
        self.assertTrue(np.all(np.isfinite(result.state.damage)))

        self.assertAlmostEqual(
            result.indicator_history[0],
            0.17831837606625445,
            places=12,
        )
        self.assertGreater(result.indicator_history[0], 0.15)
        self.assertLessEqual(result.indicator_history[1], 0.15)
        self.assertGreater(result.indicator_history[1], 0.0)
        self.assertLess(
            result.indicator_history[1],
            result.indicator_history[0],
        )
        self.assertEqual(
            result.final_indicator,
            float(result.indicator_history[-1]),
        )


if __name__ == "__main__":
    unittest.main()

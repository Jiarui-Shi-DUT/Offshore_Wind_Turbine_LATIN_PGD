# -*- coding: utf-8 -*-
"""Regression tests for LATIN relaxation and convergence control."""

import unittest

import numpy as np

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
    prescribed_displacement,
)
from fem.bar_1d import create_uniform_bar_mesh
from latin.global_stage import solve_full_order_global_stage
from latin.initialization import compute_elastic_initialization
from latin.iteration_control import (
    relative_latin_indicator,
    relax_global_state,
)
from latin.local_stage import solve_local_stage
from latin.search_directions import (
    DescentSearchDirections,
    compute_descent_search_directions,
)
from latin.state import LatinState


class TestLatinIterationControl(unittest.TestCase):
    """Verify the first relaxed LATIN update and relative indicator."""

    def test_first_full_relaxed_iteration(self) -> None:
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
        global_state_0 = LatinState.from_elastic_initialization(
            initialization=elastic,
            materials=materials,
        )
        local_state = solve_local_stage(
            global_state=global_state_0,
            materials=materials,
        )
        directions = compute_descent_search_directions(
            local_state=local_state,
            materials=materials,
        )
        candidate_state = solve_full_order_global_stage(
            global_state=global_state_0,
            local_state=local_state,
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
        ).state

        global_state_1 = relax_global_state(
            previous_state=global_state_0,
            candidate_state=candidate_state,
            relaxation=0.8,
        )
        indicator = relative_latin_indicator(
            local_state=local_state,
            global_state=global_state_1,
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
        )

        expected_stress = (
            0.2 * global_state_0.stress
            + 0.8 * candidate_state.stress
        )
        expected_damage = (
            0.2 * global_state_0.damage
            + 0.8 * candidate_state.damage
        )

        np.testing.assert_allclose(
            global_state_1.stress,
            expected_stress,
            rtol=0.0,
            atol=3.0e-14,
        )
        np.testing.assert_allclose(
            global_state_1.damage,
            expected_damage,
            rtol=0.0,
            atol=1.0e-15,
        )

        self.assertTrue(np.all(np.isfinite(global_state_1.stress)))
        self.assertTrue(np.all(np.isfinite(global_state_1.damage)))
        self.assertTrue(np.isfinite(indicator))
        self.assertGreaterEqual(indicator, 0.0)
        self.assertLessEqual(indicator, 1.0)
        self.assertAlmostEqual(
            indicator,
            0.17831837606625445,
            places=12,
        )

    def test_identical_states_have_zero_indicator(self) -> None:
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
        state = LatinState.zeros(
            time=time,
            n_elements=configuration.n_elements,
        )

        shape = state.field_shape
        directions = DescentSearchDirections(
            H_sigma=np.ones(shape, dtype=np.float64),
            H_beta=np.ones(shape, dtype=np.float64),
            H_R_bar=np.ones(shape, dtype=np.float64),
            b_damage=np.zeros(shape, dtype=np.float64),
            regularization=0.15,
        )

        indicator = relative_latin_indicator(
            local_state=state,
            global_state=state.copy(),
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
        )

        self.assertEqual(indicator, 0.0)


if __name__ == "__main__":
    unittest.main()

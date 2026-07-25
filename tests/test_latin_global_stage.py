# -*- coding: utf-8 -*-
"""Regression tests for the full-order LATIN global stage."""

import unittest

import numpy as np

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
    material_region_slices,
    prescribed_displacement,
)
from fem.bar_1d import create_uniform_bar_mesh
from latin.global_stage import solve_full_order_global_stage
from latin.initialization import compute_elastic_initialization
from latin.local_stage import solve_local_stage
from latin.search_directions import compute_descent_search_directions
from latin.state import LatinState


class TestLatinGlobalStage(unittest.TestCase):
    """Verify the first full local-global LATIN iteration."""

    def test_first_full_order_global_stage(self) -> None:
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

        result = solve_full_order_global_stage(
            global_state=global_state_0,
            local_state=local_state,
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
        )
        global_state_1 = result.state

        self.assertEqual(global_state_1.field_shape, (2001, 90))

        for field in (
            global_state_1.stress,
            global_state_1.elastic_strain,
            global_state_1.plastic_strain,
            global_state_1.damage,
            global_state_1.beta,
            global_state_1.R_bar,
            global_state_1.energy_release_rate,
            result.plastic_strain_correction,
            result.plastic_strain_rate_correction,
            result.residual_strain,
            result.displacement_correction,
        ):
            self.assertTrue(np.all(np.isfinite(field)))

        np.testing.assert_allclose(
            result.displacement_correction[:, 0],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            result.displacement_correction[:, -1],
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            np.ptp(global_state_1.stress, axis=1),
            0.0,
            rtol=0.0,
            atol=0.0,
        )

        regions = material_region_slices(
            configuration.n_elements
        )
        final_regional_damage = np.array(
            [
                np.mean(global_state_1.damage[-1, region])
                for region in regions
            ],
            dtype=np.float64,
        )
        expected_damage = np.array(
            [
                0.4608126165335944,
                0.43635333604334353,
                0.41202335421488834,
            ],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            final_regional_damage,
            expected_damage,
            rtol=1.0e-10,
            atol=1.0e-12,
        )

        self.assertAlmostEqual(
            float(np.min(global_state_1.stress)),
            -132.05030255534186,
            places=10,
        )
        self.assertAlmostEqual(
            float(np.max(global_state_1.stress)),
            132.36848647392353,
            places=10,
        )
        self.assertAlmostEqual(
            float(np.max(np.abs(
                result.plastic_strain_correction
            ))),
            0.0010476867048401905,
            places=15,
        )
        self.assertAlmostEqual(
            float(np.max(np.abs(result.residual_strain))),
            0.0009801736163273973,
            places=15,
        )


if __name__ == "__main__":
    unittest.main()

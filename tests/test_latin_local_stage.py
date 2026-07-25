# -*- coding: utf-8 -*-
"""Regression tests for the nonlinear LATIN local stage."""

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
from latin.initialization import compute_elastic_initialization
from latin.local_stage import (
    solve_local_stage,
    unilateral_elastic_strain,
)
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


class TestLatinLocalStage(unittest.TestCase):
    """Verify the first full local projection on the paper benchmark."""

    def test_first_full_local_projection(self) -> None:
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

        self.assertEqual(local_state.field_shape, (2001, 90))

        for field in (
            local_state.plastic_strain,
            local_state.alpha,
            local_state.r_bar,
            local_state.damage,
            local_state.plastic_strain_rate,
            local_state.alpha_rate,
            local_state.r_bar_rate,
            local_state.damage_rate,
            local_state.elastic_strain,
        ):
            self.assertTrue(np.all(np.isfinite(field)))

        # Infinite ascent directions keep thermodynamic forces unchanged.
        np.testing.assert_allclose(
            local_state.stress,
            global_state.stress,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            local_state.beta,
            global_state.beta,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            local_state.R_bar,
            global_state.R_bar,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            local_state.energy_release_rate,
            global_state.energy_release_rate,
            rtol=0.0,
            atol=0.0,
        )

        self.assertGreater(float(np.max(local_state.damage)), 0.0)
        self.assertTrue(np.all(local_state.damage >= 0.0))
        self.assertTrue(np.all(local_state.damage < 1.0))

        regions = material_region_slices(
            configuration.n_elements
        )
        final_regional_damage = np.array(
            [
                np.mean(local_state.damage[-1, region])
                for region in regions
            ],
            dtype=np.float64,
        )
        expected = np.array(
            [
                0.46012042913511303,
                0.43568546963097937,
                0.4113731130752918,
            ],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            final_regional_damage,
            expected,
            rtol=1.0e-10,
            atol=1.0e-12,
        )

        self.assertGreater(
            final_regional_damage[0],
            final_regional_damage[1],
        )
        self.assertGreater(
            final_regional_damage[1],
            final_regional_damage[2],
        )

    def test_unilateral_elastic_strain(self) -> None:
        material = MaterialParameters()
        damage = 0.25
        stress = 100.0

        tensile = unilateral_elastic_strain(
            stress=stress,
            damage=damage,
            material=material,
        )
        compressive = unilateral_elastic_strain(
            stress=-stress,
            damage=damage,
            material=material,
        )

        expected_tensile = stress / (
            material.E * (1.0 - damage)
        )
        expected_compressive = -stress / (
            material.E * (1.0 - material.h * damage)
        )

        self.assertAlmostEqual(
            tensile,
            expected_tensile,
            places=15,
        )
        self.assertAlmostEqual(
            compressive,
            expected_compressive,
            places=15,
        )


if __name__ == "__main__":
    unittest.main()

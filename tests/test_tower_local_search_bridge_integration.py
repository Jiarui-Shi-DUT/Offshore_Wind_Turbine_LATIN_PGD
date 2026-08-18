# -*- coding: utf-8 -*-
"""Actual tower-q integration test for the local/search bridge."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_local_stage import (
    solve_tower_local_stage,
)
from latin.tower_search_directions import (
    compute_tower_descent_search_directions,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerLocalSearchBridgeIntegration(unittest.TestCase):
    """Run the bridge on the actual beam-Gauss-fiber material-point count."""

    def test_real_tower_material_point_grid(self) -> None:
        geometry = LinearTaperedTowerGeometry(
            height=12.0,
            base_outer_diameter=2.0,
            top_outer_diameter=1.6,
            base_thickness=0.050,
            top_thickness=0.035,
        )
        mesh = create_uniform_vertical_tower_mesh(
            height=geometry.height,
            n_elements=2,
        )
        material = MaterialParameters(E=210_000.0)
        system = ViscoplasticDamageTowerSystem2D(
            mesh=mesh,
            tower_geometry=geometry,
            material=material,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        operator = build_tower_equilibrium_operator(system)

        time = np.array(
            [0.0, 0.05, 0.12, 0.20],
            dtype=np.float64,
        )
        q = np.arange(
            operator.n_material_points,
            dtype=np.float64,
        )[None, :]
        t = np.arange(
            time.size,
            dtype=np.float64,
        )[:, None]

        stress = (
            95.0
            + 12.0 * t
            + 0.15 * q
        )
        beta = (
            1.0
            + 0.05 * t
            + 0.002 * q
        )
        R_bar = (
            1.5
            + 0.08 * t
            + 0.003 * q
        )
        energy = (
            0.026
            + 0.003 * t
            + 0.00002 * q
        )
        shape = stress.shape
        zeros = np.zeros(shape, dtype=np.float64)

        global_state = LatinStateTower(
            time=time,
            plastic_strain_rate=zeros,
            elastic_strain=zeros,
            alpha_rate=zeros,
            r_bar_rate=zeros,
            damage_rate=zeros,
            stress=stress,
            beta=beta,
            R_bar=R_bar,
            energy_release_rate=energy,
            plastic_strain=zeros,
            alpha=zeros,
            r_bar=zeros,
            damage=zeros,
        )

        local_state = solve_tower_local_stage(
            global_state=global_state,
            materials=material,
        )
        directions = compute_tower_descent_search_directions(
            local_state=local_state,
            materials=material,
            regularization=0.15,
        )

        self.assertEqual(
            local_state.n_material_points,
            operator.n_material_points,
        )
        self.assertEqual(
            directions.field_shape,
            global_state.field_shape,
        )
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

        self.assertTrue(
            np.all(np.isfinite(local_state.plastic_strain_rate))
        )
        self.assertTrue(
            np.all(np.isfinite(local_state.damage))
        )
        self.assertTrue(
            np.all(local_state.damage >= 0.0)
        )
        self.assertTrue(
            np.all(local_state.damage < 1.0)
        )
        self.assertTrue(np.all(directions.H_sigma > 0.0))
        self.assertTrue(np.all(directions.H_beta > 0.0))
        self.assertTrue(np.all(directions.H_R_bar > 0.0))
        self.assertTrue(np.all(directions.b_damage == 0.0))


if __name__ == "__main__":
    unittest.main()

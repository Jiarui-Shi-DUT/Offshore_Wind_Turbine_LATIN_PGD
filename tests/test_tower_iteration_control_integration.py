# -*- coding: utf-8 -*-
"""Actual fiber-beam integration test for tower LATIN iteration control."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_iteration_control import (
    evaluate_tower_trial,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


def _tower_state(time, base, *, scale, damage):
    shape = base.shape
    return LatinStateTower(
        time=time,
        plastic_strain_rate=scale * 0.01 * base,
        elastic_strain=scale * 0.001 * base,
        alpha_rate=scale * 0.02 * base,
        r_bar_rate=scale * 0.03 * base,
        damage_rate=scale * 0.001 * base,
        stress=scale * 10.0 * base,
        beta=scale * 1.5 * base,
        R_bar=scale * 2.0 * base,
        energy_release_rate=scale * 0.4 * base,
        plastic_strain=scale * 0.05 * base,
        alpha=scale * 0.06 * base,
        r_bar=scale * 0.07 * base,
        damage=np.full(shape, damage, dtype=np.float64),
    )


class TestTowerIterationControlIntegration(unittest.TestCase):
    """Exercise Eqs. (60), (76)-(77) on the real tower q-metric."""

    def test_actual_tower_metric_scaling_case(self) -> None:
        geometry = LinearTaperedTowerGeometry(
            height=10.0,
            base_outer_diameter=1.8,
            top_outer_diameter=1.5,
            base_thickness=0.045,
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
        metric = operator.metric

        time = np.array(
            [0.0, 0.20, 0.55, 1.00],
            dtype=np.float64,
        )
        q = np.arange(
            operator.n_material_points,
            dtype=np.float64,
        )
        base = (
            1.0
            + 0.15 * np.arange(time.size)[:, None]
            + 0.01 * q[None, :]
        )
        shape = base.shape
        zero = np.zeros(shape, dtype=np.float64)

        baseline = LatinStateTower(
            time=time,
            plastic_strain_rate=zero,
            elastic_strain=zero,
            alpha_rate=zero,
            r_bar_rate=zero,
            damage_rate=zero,
            stress=zero,
            beta=zero,
            R_bar=zero,
            energy_release_rate=zero,
            plastic_strain=zero,
            alpha=zero,
            r_bar=zero,
            damage=zero,
        )
        local = _tower_state(
            time,
            base,
            scale=1.0,
            damage=0.10,
        )
        unrelaxed = local.copy()

        H_sigma = (
            2.0e-6
            + 1.0e-7 * np.arange(time.size)[:, None]
            + 1.0e-9 * q[None, :]
        )
        H_beta = (
            3.0e-4
            + 2.0e-5 * np.arange(time.size)[:, None]
            + 1.0e-7 * q[None, :]
        )
        H_R_bar = (
            4.0e-3
            + 3.0e-4 * np.arange(time.size)[:, None]
            + 1.0e-6 * q[None, :]
        )
        directions = DescentSearchDirections(
            H_sigma=H_sigma,
            H_beta=H_beta,
            H_R_bar=H_R_bar,
            b_damage=np.zeros(shape, dtype=np.float64),
            regularization=0.15,
        )

        evaluation = evaluate_tower_trial(
            baseline_state=baseline,
            local_state=local,
            unrelaxed_state=unrelaxed,
            directions=directions,
            metric=metric,
            reference_modulus=operator.reference_modulus,
            previous_indicator=1.0,
            relaxation=0.8,
            latin_tolerance=0.0,
        )

        self.assertGreater(metric.total_measure, 0.0)
        self.assertEqual(
            metric.n_material_points,
            operator.n_material_points,
        )
        self.assertAlmostEqual(
            evaluation.indicator,
            1.0 / 9.0,
            places=11,
        )
        self.assertAlmostEqual(
            evaluation.saturation,
            0.8,
            places=11,
        )
        self.assertFalse(evaluation.converged)
        self.assertTrue(evaluation.finite)
        np.testing.assert_allclose(
            evaluation.relaxed_state.stress,
            0.8 * local.stress,
            rtol=0.0,
            atol=1.0e-14,
        )


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""End-to-end coarse fiber-beam test for the tower LATIN-PGD solver."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.viscoplastic_tower_system_2d import ViscoplasticDamageTowerSystem2D
from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import build_tower_equilibrium_operator
from latin.tower_latin_pgd_solver import (
    TowerLatinPGDTerminationReason,
    solve_tower_latin_pgd,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerLatinPGDSolverIntegration(unittest.TestCase):
    """Wire the complete outer transaction on the frozen 10 x 2 x 16 mesh."""

    def test_zero_equilibrium_history_converges_with_valid_empty_basis(self) -> None:
        geometry = LinearTaperedTowerGeometry(
            height=20.0,
            base_outer_diameter=2.0,
            top_outer_diameter=1.5,
            base_thickness=0.05,
            top_thickness=0.035,
        )
        mesh = create_uniform_vertical_tower_mesh(
            height=geometry.height,
            n_elements=10,
        )
        material = MaterialParameters(
            E=210_000.0,
            C=5_500.0,
            R_inf=30.0,
            h=0.2,
        )
        system = ViscoplasticDamageTowerSystem2D(
            mesh=mesh,
            tower_geometry=geometry,
            material=material,
            n_gauss=2,
            n_circumferential=16,
            n_radial=1,
        )
        operator = build_tower_equilibrium_operator(system)
        self.assertEqual(operator.n_material_points, 10 * 2 * 16)

        time = np.array([0.0, 0.25, 0.60, 1.0], dtype=np.float64)
        initial_state = LatinStateTower.zeros(
            time=time,
            n_material_points=operator.n_material_points,
        )

        result = solve_tower_latin_pgd(
            initial_state=initial_state,
            materials=material,
            metric=operator.metric,
            equilibrium_operator=operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            max_iterations=3,
        )

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.CONVERGED,
        )
        self.assertTrue(result.converged)
        self.assertEqual(result.iterations, 1)
        self.assertEqual(result.trial_evaluations, 1)
        self.assertEqual(result.trial_kind_history, ("A",))
        self.assertEqual(result.commit_kind_history, ("A",))
        self.assertEqual(result.basis.n_modes, 0)
        self.assertAlmostEqual(result.final_indicator, 0.0)
        self.assertEqual(result.total_modes_added, 0)
        self.assertIsNone(result.last_enrichment_result)

        for field_name in LatinStateTower.MATERIAL_FIELD_NAMES:
            np.testing.assert_allclose(
                getattr(result.state, field_name),
                np.zeros(initial_state.field_shape),
                rtol=0.0,
                atol=0.0,
            )

    def test_zero_history_with_one_existing_pair_keeps_transaction_consistent(self) -> None:
        geometry = LinearTaperedTowerGeometry(
            height=20.0,
            base_outer_diameter=2.0,
            top_outer_diameter=1.5,
            base_thickness=0.05,
            top_thickness=0.035,
        )
        mesh = create_uniform_vertical_tower_mesh(
            height=geometry.height,
            n_elements=10,
        )
        material = MaterialParameters(
            E=210_000.0,
            C=5_500.0,
            R_inf=30.0,
            h=0.2,
        )
        system = ViscoplasticDamageTowerSystem2D(
            mesh=mesh,
            tower_geometry=geometry,
            material=material,
            n_gauss=2,
            n_circumferential=16,
            n_radial=1,
        )
        operator = build_tower_equilibrium_operator(system)

        time = np.array([0.0, 0.25, 0.60, 1.0], dtype=np.float64)
        initial_state = LatinStateTower.zeros(
            time=time,
            n_material_points=operator.n_material_points,
        )
        q = np.arange(operator.n_material_points, dtype=np.float64)
        p = np.sin(0.071 * q) + 0.35 * np.cos(0.113 * q)
        p = p / operator.metric.norm(p)
        s = operator.apply_spatial(p).stress
        initial_basis = PGDBasisTower(
            n_material_points=operator.n_material_points,
            n_time=time.size,
            modes=(
                PGDModeTower(
                    spatial_plastic_strain=p,
                    spatial_stress=s,
                    temporal_amplitude=np.zeros(time.size),
                    temporal_rate=np.zeros(time.size),
                    iteration_added=0,
                ),
            ),
        )

        result = solve_tower_latin_pgd(
            initial_state=initial_state,
            materials=material,
            metric=operator.metric,
            equilibrium_operator=operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            initial_basis=initial_basis,
            max_iterations=3,
        )

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.CONVERGED,
        )
        self.assertEqual(result.basis.n_modes, 1)
        self.assertEqual(result.trial_kind_history, ("A",))
        self.assertEqual(result.commit_kind_history, ("A",))
        self.assertEqual(result.total_modes_added, 0)
        self.assertAlmostEqual(result.final_indicator, 0.0)
        np.testing.assert_allclose(
            result.basis.spatial_plastic_strain_matrix()[:, 0],
            initial_basis.spatial_plastic_strain_matrix()[:, 0],
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            result.basis.temporal_amplitude_matrix(),
            np.zeros((time.size, 1)),
            rtol=0.0,
            atol=0.0,
        )



if __name__ == "__main__":
    unittest.main()

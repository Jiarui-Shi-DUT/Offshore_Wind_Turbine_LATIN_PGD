# -*- coding: utf-8 -*-
"""Actual fiber-beam integration test for tower global-stage finishing."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.viscoplastic_tower_system_2d import ViscoplasticDamageTowerSystem2D
from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import build_tower_equilibrium_operator
from latin.tower_global_finishing import (
    build_unrelaxed_candidate,
    prepare_frozen_global_data,
)
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerGlobalFinishingIntegration(unittest.TestCase):
    """Close one reduced global candidate on the actual tower q-space."""

    def test_complete_candidate_on_actual_fiber_beam_tower(self) -> None:
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
            n_circumferential=8,
            n_radial=1,
        )
        operator = build_tower_equilibrium_operator(system)
        metric = operator.metric
        nq = operator.n_material_points

        time = np.array([0.0, 0.20, 0.45, 0.80, 1.20], dtype=np.float64)
        nt = time.size
        q = np.arange(nq, dtype=np.float64)[None, :]
        t = np.arange(nt, dtype=np.float64)[:, None]

        base_stress = 35.0 + 2.0 * t + 0.1 * np.sin(0.17 * q)
        base_elastic = base_stress / material.E
        zero = np.zeros((nt, nq), dtype=np.float64)
        baseline = LatinStateTower(
            time=time,
            plastic_strain_rate=zero,
            elastic_strain=base_elastic,
            alpha_rate=zero,
            r_bar_rate=zero,
            damage_rate=zero,
            stress=base_stress,
            beta=zero,
            R_bar=zero,
            energy_release_rate=zero,
            plastic_strain=zero,
            alpha=zero,
            r_bar=zero,
            damage=zero,
        )

        local_plastic_rate = 1.0e-5 * t + 1.0e-7 * q
        local_alpha_rate = 2.0e-5 + 1.0e-6 * t + 1.0e-8 * q
        local_r_rate = 1.5e-5 + 8.0e-7 * t + 1.0e-8 * q
        local_damage = 2.0e-4 * t + 1.0e-6 * q
        local_damage_rate = 2.0e-4 + 1.0e-5 * t + 1.0e-7 * q
        local = LatinStateTower(
            time=time,
            plastic_strain_rate=local_plastic_rate,
            elastic_strain=base_elastic * (1.0 + 0.03 * t),
            alpha_rate=local_alpha_rate,
            r_bar_rate=local_r_rate,
            damage_rate=local_damage_rate,
            stress=base_stress,
            beta=zero,
            R_bar=zero,
            energy_release_rate=zero,
            plastic_strain=2.0e-6 * t + 1.0e-8 * q,
            alpha=1.0e-6 * t + np.zeros_like(q),
            r_bar=8.0e-7 * t + np.zeros_like(q),
            damage=local_damage,
        )

        H_sigma = 2.0e-6 + 1.0e-7 * t + 1.0e-9 * q
        directions = DescentSearchDirections(
            H_sigma=H_sigma,
            H_beta=np.full((nt, nq), 2.0e-4, dtype=np.float64),
            H_R_bar=np.full((nt, nq), 3.0e-4, dtype=np.float64),
            b_damage=np.zeros((nt, nq), dtype=np.float64),
            regularization=0.15,
        )

        frozen = prepare_frozen_global_data(
            baseline,
            local,
            directions,
            operator,
        )

        p = 0.7 * np.sin(0.31 * q.ravel()) + 0.2 * np.cos(0.17 * q.ravel())
        p = p / metric.norm(p)
        s = operator.apply_spatial(p).stress
        basis = PGDBasisTower(
            n_material_points=nq,
            n_time=nt,
            modes=(
                PGDModeTower(
                    p,
                    s,
                    np.zeros(nt),
                    np.zeros(nt),
                ),
            ),
        )
        fixed = update_tower_pgd_time_functions(
            basis=basis,
            time=time,
            forcing=frozen.full_plastic_forcing,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
        )

        candidate = build_unrelaxed_candidate(
            baseline_state=baseline,
            local_state=local,
            directions=directions,
            frozen_data=frozen,
            fixed_basis_result=fixed,
            materials=material,
        )

        self.assertEqual(candidate.state.field_shape, (nt, nq))
        np.testing.assert_allclose(
            operator.equilibrium_residual(
                fixed.plastic_projection.stress
                + frozen.damage_projection.stress
            ),
            np.zeros((nt, operator.n_free_dofs)),
            rtol=0.0,
            atol=1.0e-7,
        )
        np.testing.assert_array_equal(candidate.state.damage, local.damage)
        np.testing.assert_array_equal(
            candidate.state.damage_rate,
            local.damage_rate,
        )
        np.testing.assert_allclose(
            candidate.state.beta,
            material.C * candidate.state.alpha,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            candidate.state.R_bar,
            material.R_inf * candidate.state.r_bar,
            rtol=0.0,
            atol=1.0e-12,
        )
        self.assertTrue(
            np.all(np.isfinite(candidate.state.energy_release_rate))
        )
        np.testing.assert_allclose(
            candidate.total_displacement_correction,
            fixed.plastic_projection.displacement_free
            + frozen.damage_projection.displacement_free,
            rtol=0.0,
            atol=1.0e-12,
        )


if __name__ == "__main__":
    unittest.main()

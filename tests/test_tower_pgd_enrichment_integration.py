# -*- coding: utf-8 -*-
"""Actual fiber-beam integration test for tower Eq. (61)-(72) enrichment."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.pgd_basis import PGDBasisTower
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_pgd_enrichment import (
    enrich_tower_pgd_basis_once,
)
from latin.tower_pgd_time_update import (
    update_tower_pgd_time_functions,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerPGDEnrichmentIntegration(unittest.TestCase):
    """Recover one residual-driven mode on the real tower material-point space."""

    def test_exact_rank_one_forcing_on_actual_tower(self) -> None:
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
            [0.0, 0.20, 0.45, 0.80, 1.20],
            dtype=np.float64,
        )
        q = np.arange(
            operator.n_material_points,
            dtype=np.float64,
        )
        p = (
            0.7 * np.sin(0.31 * q)
            + 0.2 * np.cos(0.17 * q)
        )
        p = p / metric.norm(p)
        s = operator.apply_spatial(p).stress

        amplitude = np.array(
            [0.0, 0.08, 0.20, 0.14, 0.30],
            dtype=np.float64,
        )
        rate = np.zeros_like(amplitude)
        rate[0] = (
            amplitude[1] - amplitude[0]
        ) / (
            time[1] - time[0]
        )
        rate[1:] = (
            np.diff(amplitude)
            / np.diff(time)
        )

        H_sigma = (
            2.0e-6
            + 1.0e-7
            * np.arange(time.size)[:, None]
            + 1.0e-9
            * q[None, :]
        )
        forcing = (
            np.outer(rate, p)
            - H_sigma * np.outer(amplitude, s)
        )

        empty_basis = PGDBasisTower(
            n_material_points=operator.n_material_points,
            n_time=time.size,
        )
        trial_a = update_tower_pgd_time_functions(
            basis=empty_basis,
            time=time,
            forcing=forcing,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
        )

        enrichment = enrich_tower_pgd_basis_once(
            fixed_basis_result=trial_a,
            time=time,
            full_forcing=forcing,
            shifted_defect=trial_a.mechanical_residual,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
            fixed_point_tolerance=1.0e-7,
            max_fixed_point_iterations=100,
        )

        self.assertTrue(
            enrichment.accepted,
            msg=enrichment.failure_reason,
        )
        self.assertEqual(enrichment.n_modes, 1)
        self.assertTrue(enrichment.fixed_point_converged)
        self.assertGreater(
            enrichment.residual_benefit,
            0.99,
        )
        self.assertLess(
            enrichment.orthogonality_error,
            1.0e-10,
        )
        np.testing.assert_allclose(
            operator.equilibrium_residual(
                enrichment.candidate_fixed_basis_result
                .plastic_projection.stress
            ),
            np.zeros(
                (
                    time.size,
                    operator.n_free_dofs,
                )
            ),
            rtol=0.0,
            atol=1.0e-7,
        )


if __name__ == "__main__":
    unittest.main()

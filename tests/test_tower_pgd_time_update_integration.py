# -*- coding: utf-8 -*-
"""Integration test for tower Eq. (58)-(59) on the actual fiber-beam FEM."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_pgd_time_update import (
    update_tower_pgd_time_functions,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerPGDTimeUpdateIntegration(unittest.TestCase):
    """Recover an exact two-mode history on the real tower q-space."""

    def test_exact_history_on_actual_tower_material_points(self) -> None:
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
        operator = build_tower_equilibrium_operator(
            system
        )
        metric = operator.metric

        time = np.array(
            [0.0, 0.20, 0.45, 0.80, 1.20],
            dtype=np.float64,
        )
        q = np.arange(
            operator.n_material_points,
            dtype=np.float64,
        )
        p1 = (
            0.7 * np.sin(0.31 * q)
            + 0.2 * np.cos(0.17 * q)
        )
        p2 = (
            0.5 * np.cos(0.23 * q)
            - 0.3 * np.sin(0.11 * q)
        )
        s1 = operator.apply_spatial(p1).stress
        s2 = operator.apply_spatial(p2).stress

        amplitudes = np.array(
            [
                [0.0, 0.0],
                [0.08, -0.04],
                [0.20, 0.09],
                [0.14, 0.16],
                [0.30, 0.05],
            ],
            dtype=np.float64,
        )
        rates = np.zeros_like(amplitudes)
        rates[0] = (
            amplitudes[1] - amplitudes[0]
        ) / (
            time[1] - time[0]
        )
        rates[1:] = (
            np.diff(amplitudes, axis=0)
            / np.diff(time)[:, np.newaxis]
        )

        P = np.column_stack((p1, p2))
        S = np.column_stack((s1, s2))
        H_sigma = (
            2.0e-6
            + 1.0e-7
            * np.arange(time.size)[:, None]
            + 1.0e-9
            * q[None, :]
        )
        forcing = (
            rates @ P.T
            - H_sigma
            * (amplitudes @ S.T)
        )

        basis = PGDBasisTower(
            n_material_points=operator.n_material_points,
            n_time=time.size,
            modes=(
                PGDModeTower(
                    p1,
                    s1,
                    np.zeros(time.size),
                    np.zeros(time.size),
                ),
                PGDModeTower(
                    p2,
                    s2,
                    np.zeros(time.size),
                    np.zeros(time.size),
                ),
            ),
        )

        result = update_tower_pgd_time_functions(
            basis=basis,
            time=time,
            forcing=forcing,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
        )

        np.testing.assert_allclose(
            result.basis.temporal_amplitude_matrix(),
            amplitudes,
            rtol=0.0,
            atol=1.0e-10,
        )
        np.testing.assert_allclose(
            result.basis.temporal_rate_matrix(),
            rates,
            rtol=0.0,
            atol=1.0e-10,
        )
        self.assertLess(
            result.relative_residual,
            1.0e-9,
        )
        np.testing.assert_allclose(
            operator.equilibrium_residual(
                result.plastic_projection.stress
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

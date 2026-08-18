# -*- coding: utf-8 -*-
"""Unit tests for tower fixed-spatial-basis temporal updating."""

import unittest

import numpy as np

from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_pgd_time_update import (
    update_tower_pgd_time_functions,
)
from latin.tower_state import MaterialPointLayout


def _operator():
    layout = MaterialPointLayout(
        n_elements=1,
        n_gauss=2,
        n_fibers=2,
    )
    metric = MaterialPointMetric(
        np.array([1.0, 2.0, 1.5, 0.5]),
    )
    compatibility = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, -1.0],
        ],
        dtype=np.float64,
    )
    operator = TowerEquilibriumOperator(
        layout=layout,
        metric=metric,
        reference_modulus=np.array(
            [10.0, 20.0, 15.0, 25.0],
        ),
        compatibility_matrix=compatibility,
        free_dofs=np.array([1, 3]),
        n_dof=4,
    )
    return metric, operator


class TestTowerPGDTimeUpdate(unittest.TestCase):
    """Verify Eq. (58)-(59) temporal reuse on the tower metric."""

    def test_empty_basis_returns_minus_forcing(self) -> None:
        metric, operator = _operator()
        time = np.array([0.0, 0.5, 1.0])
        basis = PGDBasisTower(
            n_material_points=4,
            n_time=time.size,
        )
        forcing = np.array(
            [
                [1.0, -1.0, 0.5, 0.2],
                [0.5, 0.1, -0.2, 0.4],
                [0.2, -0.3, 0.1, 0.7],
            ]
        )
        H_sigma = np.full(
            forcing.shape,
            0.25,
        )

        result = update_tower_pgd_time_functions(
            basis=basis,
            time=time,
            forcing=forcing,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
        )

        self.assertEqual(result.n_modes, 0)
        self.assertFalse(result.time_functions_updated)
        self.assertFalse(result.reduced_converged)
        self.assertAlmostEqual(
            result.relative_residual,
            1.0,
            places=12,
        )
        np.testing.assert_allclose(
            result.mechanical_residual,
            -forcing,
        )
        np.testing.assert_allclose(
            result.plastic_strain_correction,
            np.zeros_like(forcing),
        )
        np.testing.assert_allclose(
            result.plastic_projection.stress,
            np.zeros_like(forcing),
        )

    def test_empty_basis_zero_forcing_is_converged(self) -> None:
        metric, operator = _operator()
        time = np.array([0.0, 0.5, 1.0])
        basis = PGDBasisTower(
            n_material_points=4,
            n_time=time.size,
        )
        forcing = np.zeros((time.size, 4))
        H_sigma = np.full(
            forcing.shape,
            0.25,
        )

        result = update_tower_pgd_time_functions(
            basis=basis,
            time=time,
            forcing=forcing,
            H_sigma=H_sigma,
            metric=metric,
            equilibrium_operator=operator,
        )

        self.assertTrue(result.reduced_converged)
        self.assertEqual(result.relative_residual, 0.0)
        self.assertEqual(result.forcing_norm, 0.0)

    def test_exact_single_mode_temporal_function_is_recovered(self) -> None:
        metric, operator = _operator()
        time = np.array(
            [0.0, 0.25, 0.50, 0.75, 1.00]
        )

        p = np.array(
            [0.6, -0.4, 0.9, -0.2],
        )
        s = operator.apply_spatial(p).stress

        exact_amplitude = np.array(
            [0.0, 0.1, 0.3, 0.2, 0.5]
        )
        exact_rate = np.zeros_like(
            exact_amplitude
        )
        exact_rate[0] = (
            exact_amplitude[1]
            - exact_amplitude[0]
        ) / (
            time[1] - time[0]
        )
        exact_rate[1:] = (
            np.diff(exact_amplitude)
            / np.diff(time)
        )

        H_sigma = np.empty((time.size, 4))
        for step in range(time.size):
            H_sigma[step] = (
                0.20
                + 0.02 * step
                + 0.01 * np.arange(4)
            )

        forcing = (
            np.outer(exact_rate, p)
            - H_sigma
            * np.outer(exact_amplitude, s)
        )

        basis = PGDBasisTower(
            n_material_points=4,
            n_time=time.size,
            modes=(
                PGDModeTower(
                    spatial_plastic_strain=p,
                    spatial_stress=s,
                    temporal_amplitude=np.zeros(
                        time.size
                    ),
                    temporal_rate=np.zeros(
                        time.size
                    ),
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

        recovered = result.basis.modes[0]
        np.testing.assert_allclose(
            recovered.temporal_amplitude,
            exact_amplitude,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            recovered.temporal_rate,
            exact_rate,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            result.mechanical_residual,
            np.zeros_like(forcing),
            rtol=0.0,
            atol=1.0e-11,
        )
        self.assertTrue(result.reduced_converged)
        self.assertTrue(result.time_functions_updated)
        self.assertLess(
            result.relative_residual,
            1.0e-11,
        )

        # Input basis is a persistent value and remains unchanged.
        np.testing.assert_allclose(
            basis.modes[0].temporal_amplitude,
            np.zeros(time.size),
        )

    def test_exact_two_mode_history_is_recovered(self) -> None:
        metric, operator = _operator()
        time = np.array(
            [0.0, 0.2, 0.5, 0.9, 1.4]
        )

        p1 = np.array([0.6, -0.4, 0.9, -0.2])
        p2 = np.array([-0.3, 0.8, 0.1, 0.7])
        s1 = operator.apply_spatial(p1).stress
        s2 = operator.apply_spatial(p2).stress

        amplitudes = np.array(
            [
                [0.0, 0.0],
                [0.1, -0.05],
                [0.25, 0.12],
                [0.15, 0.20],
                [0.35, 0.08],
            ]
        )
        rates = np.zeros_like(amplitudes)
        first_dt = time[1] - time[0]
        rates[0] = (
            amplitudes[1] - amplitudes[0]
        ) / first_dt
        rates[1:] = (
            np.diff(amplitudes, axis=0)
            / np.diff(time)[:, np.newaxis]
        )

        P = np.column_stack((p1, p2))
        S = np.column_stack((s1, s2))
        H_sigma = (
            0.22
            + 0.01
            * np.arange(time.size)[:, None]
            + 0.005
            * np.arange(4)[None, :]
        )
        forcing = (
            rates @ P.T
            - H_sigma * (amplitudes @ S.T)
        )

        basis = PGDBasisTower(
            n_material_points=4,
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
            atol=2.0e-12,
        )
        np.testing.assert_allclose(
            result.basis.temporal_rate_matrix(),
            rates,
            rtol=0.0,
            atol=2.0e-12,
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
            atol=1.0e-10,
        )
        self.assertLess(
            result.relative_residual,
            1.0e-10,
        )

    def test_inconsistent_stress_mode_is_rejected(self) -> None:
        metric, operator = _operator()
        time = np.array([0.0, 0.5, 1.0])
        p = np.array([0.6, -0.4, 0.9, -0.2])
        basis = PGDBasisTower(
            n_material_points=4,
            n_time=time.size,
            modes=(
                PGDModeTower(
                    spatial_plastic_strain=p,
                    spatial_stress=np.ones(4),
                    temporal_amplitude=np.zeros(
                        time.size
                    ),
                    temporal_rate=np.zeros(
                        time.size
                    ),
                ),
            ),
        )

        with self.assertRaises(ValueError):
            update_tower_pgd_time_functions(
                basis=basis,
                time=time,
                forcing=np.zeros(
                    (time.size, 4)
                ),
                H_sigma=np.ones(
                    (time.size, 4)
                ),
                metric=metric,
                equilibrium_operator=operator,
            )


if __name__ == "__main__":
    unittest.main()

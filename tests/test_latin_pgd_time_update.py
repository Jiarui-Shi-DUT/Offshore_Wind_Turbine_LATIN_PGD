# -*- coding: utf-8 -*-
"""Regression test for PGD temporal-function updating."""

import unittest

import numpy as np

from fem.bar_1d import create_uniform_bar_mesh
from latin.pgd_basis import PGDBasis1D, PGDMode1D
from latin.pgd_time_update import update_pgd_time_functions
from latin.search_directions import DescentSearchDirections


class TestPGDTimeUpdate(unittest.TestCase):
    """Verify recovery of an exactly representable temporal function."""

    def test_exact_single_mode_time_update(self) -> None:
        mesh = create_uniform_bar_mesh(
            length=3.0,
            n_elements=3,
        )
        area = 2.0
        time = np.array(
            [0.0, 0.25, 0.50, 0.75, 1.00],
            dtype=np.float64,
        )

        spatial_plastic_strain = np.array(
            [1.0, 1.5, 2.0],
            dtype=np.float64,
        )
        spatial_stress = np.array(
            [-0.4, -0.6, -0.8],
            dtype=np.float64,
        )

        exact_amplitude = np.array(
            [0.0, 0.1, 0.3, 0.2, 0.5],
            dtype=np.float64,
        )
        exact_rate = np.zeros_like(exact_amplitude)
        exact_rate[0] = (
            exact_amplitude[1] - exact_amplitude[0]
        ) / (time[1] - time[0])
        exact_rate[1:] = np.diff(exact_amplitude) / np.diff(time)

        shape = (time.size, mesh.n_elements)
        H_sigma = np.full(shape, 0.25, dtype=np.float64)
        directions = DescentSearchDirections(
            H_sigma=H_sigma,
            H_beta=np.ones(shape, dtype=np.float64),
            H_R_bar=np.ones(shape, dtype=np.float64),
            b_damage=np.zeros(shape, dtype=np.float64),
            regularization=0.15,
        )

        forcing = (
            np.outer(exact_rate, spatial_plastic_strain)
            - H_sigma
            * np.outer(exact_amplitude, spatial_stress)
        )

        basis = PGDBasis1D(
            n_elements=mesh.n_elements,
            n_time=time.size,
            modes=[
                PGDMode1D(
                    spatial_plastic_strain=spatial_plastic_strain,
                    spatial_stress=spatial_stress,
                    temporal_amplitude=np.zeros(time.size),
                    temporal_rate=np.zeros(time.size),
                )
            ],
        )

        result = update_pgd_time_functions(
            basis=basis,
            time=time,
            directions=directions,
            forcing=forcing,
            mesh=mesh,
            area=area,
        )

        recovered_mode = result.basis.modes[0]

        np.testing.assert_allclose(
            recovered_mode.temporal_amplitude,
            exact_amplitude,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            recovered_mode.temporal_rate,
            exact_rate,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            result.residual,
            np.zeros(shape, dtype=np.float64),
            rtol=0.0,
            atol=1.0e-12,
        )

        self.assertEqual(result.n_modes, 1)
        self.assertLess(result.weighted_residual_norm, 1.0e-11)
        self.assertLess(result.relative_residual, 1.0e-11)
        self.assertTrue(
            np.all(np.isfinite(result.condition_history))
        )

        # The input basis must remain unchanged.
        np.testing.assert_allclose(
            basis.modes[0].temporal_amplitude,
            np.zeros(time.size),
            rtol=0.0,
            atol=0.0,
        )


if __name__ == "__main__":
    unittest.main()

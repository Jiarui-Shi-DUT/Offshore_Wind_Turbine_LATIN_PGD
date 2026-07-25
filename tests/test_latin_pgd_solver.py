# -*- coding: utf-8 -*-
"""Regression tests for the adaptive one-dimensional LATIN-PGD solver."""

import unittest

import numpy as np

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
    prescribed_displacement,
)
from fem.bar_1d import create_uniform_bar_mesh
from latin.initialization import compute_elastic_initialization
from latin.pgd_solver import (
    PGDTerminationReason,
    solve_latin_pgd,
)
from latin.state import LatinState


class TestLatinPGDSolver(unittest.TestCase):
    """Validate the first adaptive LATIN-PGD outer iteration."""

    def test_one_iteration_advances_after_positive_saturation_gain(
        self,
    ) -> None:
        """
        A positive saturation gain above 0.1 must advance LATIN.

        The first reduced global stage starts from the normalised reference
        indicator xi_0 = 1.0.  Its first PGD mode gives xi_1 about 0.1784,
        hence zeta is positive and larger than the enrichment threshold.
        The candidate must therefore be accepted and the solver must stop
        only because max_iterations is one.
        """
        configuration = BenchmarkConfiguration()
        mesh = create_uniform_bar_mesh(
            configuration.length,
            configuration.n_elements,
        )
        materials = create_three_material_distribution(
            configuration.n_elements
        )
        time = create_time_grid(
            configuration.total_time,
            configuration.time_step,
        )
        displacement = prescribed_displacement(
            time,
            configuration.displacement_amplitude,
            configuration.period,
        )

        elastic = compute_elastic_initialization(
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            time=time,
            right_displacement=displacement,
        )
        initial_state = LatinState.from_elastic_initialization(
            elastic,
            materials,
        )

        result = solve_latin_pgd(
            initial_state=initial_state,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            tolerance=1.0e-12,
            max_iterations=1,
            relaxation=0.8,
            reduced_tolerance=1.0e-3,
            max_enrichments_per_iteration=1,
        )

        self.assertEqual(
            result.termination_reason,
            PGDTerminationReason.MAX_ITERATIONS,
        )
        self.assertFalse(result.converged)
        self.assertFalse(result.saturated)
        self.assertEqual(result.iterations, 1)
        self.assertEqual(result.trial_evaluations, 1)

        self.assertEqual(result.basis.n_modes, 1)
        self.assertEqual(result.total_modes_added, 1)
        np.testing.assert_array_equal(
            result.modes_added_history,
            np.array([1], dtype=np.int64),
        )
        np.testing.assert_array_equal(
            result.trial_basis_size_history,
            np.array([1], dtype=np.int64),
        )

        np.testing.assert_allclose(
            result.baseline_indicator_history,
            np.array([1.0]),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            result.trial_reduced_residual_history,
            np.array([0.03120500860641662]),
            rtol=0.0,
            atol=1.0e-10,
        )
        np.testing.assert_allclose(
            result.trial_indicator_history,
            np.array([0.17839546]),
            rtol=0.0,
            atol=1.0e-8,
        )
        np.testing.assert_allclose(
            result.indicator_history,
            result.trial_indicator_history,
            rtol=0.0,
            atol=0.0,
        )

        expected_zeta = (
            result.baseline_indicator_history[0]
            - result.trial_indicator_history[0]
        ) / (
            result.baseline_indicator_history[0]
            + result.trial_indicator_history[0]
        )
        self.assertGreater(expected_zeta, 0.1)
        self.assertAlmostEqual(
            result.saturation_history[0],
            expected_zeta,
            places=12,
        )

        self.assertTrue(np.all(np.isfinite(result.state.stress)))
        self.assertTrue(
            np.all(np.isfinite(result.state.elastic_strain))
        )
        self.assertTrue(
            np.all(np.isfinite(result.state.damage))
        )

        maximum_stress_spread = float(
            np.max(np.ptp(result.state.stress, axis=1))
        )
        self.assertLess(maximum_stress_spread, 1.0e-10)


if __name__ == "__main__":
    unittest.main()

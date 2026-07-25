# -*- coding: utf-8 -*-
"""Regression test for the adaptive LATIN-PGD outer solver."""

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
    """Verify one complete adaptive LATIN-PGD outer iteration."""

    def test_one_iteration_stops_at_enrichment_limit(self) -> None:
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
        initial_state = LatinState.from_elastic_initialization(
            initialization=elastic,
            materials=materials,
        )

        result = solve_latin_pgd(
            initial_state=initial_state,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            tolerance=1.0e-12,
            max_iterations=1,
            relaxation=0.8,
            reduced_tolerance=1.0e-2,
            max_enrichments_per_iteration=1,
        )

        self.assertEqual(
            result.termination_reason,
            PGDTerminationReason.MAX_ENRICHMENTS,
        )
        self.assertFalse(result.converged)
        self.assertFalse(result.saturated)
        self.assertEqual(result.iterations, 1)
        self.assertEqual(result.trial_evaluations, 1)

        self.assertEqual(result.basis.n_modes, 1)
        self.assertEqual(result.total_modes_added, 1)

        self.assertEqual(
            result.baseline_indicator_history.shape,
            (1,),
        )
        self.assertEqual(
            result.trial_indicator_history.shape,
            (1,),
        )
        self.assertEqual(
            result.saturation_history.shape,
            (1,),
        )
        self.assertEqual(
            result.indicator_history.shape,
            (1,),
        )

        self.assertAlmostEqual(
            result.baseline_indicator_history[0],
            0.13552558,
            places=7,
        )
        self.assertAlmostEqual(
            result.trial_indicator_history[0],
            0.17839546,
            places=7,
        )
        self.assertAlmostEqual(
            result.saturation_history[0],
            -0.13656262,
            places=7,
        )
        self.assertLess(result.saturation_history[0], 0.0)

        self.assertTrue(
            np.all(np.isfinite(result.state.stress))
        )
        self.assertTrue(
            np.all(np.isfinite(result.trial_indicator_history))
        )
        self.assertLessEqual(
            float(np.max(np.ptp(result.state.stress, axis=1))),
            1.0e-12,
        )


if __name__ == "__main__":
    unittest.main()

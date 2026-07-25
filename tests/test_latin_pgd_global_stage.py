# -*- coding: utf-8 -*-
"""Regression test for the adaptive one-dimensional PGD global stage."""

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
from latin.local_stage import solve_local_stage
from latin.pgd_basis import PGDBasis1D
from latin.pgd_global_stage import solve_pgd_global_stage
from latin.search_directions import (
    compute_descent_search_directions,
)
from latin.state import LatinState


class TestPGDGlobalStage(unittest.TestCase):
    """Verify one-mode enrichment of the first LATIN global stage."""

    def test_empty_basis_adds_one_equilibrated_mode(self) -> None:
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
        global_state = LatinState.from_elastic_initialization(
            initialization=elastic,
            materials=materials,
        )
        local_state = solve_local_stage(
            global_state=global_state,
            materials=materials,
        )
        directions = compute_descent_search_directions(
            local_state=local_state,
            materials=materials,
        )
        empty_basis = PGDBasis1D(
            n_elements=mesh.n_elements,
            n_time=time.size,
        )

        result = solve_pgd_global_stage(
            global_state=global_state,
            local_state=local_state,
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            basis=empty_basis,
            reduced_tolerance=1.0e-2,
            max_new_modes=1,
            allow_enrichment=True,
            latin_iteration=1,
        )

        self.assertEqual(result.n_modes, 1)
        self.assertEqual(result.modes_added, 1)
        self.assertTrue(result.time_functions_updated)
        self.assertFalse(result.reduced_converged)

        self.assertEqual(result.state.field_shape, (2001, 90))
        self.assertEqual(result.residual_history.shape, (2,))
        self.assertAlmostEqual(
            result.residual_history[0],
            1.0,
            places=12,
        )
        self.assertAlmostEqual(
            result.relative_residual,
            0.03120500860641662,
            places=10,
        )
        self.assertLess(
            result.relative_residual,
            result.residual_history[0],
        )

        self.assertTrue(
            np.all(np.isfinite(result.state.stress))
        )
        self.assertTrue(
            np.all(np.isfinite(result.mechanical_residual))
        )
        self.assertLessEqual(
            float(np.max(np.ptp(result.state.stress, axis=1))),
            1.0e-12,
        )

        # The input basis is copied; enrichment must not mutate it.
        self.assertEqual(empty_basis.n_modes, 0)


if __name__ == "__main__":
    unittest.main()

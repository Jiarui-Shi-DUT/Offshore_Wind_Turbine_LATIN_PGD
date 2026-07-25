# -*- coding: utf-8 -*-
"""Regression tests for residual-driven PGD basis enrichment."""

import unittest

import numpy as np

from examples.three_material_bar import (
    BenchmarkConfiguration,
    create_three_material_distribution,
    create_time_grid,
    prescribed_displacement,
)
from fem.bar_1d import create_uniform_bar_mesh
from latin.equilibrium_operator import apply_equilibrium_operator
from latin.initialization import compute_elastic_initialization
from latin.local_stage import solve_local_stage
from latin.pgd_basis import PGDBasis1D
from latin.pgd_enrichment import enrich_pgd_basis_once
from latin.pgd_global_stage import solve_pgd_global_stage
from latin.search_directions import (
    DescentSearchDirections,
    compute_descent_search_directions,
)
from latin.state import LatinState


class TestPGDEnrichment(unittest.TestCase):
    """Verify residual reduction by one separated PGD pair."""

    def test_single_mode_enrichment_reduces_residual(self) -> None:
        mesh = create_uniform_bar_mesh(
            length=3.0,
            n_elements=3,
        )
        materials = create_three_material_distribution(
            mesh.n_elements
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
        spatial_stress = apply_equilibrium_operator(
            mesh=mesh,
            area=area,
            materials=materials,
            source_strain=spatial_plastic_strain,
        ).stress[0, :]

        temporal_amplitude = np.array(
            [0.0, 0.1, 0.3, 0.2, 0.5],
            dtype=np.float64,
        )
        temporal_rate = np.zeros_like(temporal_amplitude)
        temporal_rate[0] = (
            temporal_amplitude[1] - temporal_amplitude[0]
        ) / (time[1] - time[0])
        temporal_rate[1:] = (
            np.diff(temporal_amplitude) / np.diff(time)
        )

        shape = (time.size, mesh.n_elements)
        H_sigma = np.full(
            shape,
            0.25,
            dtype=np.float64,
        )
        directions = DescentSearchDirections(
            H_sigma=H_sigma,
            H_beta=np.ones(shape, dtype=np.float64),
            H_R_bar=np.ones(shape, dtype=np.float64),
            b_damage=np.zeros(shape, dtype=np.float64),
            regularization=0.15,
        )

        residual = -(
            np.outer(
                temporal_rate,
                spatial_plastic_strain,
            )
            - H_sigma
            * np.outer(
                temporal_amplitude,
                spatial_stress,
            )
        )

        basis = PGDBasis1D(
            n_elements=mesh.n_elements,
            n_time=time.size,
        )

        result = enrich_pgd_basis_once(
            basis=basis,
            time=time,
            residual=residual,
            directions=directions,
            mesh=mesh,
            area=area,
            materials=materials,
            max_fixed_point_iterations=30,
        )

        self.assertTrue(result.converged)
        self.assertTrue(result.accepted)
        self.assertGreaterEqual(result.iterations, 2)
        self.assertLessEqual(result.iterations, 10)

        self.assertEqual(
            result.mode.n_elements,
            mesh.n_elements,
        )
        self.assertEqual(
            result.mode.n_time,
            time.size,
        )
        self.assertEqual(
            result.residual.shape,
            shape,
        )

        self.assertTrue(np.all(np.isfinite(result.residual)))
        self.assertTrue(
            np.all(np.isfinite(result.fixed_point_history[1:]))
        )
        self.assertLess(
            result.residual_norm_after
            / result.residual_norm_before,
            1.0e-9,
        )
        self.assertGreater(
            result.residual_reduction,
            0.999999999,
        )

        # The enrichment routine only returns a candidate mode.
        self.assertEqual(basis.n_modes, 0)

    def test_second_mode_line_search_prevents_residual_growth(
        self,
    ) -> None:
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

        first_stage = solve_pgd_global_stage(
            global_state=global_state,
            local_state=local_state,
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            basis=PGDBasis1D(
                n_elements=mesh.n_elements,
                n_time=time.size,
            ),
            reduced_tolerance=1.0e-3,
            max_new_modes=1,
            allow_enrichment=True,
            latin_iteration=1,
        )

        second_mode = enrich_pgd_basis_once(
            basis=first_stage.basis,
            time=time,
            residual=first_stage.mechanical_residual,
            directions=directions,
            mesh=mesh,
            area=configuration.area,
            materials=materials,
            iteration_added=1,
        )

        self.assertTrue(second_mode.converged)
        self.assertTrue(second_mode.accepted)
        self.assertLess(
            second_mode.residual_norm_after,
            second_mode.residual_norm_before,
        )
        self.assertLess(second_mode.relative_residual, 1.0)
        self.assertGreater(second_mode.residual_reduction, 0.0)
        self.assertAlmostEqual(
            second_mode.relative_residual,
            0.9360938296734366,
            places=10,
        )
        self.assertEqual(first_stage.basis.n_modes, 1)


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""Actual coarse-tower integration tests for elastic LATIN initialization."""

import unittest
import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.tower_loading import create_reversed_top_force_history
from fem.tower_response_2d import recover_elastic_tower_response
from fem.tower_system_2d import (
    solve_elastic_tower_top_force,
    top_horizontal_load_vector,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from latin.tower_initialization import (
    compute_tower_elastic_initialization,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerElasticInitializationIntegration(unittest.TestCase):
    def setUp(self) -> None:
        self.geometry = LinearTaperedTowerGeometry(
            height=87.6,
            base_outer_diameter=6.0,
            top_outer_diameter=3.87,
            base_thickness=0.027,
            top_thickness=0.019,
        )
        self.mesh = create_uniform_vertical_tower_mesh(
            height=self.geometry.height,
            n_elements=10,
        )
        self.material = MaterialParameters()
        self.system = ViscoplasticDamageTowerSystem2D(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            material=self.material,
            n_gauss=2,
            n_circumferential=16,
            n_radial=1,
        )
        self.operator = build_tower_equilibrium_operator(
            self.system
        )
        self.loading = create_reversed_top_force_history(
            force_amplitude=1.0e6,
            period=10.0,
            n_cycles=1,
            increments_per_cycle=40,
        )
        self.load_vectors = np.stack(
            [
                top_horizontal_load_vector(
                    mesh=self.mesh,
                    horizontal_force=float(force),
                )
                for force in self.loading.forces
            ],
            axis=0,
        )

    def _run(self):
        return compute_tower_elastic_initialization(
            time=self.loading.times,
            load_vectors=self.load_vectors,
            materials=self.material,
            equilibrium_operator=self.operator,
            stress_to_force_factor=(
                self.system.stress_to_force_factor
            ),
        )

    def test_actual_tower_history_is_globally_admissible(self) -> None:
        result = self._run()
        state = result.state

        self.assertEqual(
            self.operator.n_material_points,
            10 * 2 * 16,
        )
        self.assertEqual(
            state.field_shape,
            (41, 10 * 2 * 16),
        )
        self.assertEqual(
            result.displacements.shape,
            (41, self.mesh.n_dof),
        )

        np.testing.assert_allclose(
            state.stress[0],
            np.zeros(state.n_material_points),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            state.elastic_strain[0],
            np.zeros(state.n_material_points),
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            state.stress[-1],
            np.zeros(state.n_material_points),
            rtol=0.0,
            atol=0.0,
        )

        self.assertLessEqual(
            result.maximum_free_equilibrium_residual,
            1.0e-9 * float(np.max(np.abs(self.loading.forces))),
        )

        np.testing.assert_allclose(
            state.elastic_strain,
            result.displacements[:, self.operator.free_dofs]
            @ self.operator.compatibility_matrix.T,
            rtol=1.0e-12,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            state.stress,
            state.elastic_strain
            * self.operator.reference_modulus[np.newaxis, :],
            rtol=1.0e-12,
            atol=1.0e-12,
        )

    def test_positive_peak_matches_independent_elastic_fem(self) -> None:
        result = self._run()
        quarter = self.loading.increments_per_cycle // 4
        peak_force = float(self.loading.forces[quarter])
        self.assertAlmostEqual(peak_force, 1.0e6)

        assembly, solution = solve_elastic_tower_top_force(
            mesh=self.mesh,
            tower_geometry=self.geometry,
            elastic_modulus=self.material.E * 1.0e6,
            horizontal_force=peak_force,
            n_gauss=2,
            n_circumferential=16,
            n_radial=1,
        )
        response = recover_elastic_tower_response(
            assembly=assembly,
            solution=solution,
            tower_geometry=self.geometry,
            elastic_modulus=self.material.E * 1.0e6,
            n_circumferential=16,
            n_radial=1,
        )

        np.testing.assert_allclose(
            result.displacements[quarter],
            solution.displacements,
            rtol=1.0e-10,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            result.state.stress[quarter],
            response.fiber_stresses.reshape(-1) / 1.0e6,
            rtol=1.0e-10,
            atol=1.0e-10,
        )


if __name__ == "__main__":
    unittest.main()

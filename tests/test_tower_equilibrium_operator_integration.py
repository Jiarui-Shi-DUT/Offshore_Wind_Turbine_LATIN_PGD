# -*- coding: utf-8 -*-
"""
Integration test against the repository's actual tower FEM discretisation.

This test verifies that the material-point operator H^T M C0 H reconstructs
the same free-DOF reference stiffness as the existing elastic fiber-beam tower
assembly when both use the same mesh, Gauss rule, fiber discretisation, and
reference modulus.
"""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    create_uniform_vertical_tower_mesh,
)
from fem.tower_system_2d import (
    assemble_elastic_tower_stiffness,
    cantilever_base_fixed_dofs,
    free_dofs_from_fixed,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
)
from latin.tower_equilibrium_operator import (
    build_tower_equilibrium_operator,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerEquilibriumOperatorIntegration(unittest.TestCase):
    """Compare the LATIN reference operator with the actual tower FEM."""

    def test_reference_stiffness_matches_elastic_fiber_beam_assembly(
        self,
    ) -> None:
        geometry = LinearTaperedTowerGeometry(
            height=12.0,
            base_outer_diameter=2.0,
            top_outer_diameter=1.6,
            base_thickness=0.05,
            top_thickness=0.04,
        )
        mesh = create_uniform_vertical_tower_mesh(
            height=geometry.height,
            n_elements=2,
        )
        material = MaterialParameters(E=210_000.0)

        n_gauss = 2
        n_circumferential = 8
        n_radial = 1

        nonlinear_system = ViscoplasticDamageTowerSystem2D(
            mesh=mesh,
            tower_geometry=geometry,
            material=material,
            n_gauss=n_gauss,
            n_circumferential=n_circumferential,
            n_radial=n_radial,
        )

        operator = build_tower_equilibrium_operator(
            nonlinear_system,
        )

        elastic_assembly = assemble_elastic_tower_stiffness(
            mesh=mesh,
            tower_geometry=geometry,
            elastic_modulus=material.E,
            n_gauss=n_gauss,
            n_circumferential=n_circumferential,
            n_radial=n_radial,
        )

        fixed = cantilever_base_fixed_dofs(mesh)
        free = free_dofs_from_fixed(
            mesh.n_dof,
            fixed,
        )
        expected_reduced = elastic_assembly.global_stiffness[
            np.ix_(free, free)
        ]

        self.assertEqual(
            operator.layout.n_elements,
            mesh.n_elements,
        )
        self.assertEqual(
            operator.layout.n_gauss,
            n_gauss,
        )
        self.assertEqual(
            operator.layout.n_fibers,
            n_circumferential * n_radial,
        )
        np.testing.assert_array_equal(
            operator.free_dofs,
            free,
        )

        np.testing.assert_allclose(
            operator.reduced_stiffness,
            expected_reduced,
            rtol=1.0e-12,
            atol=1.0e-10,
        )

        expected_measure = 0.0
        for element in nonlinear_system.elements:
            for gauss_index, section_wrapper in enumerate(
                element.sections
            ):
                expected_measure += (
                    element.jacobian
                    * element.gauss_weights[gauss_index]
                    * section_wrapper.section.area
                )

        self.assertAlmostEqual(
            operator.metric.total_measure,
            expected_measure,
            places=12,
        )

    def test_real_tower_projection_returns_free_dof_equilibrium(
        self,
    ) -> None:
        geometry = LinearTaperedTowerGeometry(
            height=8.0,
            base_outer_diameter=1.5,
            top_outer_diameter=1.3,
            base_thickness=0.04,
            top_thickness=0.035,
        )
        mesh = create_uniform_vertical_tower_mesh(
            height=geometry.height,
            n_elements=2,
        )
        material = MaterialParameters()

        system = ViscoplasticDamageTowerSystem2D(
            mesh=mesh,
            tower_geometry=geometry,
            material=material,
            n_gauss=2,
            n_circumferential=8,
            n_radial=1,
        )
        operator = build_tower_equilibrium_operator(system)

        q = np.arange(
            operator.n_material_points,
            dtype=np.float64,
        )
        source = 1.0e-4 * np.sin(0.37 * q)

        projection = operator.apply_spatial(source)
        residual = operator.equilibrium_residual(
            projection.stress
        )

        np.testing.assert_allclose(
            residual,
            np.zeros(operator.n_free_dofs),
            rtol=0.0,
            atol=1.0e-8,
        )
        self.assertTrue(
            np.all(np.isfinite(projection.compatible_strain))
        )
        self.assertTrue(
            np.all(np.isfinite(projection.stress))
        )


if __name__ == "__main__":
    unittest.main()

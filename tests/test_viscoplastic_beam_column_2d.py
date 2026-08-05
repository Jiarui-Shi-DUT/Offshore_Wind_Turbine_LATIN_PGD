# -*- coding: utf-8 -*-
"""Regression tests for the viscoplastic-damage 2D beam element."""

import unittest

import numpy as np

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    compute_elastic_beam_element_stiffness,
)
from fem.viscoplastic_beam_column_2d import (
    ViscoplasticDamageBeamElement2D,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestViscoplasticDamageBeamElement2D(unittest.TestCase):
    """Verify kinematics, integration, tangent, and state management."""

    def setUp(self) -> None:
        self.length = 2.0
        self.coordinates = np.array(
            [
                [0.0, 0.0],
                [0.0, self.length],
            ],
            dtype=np.float64,
        )
        self.geometry = LinearTaperedTowerGeometry(
            height=self.length,
            base_outer_diameter=1.0,
            top_outer_diameter=1.0,
            base_thickness=0.1,
            top_thickness=0.1,
        )
        self.elastic_material = MaterialParameters(
            sigma_y=1.0e9,
            k_damage=0.0,
        )

    def create_element(
        self,
        material: MaterialParameters = None,
    ) -> ViscoplasticDamageBeamElement2D:
        """Create one constant-section vertical beam element."""
        if material is None:
            material = self.elastic_material

        return ViscoplasticDamageBeamElement2D(
            node_coordinates=self.coordinates,
            tower_axis_start=0.0,
            tower_axis_end=self.length,
            tower_geometry=self.geometry,
            material=material,
            n_gauss=4,
            n_circumferential=8,
            n_radial=1,
        )

    @staticmethod
    def local_to_global(
        element: ViscoplasticDamageBeamElement2D,
        local_displacements: np.ndarray,
    ) -> np.ndarray:
        """Convert a local element displacement vector to global DOFs."""
        return element.transformation.T @ local_displacements

    def test_initial_element_state(self) -> None:
        """Every Gauss point must own an independent zero-state section."""
        element = self.create_element()
        response = element.committed_response

        self.assertEqual(len(element.sections), 4)
        self.assertEqual(response.generalized_strains.shape, (4, 2))
        self.assertEqual(response.section_resultants.shape, (4, 2))
        self.assertEqual(response.local_internal_force.shape, (6,))
        self.assertEqual(response.global_internal_force.shape, (6,))
        self.assertFalse(element.has_uncommitted_trial)

        np.testing.assert_allclose(
            response.global_displacements,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            response.generalized_strains,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            response.section_resultants,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

        for section in element.sections:
            np.testing.assert_allclose(
                section.committed_states,
                0.0,
                rtol=0.0,
                atol=0.0,
            )

    def test_uniform_elastic_axial_response(self) -> None:
        """Uniform extension must recover constant strain and axial force."""
        element = self.create_element()
        axial_strain = 1.0e-4

        local_displacements = np.array(
            [
                0.0,
                0.0,
                0.0,
                axial_strain * self.length,
                0.0,
                0.0,
            ],
            dtype=np.float64,
        )
        global_displacements = self.local_to_global(
            element,
            local_displacements,
        )

        response = element.set_trial_displacements(
            time=0.1,
            global_displacements=global_displacements,
            compute_tangent=False,
        )

        area = element.sections[0].section.area
        expected_force = (
            self.elastic_material.E
            * 1.0e6
            * area
            * axial_strain
        )

        np.testing.assert_allclose(
            response.axial_strains,
            axial_strain,
            rtol=0.0,
            atol=1.0e-16,
        )
        np.testing.assert_allclose(
            response.curvatures,
            0.0,
            rtol=0.0,
            atol=1.0e-16,
        )
        np.testing.assert_allclose(
            response.axial_forces,
            expected_force,
            rtol=1.0e-12,
            atol=1.0e-6,
        )
        np.testing.assert_allclose(
            response.bending_moments,
            0.0,
            rtol=0.0,
            atol=1.0e-6,
        )
        np.testing.assert_allclose(
            response.local_internal_force,
            np.array(
                [
                    -expected_force,
                    0.0,
                    0.0,
                    expected_force,
                    0.0,
                    0.0,
                ],
                dtype=np.float64,
            ),
            rtol=1.0e-12,
            atol=1.0e-6,
        )

    def test_constant_elastic_curvature_response(self) -> None:
        """A quadratic transverse field must recover constant curvature."""
        element = self.create_element()
        curvature = 5.0e-4

        local_displacements = np.array(
            [
                0.0,
                0.0,
                0.0,
                0.0,
                0.5 * curvature * self.length ** 2,
                curvature * self.length,
            ],
            dtype=np.float64,
        )
        global_displacements = self.local_to_global(
            element,
            local_displacements,
        )

        response = element.set_trial_displacements(
            time=0.1,
            global_displacements=global_displacements,
            compute_tangent=False,
        )

        second_moment = (
            element.sections[0].section.second_moment_x
        )
        expected_moment = (
            self.elastic_material.E
            * 1.0e6
            * second_moment
            * curvature
        )

        np.testing.assert_allclose(
            response.axial_strains,
            0.0,
            rtol=0.0,
            atol=1.0e-16,
        )
        np.testing.assert_allclose(
            response.curvatures,
            curvature,
            rtol=0.0,
            atol=1.0e-16,
        )
        np.testing.assert_allclose(
            response.axial_forces,
            0.0,
            rtol=0.0,
            atol=1.0e-6,
        )
        np.testing.assert_allclose(
            response.bending_moments,
            expected_moment,
            rtol=1.0e-12,
            atol=1.0e-6,
        )
        np.testing.assert_allclose(
            response.local_internal_force,
            np.array(
                [
                    0.0,
                    0.0,
                    -expected_moment,
                    0.0,
                    0.0,
                    expected_moment,
                ],
                dtype=np.float64,
            ),
            rtol=1.0e-12,
            atol=1.0e-6,
        )

    def test_elastic_tangent_matches_existing_beam_element(self) -> None:
        """The nonlinear element tangent must recover the elastic baseline."""
        element = self.create_element()

        local_displacements = np.array(
            [
                0.0,
                0.0,
                0.0,
                2.0e-4,
                3.0e-4,
                1.0e-4,
            ],
            dtype=np.float64,
        )
        global_displacements = self.local_to_global(
            element,
            local_displacements,
        )

        response = element.set_trial_displacements(
            time=0.1,
            global_displacements=global_displacements,
            compute_tangent=True,
        )

        elastic_result = compute_elastic_beam_element_stiffness(
            node_coordinates=self.coordinates,
            tower_axis_start=0.0,
            tower_axis_end=self.length,
            tower_geometry=self.geometry,
            elastic_modulus=self.elastic_material.E * 1.0e6,
            n_gauss=4,
            n_circumferential=8,
            n_radial=1,
        )

        self.assertIsNotNone(response.local_tangent)
        self.assertIsNotNone(response.global_tangent)

        np.testing.assert_allclose(
            response.local_tangent,
            elastic_result.local_stiffness,
            rtol=1.0e-7,
            atol=1.0e-2,
        )
        np.testing.assert_allclose(
            response.global_tangent,
            elastic_result.global_stiffness,
            rtol=1.0e-7,
            atol=1.0e-2,
        )
        np.testing.assert_allclose(
            response.local_tangent,
            response.local_tangent.T,
            rtol=0.0,
            atol=1.0e-8,
        )

    def test_repeated_trial_does_not_accumulate_damage(self) -> None:
        """Repeated nonlinear trials must restart from the last commit."""
        element = self.create_element(MaterialParameters())
        axial_strain = 1.2e-3

        local_displacements = np.array(
            [
                0.0,
                0.0,
                0.0,
                axial_strain * self.length,
                0.0,
                0.0,
            ],
            dtype=np.float64,
        )
        global_displacements = self.local_to_global(
            element,
            local_displacements,
        )

        first = element.set_trial_displacements(
            time=0.1,
            global_displacements=global_displacements,
            compute_tangent=False,
        )
        first_states = np.stack(
            [
                response.fiber_states
                for response in first.section_responses
            ],
            axis=0,
        )

        second = element.set_trial_displacements(
            time=0.1,
            global_displacements=global_displacements,
            compute_tangent=False,
        )
        second_states = np.stack(
            [
                response.fiber_states
                for response in second.section_responses
            ],
            axis=0,
        )

        np.testing.assert_allclose(
            second_states,
            first_states,
            rtol=0.0,
            atol=0.0,
        )
        self.assertGreater(first.maximum_damage, 0.0)

        for section in element.sections:
            np.testing.assert_allclose(
                section.committed_states,
                0.0,
                rtol=0.0,
                atol=0.0,
            )

    def test_commit_rollback_and_restart(self) -> None:
        """Element state operations must act at every Gauss point."""
        element = self.create_element(MaterialParameters())

        local_first = np.array(
            [
                0.0,
                0.0,
                0.0,
                1.2e-3 * self.length,
                0.0,
                0.0,
            ],
            dtype=np.float64,
        )
        global_first = self.local_to_global(
            element,
            local_first,
        )

        element.set_trial_displacements(
            time=0.1,
            global_displacements=global_first,
            compute_tangent=False,
        )
        committed = element.commit_state()
        committed_states = [
            section.committed_states
            for section in element.sections
        ]

        self.assertFalse(element.has_uncommitted_trial)
        self.assertGreater(committed.maximum_damage, 0.0)

        local_second = np.array(
            [
                0.0,
                0.0,
                0.0,
                -1.2e-3 * self.length,
                0.0,
                0.0,
            ],
            dtype=np.float64,
        )
        global_second = self.local_to_global(
            element,
            local_second,
        )

        element.set_trial_displacements(
            time=0.2,
            global_displacements=global_second,
            compute_tangent=False,
        )
        self.assertTrue(element.has_uncommitted_trial)

        element.revert_to_last_commit()

        self.assertFalse(element.has_uncommitted_trial)
        for section, expected in zip(
            element.sections,
            committed_states,
        ):
            np.testing.assert_allclose(
                section.trial_states,
                expected,
                rtol=0.0,
                atol=0.0,
            )

        restarted = element.revert_to_start()

        self.assertEqual(element.committed_time, 0.0)
        self.assertFalse(element.has_uncommitted_trial)
        np.testing.assert_allclose(
            restarted.global_displacements,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        for section in element.sections:
            np.testing.assert_allclose(
                section.committed_states,
                0.0,
                rtol=0.0,
                atol=0.0,
            )

    def test_invalid_operations_are_rejected(self) -> None:
        """Invalid displacement and commit operations must raise errors."""
        element = self.create_element()

        with self.assertRaises(RuntimeError):
            element.commit_state()

        with self.assertRaises(ValueError):
            element.set_trial_displacements(
                time=0.1,
                global_displacements=np.zeros(5),
            )

        with self.assertRaises(ValueError):
            element.set_trial_displacements(
                time=0.0,
                global_displacements=np.zeros(6),
            )


if __name__ == "__main__":
    unittest.main()

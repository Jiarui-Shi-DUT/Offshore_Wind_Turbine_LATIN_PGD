# -*- coding: utf-8 -*-
"""Regression tests for the viscoplastic-damage annular fiber section."""

import unittest

import numpy as np

from fem.fiber_section import create_annular_fiber_section
from fem.viscoplastic_fiber_section import (
    ViscoplasticDamageFiberSection,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestViscoplasticDamageFiberSection(unittest.TestCase):
    """Verify section integration and trial-state management."""

    def setUp(self) -> None:
        self.geometry = create_annular_fiber_section(
            outer_diameter=1.0,
            thickness=0.1,
            n_circumferential=8,
            n_radial=1,
        )
        self.elastic_material = MaterialParameters(
            sigma_y=1.0e9,
            k_damage=0.0,
        )

    def create_elastic_section(
        self,
    ) -> ViscoplasticDamageFiberSection:
        """Create a section that remains elastic in the test range."""
        return ViscoplasticDamageFiberSection(
            section=self.geometry,
            material=self.elastic_material,
        )

    def test_initial_state_is_zero_and_independent(self) -> None:
        """Every fiber must start with its own four-component state."""
        section = self.create_elastic_section()

        self.assertEqual(
            section.committed_states.shape,
            (self.geometry.n_fibers, 4),
        )
        np.testing.assert_allclose(
            section.committed_states,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        self.assertFalse(section.has_uncommitted_trial)

        returned_states = section.committed_states
        returned_states[0, 0] = 1.0

        np.testing.assert_allclose(
            section.committed_states,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

    def test_uniform_elastic_axial_response(self) -> None:
        """Zero curvature must produce uniform stress and no moment."""
        section = self.create_elastic_section()
        axial_strain = 1.0e-4

        response = section.set_trial_deformation(
            time=0.1,
            axial_strain=axial_strain,
            curvature=0.0,
        )

        expected_stress = self.elastic_material.E * axial_strain
        expected_force = (
            expected_stress * 1.0e6 * self.geometry.area
        )

        np.testing.assert_allclose(
            response.fiber_strains,
            axial_strain,
            rtol=0.0,
            atol=1.0e-16,
        )
        np.testing.assert_allclose(
            response.fiber_stresses,
            expected_stress,
            rtol=1.0e-13,
            atol=1.0e-13,
        )
        self.assertAlmostEqual(
            response.axial_force / expected_force,
            1.0,
            places=12,
        )
        self.assertAlmostEqual(
            response.bending_moment,
            0.0,
            places=6,
        )
        np.testing.assert_allclose(
            response.fiber_states,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

    def test_pure_elastic_bending_response(self) -> None:
        """Pure curvature must recover N = 0 and M = E I kappa."""
        section = self.create_elastic_section()
        curvature = 5.0e-4

        response = section.set_trial_deformation(
            time=0.1,
            axial_strain=0.0,
            curvature=curvature,
        )

        expected_strains = (
            -curvature * self.geometry.y_coordinates
        )
        expected_stresses = (
            self.elastic_material.E * expected_strains
        )
        expected_moment = (
            self.elastic_material.E
            * 1.0e6
            * self.geometry.second_moment_x
            * curvature
        )

        np.testing.assert_allclose(
            response.fiber_strains,
            expected_strains,
            rtol=0.0,
            atol=1.0e-16,
        )
        np.testing.assert_allclose(
            response.fiber_stresses,
            expected_stresses,
            rtol=1.0e-13,
            atol=1.0e-13,
        )
        self.assertAlmostEqual(
            response.axial_force,
            0.0,
            places=6,
        )
        self.assertAlmostEqual(
            response.bending_moment / expected_moment,
            1.0,
            places=12,
        )

    def test_repeated_trial_does_not_accumulate_state(self) -> None:
        """Repeated Newton/LATIN trials must restart from the last commit."""
        section = ViscoplasticDamageFiberSection(
            section=self.geometry,
            material=MaterialParameters(),
        )

        first = section.set_trial_deformation(
            time=0.1,
            axial_strain=1.2e-3,
            curvature=0.0,
        )
        second = section.set_trial_deformation(
            time=0.1,
            axial_strain=1.2e-3,
            curvature=0.0,
        )

        np.testing.assert_allclose(
            second.fiber_states,
            first.fiber_states,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            second.fiber_stresses,
            first.fiber_stresses,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            section.committed_states,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        self.assertGreater(
            float(np.max(np.abs(first.plastic_strains))),
            0.0,
        )
        self.assertGreater(first.maximum_damage, 0.0)

    def test_commit_rollback_and_restart(self) -> None:
        """Commit, rollback, and restart must preserve state semantics."""
        section = ViscoplasticDamageFiberSection(
            section=self.geometry,
            material=MaterialParameters(),
        )

        trial = section.set_trial_deformation(
            time=0.1,
            axial_strain=1.2e-3,
            curvature=0.0,
        )
        committed = section.commit_state()
        committed_states = committed.fiber_states.copy()

        self.assertFalse(section.has_uncommitted_trial)
        np.testing.assert_allclose(
            section.committed_states,
            committed_states,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            committed_states,
            trial.fiber_states,
            rtol=0.0,
            atol=0.0,
        )

        section.set_trial_deformation(
            time=0.2,
            axial_strain=-1.2e-3,
            curvature=0.0,
        )
        self.assertTrue(section.has_uncommitted_trial)
        self.assertFalse(
            np.array_equal(
                section.trial_states,
                committed_states,
            )
        )

        rolled_back = section.revert_to_last_commit()

        self.assertFalse(section.has_uncommitted_trial)
        np.testing.assert_allclose(
            rolled_back.fiber_states,
            committed_states,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            section.trial_states,
            committed_states,
            rtol=0.0,
            atol=0.0,
        )

        restarted = section.revert_to_start()

        self.assertEqual(section.committed_time, 0.0)
        self.assertFalse(section.has_uncommitted_trial)
        np.testing.assert_allclose(
            restarted.fiber_states,
            0.0,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            section.committed_states,
            0.0,
            rtol=0.0,
            atol=0.0,
        )

    def test_elastic_numerical_section_tangent(self) -> None:
        """The numerical tangent must recover EA and EI in elasticity."""
        section = self.create_elastic_section()

        response = section.set_trial_deformation(
            time=0.1,
            axial_strain=1.0e-4,
            curvature=2.0e-4,
            compute_tangent=True,
        )

        self.assertIsNotNone(response.tangent)
        expected = np.array(
            [
                [
                    self.elastic_material.E
                    * 1.0e6
                    * self.geometry.area,
                    0.0,
                ],
                [
                    0.0,
                    self.elastic_material.E
                    * 1.0e6
                    * self.geometry.second_moment_x,
                ],
            ],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            response.tangent,
            expected,
            rtol=1.0e-7,
            atol=1.0e-2,
        )

    def test_invalid_state_operations_are_rejected(self) -> None:
        """Invalid time progression and empty commits must raise errors."""
        section = self.create_elastic_section()

        with self.assertRaises(RuntimeError):
            section.commit_state()

        with self.assertRaises(ValueError):
            section.set_trial_deformation(
                time=0.0,
                axial_strain=0.0,
                curvature=0.0,
            )

        with self.assertRaises(ValueError):
            section.set_trial_deformation(
                time=-0.1,
                axial_strain=0.0,
                curvature=0.0,
            )


if __name__ == "__main__":
    unittest.main()

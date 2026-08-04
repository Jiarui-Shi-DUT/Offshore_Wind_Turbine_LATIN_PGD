# -*- coding: utf-8 -*-
"""Regression tests for circular annular fiber-section response."""

import unittest

import numpy as np

from fem.fiber_section import (
    compute_fiber_strains,
    create_annular_fiber_section,
    evaluate_linear_elastic_section,
    integrate_fiber_stresses,
)


class TestAnnularFiberSection(unittest.TestCase):
    """Verify annular-sector fiber geometry and section response."""

    def setUp(self) -> None:
        self.section = create_annular_fiber_section(
            outer_diameter=6.0,
            thickness=0.027,
            n_circumferential=32,
            n_radial=2,
        )
        self.elastic_modulus = 210.0e9

    def test_default_discretization(self) -> None:
        section = create_annular_fiber_section(
            outer_diameter=6.0,
            thickness=0.027,
        )

        self.assertEqual(section.n_circumferential, 32)
        self.assertEqual(section.n_radial, 2)
        self.assertEqual(section.n_fibers, 64)
        self.assertEqual(section.coordinates.shape, (64, 2))

        self.assertAlmostEqual(
            section.outer_radius,
            3.0,
            places=15,
        )
        self.assertAlmostEqual(
            section.inner_radius,
            2.973,
            places=15,
        )

    def test_area_and_centroid(self) -> None:
        np.testing.assert_allclose(
            self.section.area,
            self.section.exact_area,
            rtol=0.0,
            atol=1.0e-14,
        )
        self.assertAlmostEqual(
            self.section.centroid_x,
            0.0,
            places=14,
        )
        self.assertAlmostEqual(
            self.section.centroid_y,
            0.0,
            places=14,
        )
        self.assertAlmostEqual(
            self.section.first_moment_y,
            0.0,
            places=14,
        )

    def test_centroids_align_with_coordinate_axes(self) -> None:
        first_layer = self.section.radial_indices == 0
        x = self.section.x_coordinates[first_layer]
        y = self.section.y_coordinates[first_layer]

        self.assertAlmostEqual(y[0], 0.0, places=14)
        self.assertGreater(x[0], 0.0)

        self.assertAlmostEqual(x[8], 0.0, places=14)
        self.assertGreater(y[8], 0.0)

        self.assertAlmostEqual(y[16], 0.0, places=14)
        self.assertLess(x[16], 0.0)

        self.assertAlmostEqual(x[24], 0.0, places=14)
        self.assertLess(y[24], 0.0)

    def test_second_moment_symmetry_and_accuracy(self) -> None:
        np.testing.assert_allclose(
            self.section.second_moment_x,
            self.section.second_moment_y,
            rtol=0.0,
            atol=1.0e-13,
        )
        self.assertAlmostEqual(
            self.section.product_moment_xy,
            0.0,
            places=14,
        )

        relative_error = abs(
            self.section.second_moment_x
            - self.section.exact_second_moment
        ) / self.section.exact_second_moment

        self.assertLess(relative_error, 5.0e-3)

    def test_fiber_strain_field(self) -> None:
        axial_strain = 2.0e-4
        curvature = 1.5e-4

        strains = compute_fiber_strains(
            section=self.section,
            axial_strain=axial_strain,
            curvature=curvature,
        )

        expected = (
            axial_strain
            - curvature * self.section.y_coordinates
        )

        np.testing.assert_allclose(
            strains,
            expected,
            rtol=0.0,
            atol=1.0e-16,
        )

    def test_uniform_stress_integration(self) -> None:
        uniform_stress = 125.0e6
        stresses = np.full(
            self.section.n_fibers,
            uniform_stress,
            dtype=np.float64,
        )

        resultants = integrate_fiber_stresses(
            section=self.section,
            fiber_stresses=stresses,
        )

        np.testing.assert_allclose(
            resultants[0],
            uniform_stress * self.section.area,
            rtol=1.0e-14,
            atol=0.0,
        )
        self.assertAlmostEqual(
            resultants[1],
            0.0,
            places=7,
        )

    def test_pure_axial_linear_elastic_response(self) -> None:
        axial_strain = 3.0e-4

        response = evaluate_linear_elastic_section(
            section=self.section,
            axial_strain=axial_strain,
            curvature=0.0,
            elastic_modulus=self.elastic_modulus,
        )

        expected_stress = self.elastic_modulus * axial_strain
        expected_force = (
            self.elastic_modulus
            * self.section.area
            * axial_strain
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
            rtol=0.0,
            atol=1.0e-5,
        )
        np.testing.assert_allclose(
            response.axial_force,
            expected_force,
            rtol=1.0e-14,
            atol=0.0,
        )
        self.assertAlmostEqual(
            response.bending_moment,
            0.0,
            places=6,
        )

    def test_pure_bending_linear_elastic_response(self) -> None:
        curvature = 2.0e-4

        response = evaluate_linear_elastic_section(
            section=self.section,
            axial_strain=0.0,
            curvature=curvature,
            elastic_modulus=self.elastic_modulus,
        )

        expected_moment = (
            self.elastic_modulus
            * self.section.second_moment_x
            * curvature
        )

        self.assertAlmostEqual(
            response.axial_force,
            0.0,
            places=6,
        )
        np.testing.assert_allclose(
            response.bending_moment,
            expected_moment,
            rtol=1.0e-14,
            atol=0.0,
        )

        upper = np.argmax(self.section.y_coordinates)
        lower = np.argmin(self.section.y_coordinates)

        self.assertLess(response.fiber_strains[upper], 0.0)
        self.assertGreater(response.fiber_strains[lower], 0.0)

    def test_combined_response_and_tangent(self) -> None:
        axial_strain = 1.5e-4
        curvature = 8.0e-5

        response = evaluate_linear_elastic_section(
            section=self.section,
            axial_strain=axial_strain,
            curvature=curvature,
            elastic_modulus=self.elastic_modulus,
        )

        expected_tangent = np.array(
            [
                [
                    self.elastic_modulus * self.section.area,
                    -self.elastic_modulus
                    * self.section.first_moment_y,
                ],
                [
                    -self.elastic_modulus
                    * self.section.first_moment_y,
                    self.elastic_modulus
                    * self.section.second_moment_x,
                ],
            ],
            dtype=np.float64,
        )
        generalized_strain = np.array(
            [axial_strain, curvature],
            dtype=np.float64,
        )

        np.testing.assert_allclose(
            response.tangent,
            expected_tangent,
            rtol=0.0,
            atol=1.0e-8,
        )
        np.testing.assert_allclose(
            response.resultants,
            expected_tangent @ generalized_strain,
            rtol=1.0e-14,
            atol=1.0e-6,
        )

    def test_invalid_geometry_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            create_annular_fiber_section(
                outer_diameter=0.0,
                thickness=0.027,
            )

        with self.assertRaises(ValueError):
            create_annular_fiber_section(
                outer_diameter=6.0,
                thickness=0.0,
            )

        with self.assertRaises(ValueError):
            create_annular_fiber_section(
                outer_diameter=6.0,
                thickness=3.0,
            )

    def test_invalid_discretization_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            create_annular_fiber_section(
                outer_diameter=6.0,
                thickness=0.027,
                n_circumferential=30,
            )

        with self.assertRaises(ValueError):
            create_annular_fiber_section(
                outer_diameter=6.0,
                thickness=0.027,
                n_radial=0,
            )

        with self.assertRaises(TypeError):
            create_annular_fiber_section(
                outer_diameter=6.0,
                thickness=0.027,
                n_circumferential=32.0,
            )

    def test_invalid_response_inputs_are_rejected(self) -> None:
        with self.assertRaises(ValueError):
            evaluate_linear_elastic_section(
                section=self.section,
                axial_strain=0.0,
                curvature=0.0,
                elastic_modulus=0.0,
            )

        with self.assertRaises(ValueError):
            integrate_fiber_stresses(
                section=self.section,
                fiber_stresses=np.zeros(
                    self.section.n_fibers - 1,
                    dtype=np.float64,
                ),
            )

        with self.assertRaises(ValueError):
            compute_fiber_strains(
                section=self.section,
                axial_strain=np.nan,
                curvature=0.0,
            )


if __name__ == "__main__":
    unittest.main()

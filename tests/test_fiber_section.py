# -*- coding: utf-8 -*-
"""Regression tests for circular annular fiber-section geometry."""

import unittest

import numpy as np

from fem.fiber_section import create_annular_fiber_section


class TestAnnularFiberSection(unittest.TestCase):
    """Verify annular-sector fiber geometry and section properties."""

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
        section = create_annular_fiber_section(
            outer_diameter=6.0,
            thickness=0.027,
            n_circumferential=32,
            n_radial=2,
        )

        np.testing.assert_allclose(
            section.area,
            section.exact_area,
            rtol=0.0,
            atol=1.0e-14,
        )
        self.assertAlmostEqual(
            section.centroid_x,
            0.0,
            places=14,
        )
        self.assertAlmostEqual(
            section.centroid_y,
            0.0,
            places=14,
        )

    def test_centroids_align_with_coordinate_axes(self) -> None:
        section = create_annular_fiber_section(
            outer_diameter=6.0,
            thickness=0.027,
            n_circumferential=32,
            n_radial=2,
        )

        first_layer = section.radial_indices == 0
        x = section.x_coordinates[first_layer]
        y = section.y_coordinates[first_layer]

        self.assertAlmostEqual(y[0], 0.0, places=14)
        self.assertGreater(x[0], 0.0)

        self.assertAlmostEqual(x[8], 0.0, places=14)
        self.assertGreater(y[8], 0.0)

        self.assertAlmostEqual(y[16], 0.0, places=14)
        self.assertLess(x[16], 0.0)

        self.assertAlmostEqual(x[24], 0.0, places=14)
        self.assertLess(y[24], 0.0)

    def test_second_moment_symmetry_and_accuracy(self) -> None:
        section = create_annular_fiber_section(
            outer_diameter=6.0,
            thickness=0.027,
            n_circumferential=32,
            n_radial=2,
        )

        np.testing.assert_allclose(
            section.second_moment_x,
            section.second_moment_y,
            rtol=0.0,
            atol=1.0e-13,
        )
        self.assertAlmostEqual(
            section.product_moment_xy,
            0.0,
            places=14,
        )

        relative_error = abs(
            section.second_moment_x
            - section.exact_second_moment
        ) / section.exact_second_moment

        self.assertLess(relative_error, 5.0e-3)

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


if __name__ == "__main__":
    unittest.main()

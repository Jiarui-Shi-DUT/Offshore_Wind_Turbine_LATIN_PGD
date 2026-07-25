# -*- coding: utf-8 -*-
"""Regression tests for the one-dimensional PGD basis representation."""

import unittest

import numpy as np

from latin.pgd_basis import PGDBasis1D, PGDMode1D


class TestPGDBasis1D(unittest.TestCase):
    """Verify empty and single-mode space-time reconstruction."""

    def test_empty_basis_returns_zero_fields(self) -> None:
        basis = PGDBasis1D(
            n_elements=90,
            n_time=2001,
        )

        self.assertEqual(basis.n_modes, 0)
        self.assertEqual(basis.field_shape, (2001, 90))

        for field in (
            basis.plastic_strain_correction(),
            basis.plastic_strain_rate_correction(),
            basis.stress_correction(),
        ):
            self.assertEqual(field.shape, (2001, 90))
            self.assertTrue(np.all(np.isfinite(field)))
            self.assertTrue(np.all(field == 0.0))

    def test_single_mode_reconstruction(self) -> None:
        spatial_plastic_strain = np.linspace(
            1.0,
            2.0,
            90,
            dtype=np.float64,
        )
        spatial_stress = -3.0 * spatial_plastic_strain
        temporal_amplitude = np.linspace(
            0.0,
            1.0,
            2001,
            dtype=np.float64,
        )
        temporal_rate = np.full(
            2001,
            0.1,
            dtype=np.float64,
        )

        mode = PGDMode1D(
            spatial_plastic_strain=spatial_plastic_strain,
            spatial_stress=spatial_stress,
            temporal_amplitude=temporal_amplitude,
            temporal_rate=temporal_rate,
            iteration_added=1,
        )
        basis = PGDBasis1D(
            n_elements=90,
            n_time=2001,
        )
        basis.append(mode)

        expected_strain = np.outer(
            temporal_amplitude,
            spatial_plastic_strain,
        )
        expected_rate = np.outer(
            temporal_rate,
            spatial_plastic_strain,
        )
        expected_stress = np.outer(
            temporal_amplitude,
            spatial_stress,
        )

        self.assertEqual(basis.n_modes, 1)
        np.testing.assert_allclose(
            basis.plastic_strain_correction(),
            expected_strain,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            basis.plastic_strain_rate_correction(),
            expected_rate,
            rtol=0.0,
            atol=0.0,
        )
        np.testing.assert_allclose(
            basis.stress_correction(),
            expected_stress,
            rtol=0.0,
            atol=0.0,
        )

        # PGDBasis1D stores a deep copy of an appended mode.
        mode.spatial_stress[:] = 999.0
        np.testing.assert_allclose(
            basis.stress_correction(),
            expected_stress,
            rtol=0.0,
            atol=0.0,
        )


if __name__ == "__main__":
    unittest.main()

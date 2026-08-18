# -*- coding: utf-8 -*-
"""Tests for immutable tower PGD mode and basis representations."""

import unittest

import numpy as np

from latin.pgd_basis import PGDBasisTower, PGDModeTower


class TestPGDModeTower(unittest.TestCase):
    """Verify tower PGD mode ownership and reconstruction."""

    def test_mode_owns_read_only_detached_arrays(self) -> None:
        p = np.array([1.0, 2.0, 3.0])
        s = np.array([-2.0, -4.0, -6.0])
        amplitude = np.array([0.0, 0.2, 0.4])
        rate = np.array([0.1, 0.2, 0.3])

        mode = PGDModeTower(
            spatial_plastic_strain=p,
            spatial_stress=s,
            temporal_amplitude=amplitude,
            temporal_rate=rate,
            iteration_added=2,
        )

        p[:] = 99.0
        s[:] = 99.0
        amplitude[:] = 99.0
        rate[:] = 99.0

        np.testing.assert_allclose(
            mode.spatial_plastic_strain,
            np.array([1.0, 2.0, 3.0]),
        )
        np.testing.assert_allclose(
            mode.spatial_stress,
            np.array([-2.0, -4.0, -6.0]),
        )
        self.assertFalse(
            mode.spatial_plastic_strain.flags.writeable
        )
        self.assertFalse(
            mode.temporal_amplitude.flags.writeable
        )

        with self.assertRaises(ValueError):
            mode.temporal_amplitude[0] = 1.0

    def test_mode_reconstructs_rank_one_fields(self) -> None:
        p = np.array([1.0, -0.5])
        s = np.array([-3.0, 1.5])
        amplitude = np.array([0.0, 0.2, -0.1])
        rate = np.array([0.4, 0.4, -0.6])

        mode = PGDModeTower(
            spatial_plastic_strain=p,
            spatial_stress=s,
            temporal_amplitude=amplitude,
            temporal_rate=rate,
        )

        np.testing.assert_allclose(
            mode.plastic_strain_correction(),
            np.outer(amplitude, p),
        )
        np.testing.assert_allclose(
            mode.plastic_strain_rate_correction(),
            np.outer(rate, p),
        )
        np.testing.assert_allclose(
            mode.stress_correction(),
            np.outer(amplitude, s),
        )


class TestPGDBasisTower(unittest.TestCase):
    """Verify immutable value-style tower PGD basis semantics."""

    def _mode(self, scale: float) -> PGDModeTower:
        return PGDModeTower(
            spatial_plastic_strain=scale
            * np.array([1.0, 2.0, 3.0]),
            spatial_stress=scale
            * np.array([-1.0, -2.0, -3.0]),
            temporal_amplitude=np.array(
                [0.0, 0.1, 0.2, 0.3]
            ),
            temporal_rate=np.array(
                [0.1, 0.1, 0.1, 0.1]
            ),
        )

    def test_empty_basis_has_zero_reconstructions(self) -> None:
        basis = PGDBasisTower(
            n_material_points=3,
            n_time=4,
        )

        self.assertEqual(basis.n_modes, 0)
        self.assertEqual(basis.field_shape, (4, 3))
        self.assertEqual(
            basis.spatial_plastic_strain_matrix().shape,
            (3, 0),
        )
        self.assertEqual(
            basis.temporal_amplitude_matrix().shape,
            (4, 0),
        )
        np.testing.assert_allclose(
            basis.plastic_strain_correction(),
            np.zeros((4, 3)),
        )
        np.testing.assert_allclose(
            basis.plastic_strain_rate_correction(),
            np.zeros((4, 3)),
        )
        np.testing.assert_allclose(
            basis.stress_correction(),
            np.zeros((4, 3)),
        )

    def test_with_appended_does_not_mutate_input_basis(self) -> None:
        empty = PGDBasisTower(
            n_material_points=3,
            n_time=4,
        )
        mode = self._mode(1.0)

        one_mode = empty.with_appended(mode)
        two_modes = one_mode.with_appended(
            self._mode(0.5)
        )

        self.assertEqual(empty.n_modes, 0)
        self.assertEqual(one_mode.n_modes, 1)
        self.assertEqual(two_modes.n_modes, 2)
        self.assertFalse(hasattr(one_mode, "append"))

        np.testing.assert_allclose(
            one_mode.plastic_strain_correction(),
            mode.plastic_strain_correction(),
        )

    def test_temporal_replacement_preserves_spatial_basis(self) -> None:
        basis = PGDBasisTower(
            n_material_points=3,
            n_time=4,
            modes=(
                self._mode(1.0),
                self._mode(0.5),
            ),
        )
        p_before = basis.spatial_plastic_strain_matrix()
        s_before = basis.spatial_stress_matrix()
        amplitude_before = basis.temporal_amplitude_matrix()

        new_amplitudes = np.array(
            [
                [0.0, 0.0],
                [0.2, 0.1],
                [0.4, -0.1],
                [0.3, 0.2],
            ]
        )
        new_rates = np.array(
            [
                [0.8, 0.4],
                [0.8, 0.4],
                [0.8, -0.8],
                [-0.4, 1.2],
            ]
        )

        updated = basis.with_temporal_coordinates(
            temporal_amplitudes=new_amplitudes,
            temporal_rates=new_rates,
        )

        np.testing.assert_allclose(
            updated.spatial_plastic_strain_matrix(),
            p_before,
        )
        np.testing.assert_allclose(
            updated.spatial_stress_matrix(),
            s_before,
        )
        np.testing.assert_allclose(
            updated.temporal_amplitude_matrix(),
            new_amplitudes,
        )
        np.testing.assert_allclose(
            updated.temporal_rate_matrix(),
            new_rates,
        )
        np.testing.assert_allclose(
            basis.temporal_amplitude_matrix(),
            amplitude_before,
        )


if __name__ == "__main__":
    unittest.main()

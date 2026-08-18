# -*- coding: utf-8 -*-
"""Unit tests for tower LATIN-PGD relaxation and Eqs. (60), (76)-(77)."""

import unittest

import numpy as np

from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import MaterialPointMetric
from latin.tower_iteration_control import (
    evaluate_tower_trial,
    relax_tower_global_state,
    tower_plastic_state_norm_squared,
    tower_relative_latin_indicator,
)
from latin.tower_state import LatinStateTower


def _state(
    time,
    *,
    stress,
    beta,
    R_bar,
    plastic_strain_rate,
    elastic_strain,
    alpha_rate,
    r_bar_rate,
    damage_rate=None,
    energy_release_rate=None,
    plastic_strain=None,
    alpha=None,
    r_bar=None,
    damage=None,
):
    shape = np.asarray(stress, dtype=np.float64).shape

    def field(values, default=0.0):
        if values is None:
            return np.full(shape, default, dtype=np.float64)
        return np.asarray(values, dtype=np.float64)

    return LatinStateTower(
        time=np.asarray(time, dtype=np.float64),
        plastic_strain_rate=field(plastic_strain_rate),
        elastic_strain=field(elastic_strain),
        alpha_rate=field(alpha_rate),
        r_bar_rate=field(r_bar_rate),
        damage_rate=field(damage_rate),
        stress=field(stress),
        beta=field(beta),
        R_bar=field(R_bar),
        energy_release_rate=field(energy_release_rate),
        plastic_strain=field(plastic_strain),
        alpha=field(alpha),
        r_bar=field(r_bar),
        damage=field(damage),
    )


def _directions(shape):
    rows, cols = shape
    q = np.arange(cols, dtype=np.float64)[None, :]
    t = np.arange(rows, dtype=np.float64)[:, None]
    return DescentSearchDirections(
        H_sigma=2.0e-6 + 1.0e-7 * t + 1.0e-8 * q,
        H_beta=3.0e-4 + 2.0e-5 * t + 1.0e-6 * q,
        H_R_bar=4.0e-3 + 3.0e-4 * t + 1.0e-5 * q,
        b_damage=np.zeros(shape, dtype=np.float64),
        regularization=0.15,
    )


class TestTowerIterationControl(unittest.TestCase):
    """Check the frozen Stage-6 iteration-control contract."""

    def test_relaxation_blends_all_fields_without_reprojecting_y(self) -> None:
        time = np.array([0.0, 0.4, 1.0], dtype=np.float64)
        shape = (time.size, 2)

        previous_fields = {}
        candidate_fields = {}
        for index, name in enumerate(
            LatinStateTower.MATERIAL_FIELD_NAMES
        ):
            previous_fields[name] = np.full(
                shape,
                0.01 * index,
                dtype=np.float64,
            )
            candidate_fields[name] = np.full(
                shape,
                0.20 + 0.02 * index,
                dtype=np.float64,
            )

        previous_fields["damage"][:] = 0.10
        candidate_fields["damage"][:] = 0.30
        previous_fields["stress"][:] = 20.0
        candidate_fields["stress"][:] = 60.0
        previous_fields["energy_release_rate"][:] = 5.0
        candidate_fields["energy_release_rate"][:] = 9.0

        previous = LatinStateTower(time=time, **previous_fields)
        candidate = LatinStateTower(time=time, **candidate_fields)

        relaxed = relax_tower_global_state(
            previous,
            candidate,
            relaxation=0.25,
        )

        for name in LatinStateTower.MATERIAL_FIELD_NAMES:
            expected = (
                0.75 * getattr(previous, name)
                + 0.25 * getattr(candidate, name)
            )
            np.testing.assert_allclose(
                getattr(relaxed, name),
                expected,
                rtol=0.0,
                atol=0.0,
            )

        expected_y = 0.75 * 5.0 + 0.25 * 9.0
        self.assertAlmostEqual(
            float(relaxed.energy_release_rate[1, 0]),
            expected_y,
        )

        sigma_relaxed = float(relaxed.stress[1, 0])
        damage_relaxed = float(relaxed.damage[1, 0])
        constitutive_y = (
            sigma_relaxed**2
            / (2.0 * 210_000.0 * (1.0 - damage_relaxed) ** 2)
        )
        self.assertGreater(
            abs(expected_y - constitutive_y),
            1.0,
        )

        np.testing.assert_allclose(
            previous.energy_release_rate,
            np.full(shape, 5.0),
        )
        np.testing.assert_allclose(
            candidate.energy_release_rate,
            np.full(shape, 9.0),
        )

    def test_eq77_norm_matches_manual_metric_trapezoid(self) -> None:
        time = np.array([0.0, 0.5, 1.5], dtype=np.float64)
        stress = np.array(
            [[1.0, 2.0], [2.0, 1.0], [3.0, 4.0]],
            dtype=np.float64,
        )
        beta = 0.5 * stress
        R_bar = 0.25 * stress
        epsp_rate = 0.02 * stress
        eps_e = 0.001 * stress
        alpha_rate = 0.03 * stress
        rbar_rate = 0.04 * stress

        state = _state(
            time,
            stress=stress,
            beta=beta,
            R_bar=R_bar,
            plastic_strain_rate=epsp_rate,
            elastic_strain=eps_e,
            alpha_rate=alpha_rate,
            r_bar_rate=rbar_rate,
        )
        directions = _directions(state.field_shape)
        metric = MaterialPointMetric(
            np.array([2.0, 3.0], dtype=np.float64)
        )
        C0 = np.array([100.0, 250.0], dtype=np.float64)

        density = (
            stress**2 * directions.H_sigma
            + beta**2 * directions.H_beta
            + R_bar**2 * directions.H_R_bar
            + epsp_rate**2 / directions.H_sigma
            + eps_e**2 * C0[None, :]
            + alpha_rate**2 / directions.H_beta
            + rbar_rate**2 / directions.H_R_bar
        )
        space_integral = density @ metric.weights
        expected = float(
            np.trapz(space_integral, x=time)
        )

        actual = tower_plastic_state_norm_squared(
            state=state,
            directions=directions,
            metric=metric,
            reference_modulus=C0,
        )
        self.assertAlmostEqual(
            actual,
            expected,
            places=12,
        )

    def test_indicator_ignores_damage_and_support_histories_directly(self) -> None:
        time = np.array([0.0, 0.4, 1.0], dtype=np.float64)
        shape = (time.size, 2)
        ones = np.ones(shape, dtype=np.float64)

        local = _state(
            time,
            stress=10.0 * ones,
            beta=2.0 * ones,
            R_bar=3.0 * ones,
            plastic_strain_rate=0.1 * ones,
            elastic_strain=0.001 * ones,
            alpha_rate=0.2 * ones,
            r_bar_rate=0.3 * ones,
            damage_rate=0.01 * ones,
            energy_release_rate=5.0 * ones,
            plastic_strain=0.4 * ones,
            alpha=0.5 * ones,
            r_bar=0.6 * ones,
            damage=0.1 * ones,
        )
        global_state = _state(
            time,
            stress=10.0 * ones,
            beta=2.0 * ones,
            R_bar=3.0 * ones,
            plastic_strain_rate=0.1 * ones,
            elastic_strain=0.001 * ones,
            alpha_rate=0.2 * ones,
            r_bar_rate=0.3 * ones,
            damage_rate=0.07 * ones,
            energy_release_rate=20.0 * ones,
            plastic_strain=0.9 * ones,
            alpha=1.1 * ones,
            r_bar=1.2 * ones,
            damage=0.7 * ones,
        )
        directions = _directions(shape)
        metric = MaterialPointMetric(
            np.array([1.2, 0.8], dtype=np.float64)
        )

        indicator = tower_relative_latin_indicator(
            local_state=local,
            global_state=global_state,
            directions=directions,
            metric=metric,
            reference_modulus=210_000.0,
        )
        self.assertEqual(indicator, 0.0)

    def test_indicator_has_expected_one_third_scaling_value(self) -> None:
        time = np.array([0.0, 0.3, 0.9], dtype=np.float64)
        shape = (time.size, 3)
        base = (
            1.0
            + np.arange(time.size, dtype=np.float64)[:, None]
            + 0.2 * np.arange(3, dtype=np.float64)[None, :]
        )

        global_state = _state(
            time,
            stress=base,
            beta=0.5 * base,
            R_bar=0.4 * base,
            plastic_strain_rate=0.01 * base,
            elastic_strain=0.001 * base,
            alpha_rate=0.02 * base,
            r_bar_rate=0.03 * base,
        )
        local = _state(
            time,
            stress=2.0 * base,
            beta=1.0 * base,
            R_bar=0.8 * base,
            plastic_strain_rate=0.02 * base,
            elastic_strain=0.002 * base,
            alpha_rate=0.04 * base,
            r_bar_rate=0.06 * base,
        )
        directions = _directions(shape)
        metric = MaterialPointMetric(
            np.array([0.7, 1.1, 1.6], dtype=np.float64)
        )

        indicator = tower_relative_latin_indicator(
            local_state=local,
            global_state=global_state,
            directions=directions,
            metric=metric,
            reference_modulus=np.array(
                [100.0, 150.0, 220.0],
                dtype=np.float64,
            ),
        )
        self.assertAlmostEqual(
            indicator,
            1.0 / 3.0,
            places=12,
        )

    def test_trial_evaluation_uses_relaxed_state_and_frozen_previous_xi(self) -> None:
        time = np.array([0.0, 0.5, 1.0], dtype=np.float64)
        shape = (time.size, 2)
        base = (
            1.0
            + np.arange(time.size, dtype=np.float64)[:, None]
            + 0.1 * np.arange(2, dtype=np.float64)[None, :]
        )
        zeros = np.zeros(shape, dtype=np.float64)

        baseline = _state(
            time,
            stress=zeros,
            beta=zeros,
            R_bar=zeros,
            plastic_strain_rate=zeros,
            elastic_strain=zeros,
            alpha_rate=zeros,
            r_bar_rate=zeros,
        )
        local = _state(
            time,
            stress=base,
            beta=0.5 * base,
            R_bar=0.4 * base,
            plastic_strain_rate=0.01 * base,
            elastic_strain=0.001 * base,
            alpha_rate=0.02 * base,
            r_bar_rate=0.03 * base,
            damage=0.1 * np.ones(shape),
        )
        unrelaxed = local.copy()
        directions = _directions(shape)
        metric = MaterialPointMetric(
            np.array([1.0, 2.0], dtype=np.float64)
        )

        evaluation = evaluate_tower_trial(
            baseline_state=baseline,
            local_state=local,
            unrelaxed_state=unrelaxed,
            directions=directions,
            metric=metric,
            reference_modulus=210_000.0,
            previous_indicator=0.6,
            relaxation=0.5,
            latin_tolerance=0.34,
        )

        self.assertAlmostEqual(
            evaluation.indicator,
            1.0 / 3.0,
            places=12,
        )
        expected_zeta = (
            (0.6 - 1.0 / 3.0)
            / (0.6 + 1.0 / 3.0)
        )
        self.assertAlmostEqual(
            evaluation.saturation,
            expected_zeta,
            places=12,
        )
        self.assertTrue(evaluation.converged)
        self.assertTrue(evaluation.finite)
        self.assertEqual(
            evaluation.previous_indicator,
            0.6,
        )
        np.testing.assert_allclose(
            evaluation.relaxed_state.stress,
            0.5 * local.stress,
            rtol=0.0,
            atol=0.0,
        )


if __name__ == "__main__":
    unittest.main()

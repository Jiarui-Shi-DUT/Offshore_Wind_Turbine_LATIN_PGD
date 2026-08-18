# -*- coding: utf-8 -*-
"""Unit tests for tower LATIN-PGD global-stage finishing."""

import unittest

import numpy as np

from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_global_finishing import (
    build_unrelaxed_candidate,
    prepare_frozen_global_data,
)
from latin.tower_pgd_time_update import update_tower_pgd_time_functions
from latin.tower_state import LatinStateTower, MaterialPointLayout
from material.viscoplastic_damage_1d import MaterialParameters


def _operator():
    layout = MaterialPointLayout(
        n_elements=1,
        n_gauss=2,
        n_fibers=2,
    )
    metric = MaterialPointMetric(
        np.array([1.0, 2.0, 1.5, 0.5], dtype=np.float64)
    )
    compatibility = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
            [2.0, -1.0],
        ],
        dtype=np.float64,
    )
    operator = TowerEquilibriumOperator(
        layout=layout,
        metric=metric,
        reference_modulus=np.full(4, 200.0, dtype=np.float64),
        compatibility_matrix=compatibility,
        free_dofs=np.array([1, 3]),
        n_dof=4,
    )
    return metric, operator


def _state(time, *, stress_shift=0.0, local=False):
    nt = time.size
    nq = 4
    q = np.arange(nq, dtype=np.float64)[None, :]
    t = np.arange(nt, dtype=np.float64)[:, None]

    sign = np.array([1.0, -1.0, 1.0, -1.0], dtype=np.float64)[None, :]
    stress = sign * (10.0 + 1.5 * t + 0.8 * q) + stress_shift
    elastic = stress / 200.0
    plastic_rate = 0.002 * t + 0.0003 * q
    plastic = 0.0005 * t + 0.0001 * q
    alpha = 0.003 * t + 0.0002 * q
    r_bar = 0.002 * t + 0.00015 * q
    alpha_rate = 0.01 + 0.001 * t + 0.0002 * q
    r_bar_rate = 0.008 + 0.0008 * t + 0.0001 * q
    beta = 50.0 * alpha
    R_bar = 30.0 * r_bar
    damage = 0.01 * t + 0.001 * q
    damage_rate = 0.004 + 0.0002 * t + 0.00005 * q
    energy = np.zeros((nt, nq), dtype=np.float64)

    if local:
        # Mimic the local-stage ascent semantics: forces remain fixed while
        # integrated histories/rates and damaged elastic strain may change.
        plastic_rate = plastic_rate + 0.003 + 0.0002 * q
        plastic = plastic + 0.0007 * t
        alpha = alpha + 0.0006 * t
        r_bar = r_bar + 0.0004 * t
        alpha_rate = alpha_rate + 0.002 + 0.0001 * q
        r_bar_rate = r_bar_rate + 0.0015 + 0.0001 * q
        # Keep local stress equal to baseline stress for the project ascent
        # choice and create the nonlinear damaged-elastic mismatch via eps_e.
        elastic = elastic * (1.0 + 0.08 * t)
        damage = damage + 0.005 * t
        damage_rate = damage_rate + 0.001 * t

    return LatinStateTower(
        time=time,
        plastic_strain_rate=plastic_rate,
        elastic_strain=elastic,
        alpha_rate=alpha_rate,
        r_bar_rate=r_bar_rate,
        damage_rate=damage_rate,
        stress=stress,
        beta=beta,
        R_bar=R_bar,
        energy_release_rate=energy,
        plastic_strain=plastic,
        alpha=alpha,
        r_bar=r_bar,
        damage=damage,
    )


def _directions(shape):
    nt, nq = shape
    t = np.arange(nt, dtype=np.float64)[:, None]
    q = np.arange(nq, dtype=np.float64)[None, :]
    return DescentSearchDirections(
        H_sigma=0.01 + 0.001 * t + 0.0002 * q,
        H_beta=0.02 + 0.001 * t + 0.0001 * q,
        H_R_bar=0.03 + 0.001 * t + 0.0001 * q,
        b_damage=np.zeros(shape, dtype=np.float64),
        regularization=0.15,
    )


def _fixed_basis_result(time, forcing, directions, metric, operator):
    p = np.array([0.7, -0.4, 0.5, -0.2], dtype=np.float64)
    p /= metric.norm(p)
    s = operator.apply_spatial(p).stress
    basis = PGDBasisTower(
        n_material_points=operator.n_material_points,
        n_time=time.size,
        modes=(
            PGDModeTower(
                p,
                s,
                np.zeros(time.size),
                np.zeros(time.size),
            ),
        ),
    )
    return update_tower_pgd_time_functions(
        basis=basis,
        time=time,
        forcing=forcing,
        H_sigma=directions.H_sigma,
        metric=metric,
        equilibrium_operator=operator,
    )


class TestPrepareFrozenGlobalData(unittest.TestCase):
    """Verify the paper-separated plastic/damage branch split."""

    def test_preparation_uses_paper_separated_forcing(self):
        metric, operator = _operator()
        time = np.array([0.0, 0.2, 0.5, 0.9], dtype=np.float64)
        baseline = _state(time)
        local = _state(time, local=True)
        directions = _directions(baseline.field_shape)

        frozen = prepare_frozen_global_data(
            baseline_state=baseline,
            local_state=local,
            directions=directions,
            equilibrium_operator=operator,
        )

        expected_forcing = (
            local.plastic_strain_rate
            - baseline.plastic_strain_rate
            - directions.H_sigma * (local.stress - baseline.stress)
        )
        expected_residual_strain = (
            baseline.stress
            - local.stress
            - operator.reference_modulus[None, :]
            * (baseline.elastic_strain - local.elastic_strain)
        ) / operator.reference_modulus[None, :]

        np.testing.assert_allclose(
            frozen.full_plastic_forcing,
            expected_forcing,
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            frozen.damage_residual_strain,
            expected_residual_strain,
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            frozen.damage_projection.source_strain,
            expected_residual_strain,
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            operator.equilibrium_residual(
                frozen.damage_projection.stress
            ),
            np.zeros((time.size, operator.n_free_dofs)),
            rtol=0.0,
            atol=1.0e-12,
        )

        # This is deliberately *not* the current 1D coupled forcing.
        coupled_1d_variant = (
            expected_forcing
            + directions.H_sigma * frozen.damage_projection.stress
        )
        self.assertGreater(
            np.linalg.norm(
                frozen.full_plastic_forcing - coupled_1d_variant
            ),
            0.0,
        )


class TestBuildUnrelaxedCandidate(unittest.TestCase):
    """Verify Eq. (73)-(75), mechanics, damage inheritance, and final Y."""

    def test_complete_candidate_obeys_frozen_contract(self):
        metric, operator = _operator()
        time = np.array([0.0, 0.2, 0.5, 0.9], dtype=np.float64)
        baseline = _state(time)
        local = _state(time, local=True)
        directions = _directions(baseline.field_shape)
        material = MaterialParameters(
            E=200.0,
            C=50.0,
            R_inf=30.0,
            h=0.25,
        )

        frozen = prepare_frozen_global_data(
            baseline_state=baseline,
            local_state=local,
            directions=directions,
            equilibrium_operator=operator,
        )
        fixed = _fixed_basis_result(
            time=time,
            forcing=frozen.full_plastic_forcing,
            directions=directions,
            metric=metric,
            operator=operator,
        )

        baseline_before = baseline.copy()
        local_before = local.copy()

        candidate = build_unrelaxed_candidate(
            baseline_state=baseline,
            local_state=local,
            directions=directions,
            frozen_data=frozen,
            fixed_basis_result=fixed,
            materials=material,
        )
        state = candidate.state

        np.testing.assert_allclose(
            state.plastic_strain,
            baseline.plastic_strain + fixed.plastic_strain_correction,
        )
        np.testing.assert_allclose(
            state.plastic_strain_rate,
            baseline.plastic_strain_rate
            + fixed.plastic_strain_rate_correction,
        )
        np.testing.assert_allclose(
            state.stress,
            baseline.stress
            + fixed.plastic_projection.stress
            + frozen.damage_projection.stress,
        )
        np.testing.assert_allclose(
            state.elastic_strain,
            baseline.elastic_strain
            + fixed.plastic_projection.compatible_strain
            - fixed.plastic_strain_correction
            + frozen.damage_projection.compatible_strain,
        )

        # Eq. (73) is exact in every candidate row.
        np.testing.assert_allclose(state.beta, material.C * state.alpha)
        np.testing.assert_allclose(state.R_bar, material.R_inf * state.r_bar)

        # Eq. (74) descent relation is exact at t0 and the BE nodes.
        np.testing.assert_allclose(
            state.alpha_rate
            + directions.H_beta * state.beta,
            local.alpha_rate + directions.H_beta * local.beta,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            state.r_bar_rate
            + directions.H_R_bar * state.R_bar,
            local.r_bar_rate + directions.H_R_bar * local.R_bar,
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            state.alpha_rate[1:],
            np.diff(state.alpha, axis=0) / np.diff(time)[:, None],
            rtol=0.0,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            state.r_bar_rate[1:],
            np.diff(state.r_bar, axis=0) / np.diff(time)[:, None],
            rtol=0.0,
            atol=1.0e-12,
        )

        # Eq. (75) plus the frozen tower-v1 direct-copy history choice.
        np.testing.assert_array_equal(state.damage_rate, local.damage_rate)
        np.testing.assert_array_equal(state.damage, local.damage)

        tensile = state.stress >= 0.0
        expected_y = np.where(
            tensile,
            state.stress**2 / (2.0 * material.E * (1.0 - state.damage) ** 2),
            material.h
            * state.stress**2
            / (2.0 * material.E * (1.0 - material.h * state.damage) ** 2),
        )
        np.testing.assert_allclose(
            state.energy_release_rate,
            expected_y,
            rtol=0.0,
            atol=1.0e-12,
        )

        np.testing.assert_allclose(
            candidate.total_displacement_correction,
            candidate.plastic_displacement_correction
            + candidate.damage_displacement_correction,
        )
        np.testing.assert_allclose(
            operator.equilibrium_residual(
                fixed.plastic_projection.stress
                + frozen.damage_projection.stress
            ),
            np.zeros((time.size, operator.n_free_dofs)),
            rtol=0.0,
            atol=1.0e-12,
        )

        # Pure-function transaction semantics: inputs remain byte-for-byte equal.
        for field_name in LatinStateTower.MATERIAL_FIELD_NAMES:
            np.testing.assert_array_equal(
                getattr(baseline, field_name),
                getattr(baseline_before, field_name),
            )
            np.testing.assert_array_equal(
                getattr(local, field_name),
                getattr(local_before, field_name),
            )

    def test_trial_a_and_b_share_hardening_and_damage_but_y_tracks_stress(self):
        metric, operator = _operator()
        time = np.array([0.0, 0.25, 0.55, 0.95], dtype=np.float64)
        baseline = _state(time)
        local = _state(time, local=True)
        directions = _directions(baseline.field_shape)
        material = MaterialParameters(E=200.0, C=50.0, R_inf=30.0, h=0.25)

        frozen = prepare_frozen_global_data(
            baseline_state=baseline,
            local_state=local,
            directions=directions,
            equilibrium_operator=operator,
        )

        fixed_a = _fixed_basis_result(
            time,
            frozen.full_plastic_forcing,
            directions,
            metric,
            operator,
        )

        # Same frozen local/search/damage data and same forcing, but a
        # different reduced spatial basis produces a different mechanical
        # Trial-B-like correction.  This is the actual A/B data contract.
        empty_basis = PGDBasisTower(
            n_material_points=operator.n_material_points,
            n_time=time.size,
        )
        fixed_b = update_tower_pgd_time_functions(
            basis=empty_basis,
            time=time,
            forcing=frozen.full_plastic_forcing,
            H_sigma=directions.H_sigma,
            metric=metric,
            equilibrium_operator=operator,
        )

        trial_a = build_unrelaxed_candidate(
            baseline,
            local,
            directions,
            frozen,
            fixed_a,
            material,
        )
        trial_b = build_unrelaxed_candidate(
            baseline,
            local,
            directions,
            frozen,
            fixed_b,
            material,
        )

        np.testing.assert_allclose(trial_a.state.alpha, trial_b.state.alpha)
        np.testing.assert_allclose(trial_a.state.alpha_rate, trial_b.state.alpha_rate)
        np.testing.assert_allclose(trial_a.state.beta, trial_b.state.beta)
        np.testing.assert_allclose(trial_a.state.r_bar, trial_b.state.r_bar)
        np.testing.assert_allclose(trial_a.state.r_bar_rate, trial_b.state.r_bar_rate)
        np.testing.assert_allclose(trial_a.state.R_bar, trial_b.state.R_bar)
        np.testing.assert_array_equal(trial_a.state.damage, trial_b.state.damage)
        np.testing.assert_array_equal(trial_a.state.damage_rate, trial_b.state.damage_rate)
        self.assertGreater(
            np.linalg.norm(trial_a.state.stress - trial_b.state.stress),
            0.0,
        )
        self.assertGreater(
            np.linalg.norm(
                trial_a.state.energy_release_rate
                - trial_b.state.energy_release_rate
            ),
            0.0,
        )

    def test_mismatched_reduced_forcing_is_rejected(self):
        metric, operator = _operator()
        time = np.array([0.0, 0.2, 0.5, 0.9], dtype=np.float64)
        baseline = _state(time)
        local = _state(time, local=True)
        directions = _directions(baseline.field_shape)
        frozen = prepare_frozen_global_data(
            baseline,
            local,
            directions,
            operator,
        )
        wrong_fixed = _fixed_basis_result(
            time,
            1.2 * frozen.full_plastic_forcing,
            directions,
            metric,
            operator,
        )

        with self.assertRaisesRegex(
            ValueError,
            "FrozenGlobalData full plastic forcing",
        ):
            build_unrelaxed_candidate(
                baseline,
                local,
                directions,
                frozen,
                wrong_fixed,
                MaterialParameters(E=200.0, C=50.0, R_inf=30.0),
            )


if __name__ == "__main__":
    unittest.main()

# -*- coding: utf-8 -*-
"""Transaction-focused unit tests for the tower LATIN-PGD solver."""

import unittest
from unittest.mock import patch

import numpy as np

from latin.pgd_basis import PGDBasisTower, PGDModeTower
from latin.pgd_saturation import saturation_indicator
from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_global_finishing import TowerGlobalCandidate
from latin.tower_iteration_control import TowerTrialEvaluation
from latin.tower_latin_pgd_solver import (
    TowerLatinPGDTerminationReason,
    solve_tower_latin_pgd,
)
from latin.tower_pgd_enrichment import TowerEnrichmentResult
from latin.tower_pgd_time_update import FixedBasisPGDResult
from latin.tower_state import LatinStateTower, MaterialPointLayout
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerLatinPGDSolver(unittest.TestCase):
    """Freeze atomic commit, rollback, bootstrap, and common-baseline rules."""

    def setUp(self) -> None:
        self.time = np.array([0.0, 0.5, 1.0], dtype=np.float64)
        self.layout = MaterialPointLayout(1, 1, 2)
        self.metric = MaterialPointMetric(np.array([1.0, 1.5]))
        self.material = MaterialParameters(
            E=210_000.0,
            C=5_500.0,
            R_inf=30.0,
            h=0.2,
        )
        self.operator = TowerEquilibriumOperator(
            layout=self.layout,
            metric=self.metric,
            reference_modulus=np.full(2, self.material.E),
            compatibility_matrix=np.array([[1.0], [0.8]]),
            free_dofs=np.array([0], dtype=np.int64),
            n_dof=1,
        )
        self.initial = self._state(0.0)
        self.local = self._state(0.25)
        shape = self.initial.field_shape
        self.directions = DescentSearchDirections(
            H_sigma=np.full(shape, 2.0e-6),
            H_beta=np.full(shape, 2.0e-4),
            H_R_bar=np.full(shape, 3.0e-4),
            b_damage=np.zeros(shape),
            regularization=0.15,
        )
        self.frozen = type(
            "FrozenStub",
            (),
            {"full_plastic_forcing": np.ones(shape, dtype=np.float64)},
        )()

    def _state(self, stress_value: float) -> LatinStateTower:
        shape = (self.time.size, 2)
        zero = np.zeros(shape, dtype=np.float64)
        stress = np.full(shape, stress_value, dtype=np.float64)
        return LatinStateTower(
            time=self.time,
            plastic_strain_rate=zero,
            elastic_strain=stress / 210_000.0,
            alpha_rate=zero,
            r_bar_rate=zero,
            damage_rate=zero,
            stress=stress,
            beta=zero,
            R_bar=zero,
            energy_release_rate=zero,
            plastic_strain=zero,
            alpha=zero,
            r_bar=zero,
            damage=zero,
        )

    def _basis(self, amplitude: float, n_modes: int = 1) -> PGDBasisTower:
        if n_modes == 0:
            return PGDBasisTower(
                n_material_points=2,
                n_time=self.time.size,
            )
        modes = []
        for index in range(n_modes):
            p = np.array([1.0, -0.4 - 0.1 * index], dtype=np.float64)
            p = p / self.metric.norm(p)
            s = self.operator.apply_spatial(p).stress
            modes.append(
                PGDModeTower(
                    spatial_plastic_strain=p,
                    spatial_stress=s,
                    temporal_amplitude=np.full(
                        self.time.size,
                        amplitude + 0.1 * index,
                    ),
                    temporal_rate=np.zeros(self.time.size),
                    iteration_added=index,
                )
            )
        return PGDBasisTower(
            n_material_points=2,
            n_time=self.time.size,
            modes=tuple(modes),
        )

    def _fixed(
        self,
        basis: PGDBasisTower,
        relative_residual: float,
    ) -> FixedBasisPGDResult:
        shape = basis.field_shape
        zero = np.zeros(shape, dtype=np.float64)
        projection = self.operator.apply_history(zero)
        residual = np.full(shape, relative_residual, dtype=np.float64)
        return FixedBasisPGDResult(
            basis=basis,
            plastic_strain_correction=zero,
            plastic_strain_rate_correction=zero,
            plastic_projection=projection,
            mechanical_residual=residual,
            weighted_residual_norm=max(0.0, float(relative_residual)),
            relative_residual=max(0.0, float(relative_residual)),
            forcing_norm=1.0,
            reduced_converged=relative_residual <= 1.0e-4,
            condition_history=np.zeros(self.time.size),
            time_functions_updated=basis.n_modes > 0,
        )

    def _candidate(self, state: LatinStateTower) -> TowerGlobalCandidate:
        displacement = np.zeros((self.time.size, 1), dtype=np.float64)
        return TowerGlobalCandidate(
            state=state,
            total_displacement_correction=displacement,
            plastic_displacement_correction=displacement,
            damage_displacement_correction=displacement,
        )

    def _trial(
        self,
        unrelaxed: LatinStateTower,
        relaxed: LatinStateTower,
        indicator: float,
        previous_indicator: float,
        tolerance: float = 1.0e-4,
    ) -> TowerTrialEvaluation:
        return TowerTrialEvaluation(
            unrelaxed_state=unrelaxed,
            relaxed_state=relaxed,
            indicator=indicator,
            saturation=saturation_indicator(
                previous_indicator,
                indicator,
            ),
            converged=indicator <= tolerance,
            difference_norm=indicator,
            local_norm=1.0,
            global_norm=1.0,
            relaxation=0.8,
            latin_tolerance=tolerance,
            previous_indicator=previous_indicator,
        )

    def _accepted_enrichment(
        self,
        fixed_b: FixedBasisPGDResult,
    ) -> TowerEnrichmentResult:
        return TowerEnrichmentResult(
            accepted=True,
            failure_reason=None,
            candidate_fixed_basis_result=fixed_b,
            raw_mode=None,
            fixed_point_history=np.array([1.0e-7]),
            fixed_point_iterations=1,
            fixed_point_converged=True,
            projection_coefficients=np.zeros(
                max(0, fixed_b.basis.n_modes - 1)
            ),
            orthogonal_scale=1.0,
            spatial_novelty=1.0,
            temporal_significance=1.0,
            orthogonality_error=0.0,
            plastic_field_invariance_error=0.0,
            plastic_rate_field_invariance_error=0.0,
            stress_field_invariance_error=0.0,
            residual_norm_before=1.0,
            residual_norm_after=0.5,
            residual_benefit=0.5,
        )

    def _rejected_enrichment(self) -> TowerEnrichmentResult:
        return TowerEnrichmentResult(
            accepted=False,
            failure_reason="full_residual_benefit_insufficient",
            candidate_fixed_basis_result=None,
            raw_mode=None,
            fixed_point_history=np.array([1.0e-7]),
            fixed_point_iterations=1,
            fixed_point_converged=True,
            projection_coefficients=np.zeros(1),
            orthogonal_scale=1.0,
            spatial_novelty=1.0,
            temporal_significance=1.0,
            orthogonality_error=0.0,
            plastic_field_invariance_error=0.0,
            plastic_rate_field_invariance_error=0.0,
            stress_field_invariance_error=0.0,
            residual_norm_before=1.0,
            residual_norm_after=0.999,
            residual_benefit=0.001,
        )

    def _solve(self, **kwargs):
        controls = dict(
            initial_state=self.initial,
            materials=self.material,
            metric=self.metric,
            equilibrium_operator=self.operator,
            mode_significance_tolerance=0.0,
            acceptance_tolerance=0.0,
        )
        controls.update(kwargs)
        return solve_tower_latin_pgd(**controls)

    def test_trial_a_absolute_convergence_commits_updated_temporal_basis(self) -> None:
        persistent = self._basis(0.0)
        fixed_a = self._fixed(self._basis(0.7), 0.8)
        candidate_a = self._candidate(self._state(0.6))
        relaxed_a = self._state(0.5)
        trial_a = self._trial(candidate_a.state, relaxed_a, 5.0e-5, 1.0)

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            return_value=self.local,
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            return_value=self.directions,
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            return_value=self.frozen,
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            return_value=fixed_a,
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            return_value=candidate_a,
        ), patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            return_value=trial_a,
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once"
        ) as enrich_mock:
            result = self._solve(initial_basis=persistent)

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.CONVERGED,
        )
        self.assertEqual(result.commit_kind_history, ("A",))
        self.assertEqual(result.basis.n_modes, persistent.n_modes)
        np.testing.assert_allclose(
            result.basis.temporal_amplitude_matrix(),
            fixed_a.basis.temporal_amplitude_matrix(),
        )
        self.assertFalse(
            np.allclose(
                result.basis.temporal_amplitude_matrix(),
                persistent.temporal_amplitude_matrix(),
            )
        )
        np.testing.assert_allclose(result.state.stress, relaxed_a.stress)
        enrich_mock.assert_not_called()

    def test_reduced_residual_escape_commits_trial_a_without_enrichment(self) -> None:
        persistent = self._basis(0.0)
        fixed_a = self._fixed(self._basis(0.4), 1.0e-6)
        candidate_a = self._candidate(self._state(0.6))
        relaxed_a = self._state(0.5)
        trial_a = self._trial(candidate_a.state, relaxed_a, 0.90, 1.0)

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            return_value=self.local,
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            return_value=self.directions,
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            return_value=self.frozen,
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            return_value=fixed_a,
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            return_value=candidate_a,
        ), patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            return_value=trial_a,
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once"
        ) as enrich_mock:
            result = self._solve(
                initial_basis=persistent,
                max_iterations=1,
            )

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.MAX_ITERATIONS,
        )
        self.assertEqual(result.commit_kind_history, ("A",))
        self.assertAlmostEqual(result.final_indicator, 0.90)
        enrich_mock.assert_not_called()

    def test_enrichment_rejection_rolls_back_entire_iteration(self) -> None:
        persistent = self._basis(0.0)
        fixed_a = self._fixed(self._basis(0.6), 0.8)
        candidate_a = self._candidate(self._state(0.6))
        trial_a = self._trial(candidate_a.state, self._state(0.5), 0.90, 1.0)
        rejected = self._rejected_enrichment()

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            return_value=self.local,
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            return_value=self.directions,
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            return_value=self.frozen,
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            return_value=fixed_a,
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            return_value=candidate_a,
        ), patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            return_value=trial_a,
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once",
            return_value=rejected,
        ):
            result = self._solve(initial_basis=persistent)

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.ENRICHMENT_FAILED,
        )
        self.assertEqual(result.iterations, 0)
        self.assertEqual(result.commit_kind_history, ())
        self.assertAlmostEqual(result.final_indicator, 1.0)
        np.testing.assert_allclose(result.state.stress, self.initial.stress)
        np.testing.assert_allclose(
            result.basis.temporal_amplitude_matrix(),
            persistent.temporal_amplitude_matrix(),
        )
        self.assertFalse(
            np.allclose(
                result.basis.temporal_amplitude_matrix(),
                fixed_a.basis.temporal_amplitude_matrix(),
            )
        )

    def test_trial_b_same_baseline_commits_without_xi_b_less_than_xi_a_gate(self) -> None:
        persistent = self._basis(0.0)
        fixed_a = self._fixed(self._basis(0.5), 0.8)
        fixed_b = self._fixed(self._basis(0.8, n_modes=2), 0.4)
        enrichment = self._accepted_enrichment(fixed_b)
        candidate_a = self._candidate(self._state(0.4))
        candidate_b = self._candidate(self._state(0.8))
        trial_a = self._trial(candidate_a.state, self._state(0.35), 0.90, 1.0)
        # Deliberately worse complete LATIN indicator than Trial A.
        trial_b = self._trial(candidate_b.state, self._state(0.75), 0.95, 1.0)

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            return_value=self.local,
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            return_value=self.directions,
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            return_value=self.frozen,
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            return_value=fixed_a,
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once",
            return_value=enrichment,
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            side_effect=[candidate_a, candidate_b],
        ) as build_mock, patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            side_effect=[trial_a, trial_b],
        ) as evaluate_mock:
            result = self._solve(
                initial_basis=persistent,
                max_iterations=1,
            )

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.MAX_ITERATIONS,
        )
        self.assertEqual(result.commit_kind_history, ("B",))
        self.assertEqual(result.modes_added_history.tolist(), [1])
        self.assertAlmostEqual(result.final_indicator, 0.95)
        self.assertGreater(result.final_indicator, trial_a.indicator)
        self.assertEqual(result.basis.n_modes, 2)
        np.testing.assert_allclose(result.state.stress, trial_b.relaxed_state.stress)

        build_calls = build_mock.call_args_list
        self.assertIs(
            build_calls[0].kwargs["baseline_state"],
            build_calls[1].kwargs["baseline_state"],
        )
        self.assertIs(
            build_calls[0].kwargs["local_state"],
            build_calls[1].kwargs["local_state"],
        )
        self.assertIs(
            build_calls[0].kwargs["directions"],
            build_calls[1].kwargs["directions"],
        )
        self.assertIs(
            build_calls[0].kwargs["frozen_data"],
            build_calls[1].kwargs["frozen_data"],
        )
        eval_calls = evaluate_mock.call_args_list
        self.assertIs(
            eval_calls[0].kwargs["baseline_state"],
            eval_calls[1].kwargs["baseline_state"],
        )
        self.assertEqual(eval_calls[0].kwargs["previous_indicator"], 1.0)
        self.assertEqual(eval_calls[1].kwargs["previous_indicator"], 1.0)

    def test_empty_basis_unresolved_residual_bootstraps_first_mode(self) -> None:
        empty = self._basis(0.0, n_modes=0)
        fixed_a = self._fixed(empty, 1.0)
        fixed_b = self._fixed(self._basis(0.7), 0.2)
        enrichment = self._accepted_enrichment(fixed_b)
        candidate_a = self._candidate(self._state(0.3))
        candidate_b = self._candidate(self._state(0.7))
        # zeta = 1/3 would normally advance LATIN, but an empty unresolved
        # reduced basis must bootstrap the first pair.
        trial_a = self._trial(candidate_a.state, self._state(0.25), 0.50, 1.0)
        trial_b = self._trial(candidate_b.state, self._state(0.65), 5.0e-5, 1.0)

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            return_value=self.local,
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            return_value=self.directions,
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            return_value=self.frozen,
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            return_value=fixed_a,
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            side_effect=[candidate_a, candidate_b],
        ), patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            side_effect=[trial_a, trial_b],
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once",
            return_value=enrichment,
        ) as enrich_mock:
            result = self._solve(initial_basis=empty)

        enrich_mock.assert_called_once()
        self.assertEqual(result.commit_kind_history, ("B",))
        self.assertEqual(result.basis.n_modes, 1)
        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.CONVERGED,
        )

    def test_trial_b_hard_failure_preserves_persistent_baseline(self) -> None:
        persistent = self._basis(0.0)
        fixed_a = self._fixed(self._basis(0.5), 0.8)
        fixed_b = self._fixed(self._basis(0.8, n_modes=2), 0.4)
        enrichment = self._accepted_enrichment(fixed_b)
        candidate_a = self._candidate(self._state(0.4))
        trial_a = self._trial(candidate_a.state, self._state(0.35), 0.90, 1.0)

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            return_value=self.local,
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            return_value=self.directions,
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            return_value=self.frozen,
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            return_value=fixed_a,
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once",
            return_value=enrichment,
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            side_effect=[candidate_a, FloatingPointError("synthetic Trial-B failure")],
        ), patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            return_value=trial_a,
        ):
            result = self._solve(initial_basis=persistent)

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.TRIAL_B_FAILED,
        )
        self.assertEqual(result.iterations, 0)
        self.assertAlmostEqual(result.final_indicator, 1.0)
        np.testing.assert_allclose(result.state.stress, self.initial.stress)
        np.testing.assert_allclose(
            result.basis.temporal_amplitude_matrix(),
            persistent.temporal_amplitude_matrix(),
        )
        self.assertIsNone(result.last_trial_b)
        self.assertIsNotNone(result.last_enrichment_result)

    def test_stagnation_counter_advances_only_on_committed_iterations(self) -> None:
        persistent = self._basis(0.0)
        fixed_1 = self._fixed(self._basis(0.2), 1.0e-6)
        fixed_2 = self._fixed(self._basis(0.3), 1.0e-6)
        candidate_1 = self._candidate(self._state(0.2))
        candidate_2 = self._candidate(self._state(0.3))
        initial_xi = 8.008e-4
        trial_1 = self._trial(
            candidate_1.state,
            self._state(0.18),
            8.004e-4,
            initial_xi,
        )
        trial_2 = self._trial(
            candidate_2.state,
            self._state(0.28),
            8.001e-4,
            8.004e-4,
        )

        with patch(
            "latin.tower_latin_pgd_solver.solve_tower_local_stage",
            side_effect=[self.local, self.local],
        ), patch(
            "latin.tower_latin_pgd_solver.compute_tower_descent_search_directions",
            side_effect=[self.directions, self.directions],
        ), patch(
            "latin.tower_latin_pgd_solver.prepare_frozen_global_data",
            side_effect=[self.frozen, self.frozen],
        ), patch(
            "latin.tower_latin_pgd_solver.update_tower_pgd_time_functions",
            side_effect=[fixed_1, fixed_2],
        ), patch(
            "latin.tower_latin_pgd_solver.build_unrelaxed_candidate",
            side_effect=[candidate_1, candidate_2],
        ), patch(
            "latin.tower_latin_pgd_solver.evaluate_tower_trial",
            side_effect=[trial_1, trial_2],
        ), patch(
            "latin.tower_latin_pgd_solver.enrich_tower_pgd_basis_once"
        ) as enrich_mock:
            result = self._solve(
                initial_basis=persistent,
                initial_indicator=initial_xi,
                max_iterations=5,
                stagnation_required_iterations=2,
            )

        self.assertEqual(
            result.termination_reason,
            TowerLatinPGDTerminationReason.STAGNATED,
        )
        self.assertTrue(result.converged)
        self.assertEqual(result.iterations, 2)
        self.assertEqual(result.commit_kind_history, ("A", "A"))
        enrich_mock.assert_not_called()


if __name__ == "__main__":
    unittest.main()

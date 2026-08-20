# -*- coding: utf-8 -*-
"""Transactional outer solver for the tower LATIN-PGD formulation.

The solver is the only owner of the persistent nonlinear snapshot

    (accepted_state, accepted_basis, accepted_indicator).

Every LATIN iteration keeps that snapshot read-only while it constructs one
local state, one set of descent search directions, and one FrozenGlobalData
object.  Trial A first updates only the temporal coordinates of the accepted
spatial basis.  If enrichment is genuinely required, exactly one new PGD pair
is attempted from the provisional Trial-A basis, and Trial B is rebuilt from
the same persistent LATIN state baseline.

No failed/rejected enrichment is allowed to promote Trial A to persistent
state.  A valid Trial B, however, is committed without an additional
``xi_B < xi_A`` gate because mode acceptance has already been decided from the
full reduced mechanical residual benefit.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Optional, Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from latin.pgd_basis import PGDBasisTower
from latin.pgd_saturation import decide_pgd_saturation
from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_global_finishing import (
    build_unrelaxed_candidate,
    prepare_frozen_global_data,
)
from latin.tower_iteration_control import (
    TowerTrialEvaluation,
    evaluate_tower_trial,
)
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_pgd_enrichment import (
    TowerEnrichmentResult,
    enrich_tower_pgd_basis_once,
)
from latin.tower_pgd_time_update import (
    FixedBasisPGDResult,
    update_tower_pgd_time_functions,
)
from latin.tower_search_directions import (
    compute_tower_descent_search_directions,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]
MaterialInput = Union[
    MaterialParameters,
    Sequence[MaterialParameters],
]


class TowerLatinPGDTerminationReason(Enum):
    """Reason for terminating the transactional tower LATIN-PGD solver."""

    CONVERGED = "converged"
    STAGNATED = "stagnated"
    MAX_ITERATIONS = "max_iterations"
    ENRICHMENT_FAILED = "enrichment_failed"
    TRIAL_B_FAILED = "trial_b_failed"


def _readonly_float_history(
    values: Sequence[float],
    *,
    allow_infinity: bool = False,
) -> FloatArray:
    result = np.asarray(values, dtype=np.float64).copy()
    if result.ndim != 1:
        raise ValueError("Solver float histories must be vectors.")
    if np.any(np.isnan(result)):
        raise ValueError("Solver float histories must not contain NaN.")
    if not allow_infinity and np.any(~np.isfinite(result)):
        raise ValueError("Solver float histories must be finite.")
    result.setflags(write=False)
    return result


def _readonly_int_history(values: Sequence[int]) -> IntArray:
    result = np.asarray(values, dtype=np.int64).copy()
    if result.ndim != 1:
        raise ValueError("Solver integer histories must be vectors.")
    result.setflags(write=False)
    return result


def _copy_directions(
    directions: DescentSearchDirections,
) -> DescentSearchDirections:
    """Detach the final diagnostic search-direction snapshot."""
    return DescentSearchDirections(
        H_sigma=np.array(directions.H_sigma, dtype=np.float64, copy=True),
        H_beta=np.array(directions.H_beta, dtype=np.float64, copy=True),
        H_R_bar=np.array(directions.H_R_bar, dtype=np.float64, copy=True),
        b_damage=np.array(directions.b_damage, dtype=np.float64, copy=True),
        regularization=float(directions.regularization),
    )


@dataclass(frozen=True)
class TowerLatinPGDResult:
    """Last accepted persistent snapshot plus transactional diagnostics."""

    state: LatinStateTower
    local_state: LatinStateTower
    directions: DescentSearchDirections
    basis: PGDBasisTower
    accepted_indicator: float

    indicator_history: FloatArray
    baseline_indicator_history: FloatArray
    trial_indicator_history: FloatArray
    saturation_history: FloatArray
    trial_basis_size_history: IntArray
    trial_reduced_residual_history: FloatArray
    modes_added_history: IntArray
    trial_kind_history: Tuple[str, ...]
    commit_kind_history: Tuple[str, ...]

    converged: bool
    iterations: int
    trial_evaluations: int
    termination_reason: TowerLatinPGDTerminationReason
    failure_reason: Optional[str]

    last_trial_a: Optional[TowerTrialEvaluation]
    last_trial_b: Optional[TowerTrialEvaluation]
    last_enrichment_result: Optional[TowerEnrichmentResult]

    tolerance: float
    relaxation: float
    saturation_enrichment_tolerance: float
    saturation_stopping_tolerance: float
    reduced_tolerance: float
    mode_significance_tolerance: float
    acceptance_tolerance: float

    def __post_init__(self) -> None:
        if not isinstance(self.state, LatinStateTower):
            raise TypeError("state must be a LatinStateTower.")
        if not isinstance(self.local_state, LatinStateTower):
            raise TypeError("local_state must be a LatinStateTower.")
        if not isinstance(self.directions, DescentSearchDirections):
            raise TypeError("directions must be DescentSearchDirections.")
        if not isinstance(self.basis, PGDBasisTower):
            raise TypeError("basis must be a PGDBasisTower.")
        if not isinstance(
            self.termination_reason,
            TowerLatinPGDTerminationReason,
        ):
            raise TypeError(
                "termination_reason must be TowerLatinPGDTerminationReason."
            )

        if self.state.field_shape != self.local_state.field_shape:
            raise ValueError("state and local_state field shapes must match.")
        if self.state.field_shape != self.basis.field_shape:
            raise ValueError("state and basis field shapes must match.")
        if self.directions.field_shape != self.state.field_shape:
            raise ValueError("directions and state field shapes must match.")

        accepted_indicator = float(self.accepted_indicator)
        if not np.isfinite(accepted_indicator) or accepted_indicator < 0.0:
            raise ValueError(
                "accepted_indicator must be finite and non-negative."
            )

        indicator_history = _readonly_float_history(self.indicator_history)
        baseline_history = _readonly_float_history(
            self.baseline_indicator_history
        )
        trial_history = _readonly_float_history(self.trial_indicator_history)
        saturation_history = _readonly_float_history(self.saturation_history)
        trial_basis_history = _readonly_int_history(
            self.trial_basis_size_history
        )
        reduced_history = _readonly_float_history(
            self.trial_reduced_residual_history,
            allow_infinity=True,
        )
        modes_history = _readonly_int_history(self.modes_added_history)

        trial_kind_history = tuple(str(value) for value in self.trial_kind_history)
        commit_kind_history = tuple(
            str(value) for value in self.commit_kind_history
        )
        if any(value not in ("A", "B") for value in trial_kind_history):
            raise ValueError("trial_kind_history may contain only 'A' or 'B'.")
        if any(value not in ("A", "B") for value in commit_kind_history):
            raise ValueError("commit_kind_history may contain only 'A' or 'B'.")

        n_trials = trial_history.size
        if saturation_history.size != n_trials:
            raise ValueError("saturation_history must match trial history.")
        if trial_basis_history.size != n_trials:
            raise ValueError("trial_basis_size_history must match trial history.")
        if reduced_history.size != n_trials:
            raise ValueError(
                "trial_reduced_residual_history must match trial history."
            )
        if len(trial_kind_history) != n_trials:
            raise ValueError("trial_kind_history must match trial history.")

        iterations = int(self.iterations)
        if iterations < 0:
            raise ValueError("iterations must be non-negative.")
        if indicator_history.size != iterations:
            raise ValueError(
                "indicator_history must contain one entry per commit."
            )
        if modes_history.size != iterations:
            raise ValueError(
                "modes_added_history must contain one entry per commit."
            )
        if len(commit_kind_history) != iterations:
            raise ValueError(
                "commit_kind_history must contain one entry per commit."
            )
        if int(self.trial_evaluations) != n_trials:
            raise ValueError(
                "trial_evaluations must equal the trial history length."
            )
        if iterations > 0:
            scale = max(1.0, abs(accepted_indicator))
            if abs(float(indicator_history[-1]) - accepted_indicator) > 1.0e-12 * scale:
                raise ValueError(
                    "accepted_indicator must equal the last committed indicator."
                )

        if self.last_trial_a is not None and not isinstance(
            self.last_trial_a,
            TowerTrialEvaluation,
        ):
            raise TypeError("last_trial_a must be TowerTrialEvaluation or None.")
        if self.last_trial_b is not None and not isinstance(
            self.last_trial_b,
            TowerTrialEvaluation,
        ):
            raise TypeError("last_trial_b must be TowerTrialEvaluation or None.")
        if self.last_enrichment_result is not None and not isinstance(
            self.last_enrichment_result,
            TowerEnrichmentResult,
        ):
            raise TypeError(
                "last_enrichment_result must be TowerEnrichmentResult or None."
            )

        for name in (
            "tolerance",
            "relaxation",
            "saturation_enrichment_tolerance",
            "saturation_stopping_tolerance",
            "reduced_tolerance",
            "mode_significance_tolerance",
            "acceptance_tolerance",
        ):
            value = float(getattr(self, name))
            if not np.isfinite(value):
                raise ValueError(name + " must be finite.")
            object.__setattr__(self, name, value)

        object.__setattr__(self, "state", self.state.copy())
        object.__setattr__(self, "local_state", self.local_state.copy())
        object.__setattr__(self, "directions", _copy_directions(self.directions))
        object.__setattr__(self, "basis", self.basis.copy())
        object.__setattr__(self, "accepted_indicator", accepted_indicator)
        object.__setattr__(self, "indicator_history", indicator_history)
        object.__setattr__(self, "baseline_indicator_history", baseline_history)
        object.__setattr__(self, "trial_indicator_history", trial_history)
        object.__setattr__(self, "saturation_history", saturation_history)
        object.__setattr__(
            self,
            "trial_basis_size_history",
            trial_basis_history,
        )
        object.__setattr__(
            self,
            "trial_reduced_residual_history",
            reduced_history,
        )
        object.__setattr__(self, "modes_added_history", modes_history)
        object.__setattr__(self, "trial_kind_history", trial_kind_history)
        object.__setattr__(self, "commit_kind_history", commit_kind_history)
        object.__setattr__(self, "converged", bool(self.converged))
        object.__setattr__(self, "iterations", iterations)
        object.__setattr__(self, "trial_evaluations", int(self.trial_evaluations))
        if self.failure_reason is not None:
            object.__setattr__(self, "failure_reason", str(self.failure_reason))

    @property
    def final_indicator(self) -> float:
        """Return the indicator belonging to the returned persistent state."""
        return float(self.accepted_indicator)

    @property
    def total_modes_added(self) -> int:
        """Return the number of new pairs committed during this solve."""
        return int(np.sum(self.modes_added_history))

    @property
    def attempted_iterations(self) -> int:
        """Return the number of outer LATIN iterations that were opened."""
        return int(self.baseline_indicator_history.size)


def _positive_finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(name + " must be positive and finite.")
    return result


def _nonnegative_finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result < 0.0:
        raise ValueError(name + " must be non-negative and finite.")
    return result


def _positive_integer(value: int, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(
        value,
        (int, np.integer),
    ):
        raise TypeError(name + " must be an integer.")
    result = int(value)
    if result < 1:
        raise ValueError(name + " must be at least one.")
    return result


def _normalise_materials(
    materials: MaterialInput,
    n_material_points: int,
) -> Tuple[MaterialParameters, ...]:
    """Return exactly one material object per canonical q-point."""
    if isinstance(materials, MaterialParameters):
        result = (materials,) * n_material_points
    else:
        result = tuple(materials)
        if len(result) == 1:
            result = result * n_material_points

    if len(result) != n_material_points:
        raise ValueError(
            "materials must contain one MaterialParameters object or one "
            "object per material point."
        )
    if any(not isinstance(material, MaterialParameters) for material in result):
        raise TypeError("materials must contain only MaterialParameters objects.")
    return result


def _validate_solver_inputs(
    initial_state: LatinStateTower,
    materials: MaterialInput,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
    initial_basis: Optional[PGDBasisTower],
    initial_indicator: float,
    tolerance: float,
    stagnation_indicator_threshold: float,
    stagnation_absolute_tolerance: float,
    stagnation_required_iterations: int,
    max_iterations: int,
    relaxation: float,
    saturation_enrichment_tolerance: float,
    saturation_stopping_tolerance: float,
    reduced_tolerance: float,
    search_regularization: float,
    mode_significance_tolerance: float,
    acceptance_tolerance: float,
    fixed_point_tolerance: float,
    max_fixed_point_iterations: int,
    minimum_spatial_norm: float,
    minimum_spatial_novelty: float,
    basis_health_tolerance: float,
    field_invariance_tolerance: float,
    reorthogonalization_passes: int,
    rcond: float,
    spatial_strategy: str,
) -> Tuple[MaterialParameters, ...]:
    if not isinstance(initial_state, LatinStateTower):
        raise TypeError("initial_state must be a LatinStateTower.")
    if not isinstance(metric, MaterialPointMetric):
        raise TypeError("metric must be a MaterialPointMetric.")
    if not isinstance(equilibrium_operator, TowerEquilibriumOperator):
        raise TypeError(
            "equilibrium_operator must be a TowerEquilibriumOperator."
        )
    if metric.n_material_points != initial_state.n_material_points:
        raise ValueError("metric and initial_state material-point counts differ.")
    if equilibrium_operator.n_material_points != initial_state.n_material_points:
        raise ValueError(
            "equilibrium_operator and initial_state material-point counts differ."
        )
    if not np.array_equal(metric.weights, equilibrium_operator.metric.weights):
        raise ValueError("metric must match equilibrium_operator.metric.")

    if initial_basis is not None:
        if not isinstance(initial_basis, PGDBasisTower):
            raise TypeError("initial_basis must be a PGDBasisTower or None.")
        if initial_basis.field_shape != initial_state.field_shape:
            raise ValueError("initial_basis and initial_state field shapes differ.")

    material_tuple = _normalise_materials(
        materials,
        initial_state.n_material_points,
    )
    elastic_modulus = np.asarray(
        [material.E for material in material_tuple],
        dtype=np.float64,
    )
    modulus_scale = max(
        1.0,
        float(np.max(np.abs(equilibrium_operator.reference_modulus))),
    )
    if not np.allclose(
        elastic_modulus,
        equilibrium_operator.reference_modulus,
        rtol=0.0,
        atol=1.0e-12 * modulus_scale,
    ):
        raise ValueError(
            "materials.E must match equilibrium_operator.reference_modulus."
        )

    _positive_finite(initial_indicator, "initial_indicator")
    tolerance_value = _positive_finite(tolerance, "tolerance")
    stagnation_threshold = _positive_finite(
        stagnation_indicator_threshold,
        "stagnation_indicator_threshold",
    )
    if tolerance_value > stagnation_threshold:
        raise ValueError(
            "tolerance must not exceed stagnation_indicator_threshold."
        )
    _positive_finite(
        stagnation_absolute_tolerance,
        "stagnation_absolute_tolerance",
    )
    _positive_integer(
        stagnation_required_iterations,
        "stagnation_required_iterations",
    )
    _positive_integer(max_iterations, "max_iterations")

    relaxation_value = _positive_finite(relaxation, "relaxation")
    if relaxation_value > 1.0:
        raise ValueError("relaxation must satisfy 0 < relaxation <= 1.")

    enrichment_value = _positive_finite(
        saturation_enrichment_tolerance,
        "saturation_enrichment_tolerance",
    )
    stopping_value = _nonnegative_finite(
        saturation_stopping_tolerance,
        "saturation_stopping_tolerance",
    )
    if stopping_value >= enrichment_value:
        raise ValueError(
            "saturation_stopping_tolerance must be smaller than "
            "saturation_enrichment_tolerance."
        )

    _positive_finite(reduced_tolerance, "reduced_tolerance")
    _positive_finite(search_regularization, "search_regularization")
    _nonnegative_finite(
        mode_significance_tolerance,
        "mode_significance_tolerance",
    )
    _nonnegative_finite(acceptance_tolerance, "acceptance_tolerance")
    _positive_finite(fixed_point_tolerance, "fixed_point_tolerance")
    _positive_integer(
        max_fixed_point_iterations,
        "max_fixed_point_iterations",
    )
    _positive_finite(minimum_spatial_norm, "minimum_spatial_norm")
    _nonnegative_finite(
        minimum_spatial_novelty,
        "minimum_spatial_novelty",
    )
    _positive_finite(basis_health_tolerance, "basis_health_tolerance")
    _positive_finite(
        field_invariance_tolerance,
        "field_invariance_tolerance",
    )
    passes = _positive_integer(
        reorthogonalization_passes,
        "reorthogonalization_passes",
    )
    if passes not in (1, 2):
        raise ValueError("reorthogonalization_passes must be 1 or 2.")
    _positive_finite(rcond, "rcond")
    if spatial_strategy not in ("paper_galerkin", "residual_ls"):
        raise ValueError(
            "spatial_strategy must be 'paper_galerkin' or 'residual_ls'."
        )

    return material_tuple


def _updated_stagnation_count(
    previous_indicator: float,
    current_indicator: float,
    current_count: int,
    indicator_threshold: float,
    absolute_tolerance: float,
) -> int:
    """Update stagnation only after an atomic persistent commit."""
    absolute_change = abs(current_indicator - previous_indicator)
    if (
        current_indicator <= indicator_threshold
        and absolute_change <= absolute_tolerance
    ):
        return current_count + 1
    return 0


def _record_trial(
    trial: TowerTrialEvaluation,
    fixed_result: FixedBasisPGDResult,
    kind: str,
    trial_indicator_values: list,
    saturation_values: list,
    trial_basis_sizes: list,
    trial_reduced_residual_values: list,
    trial_kind_values: list,
) -> None:
    trial_indicator_values.append(float(trial.indicator))
    saturation_values.append(float(trial.saturation))
    trial_basis_sizes.append(int(fixed_result.basis.n_modes))
    trial_reduced_residual_values.append(float(fixed_result.relative_residual))
    trial_kind_values.append(str(kind))


def _build_result(
    *,
    state: LatinStateTower,
    local_state: LatinStateTower,
    directions: DescentSearchDirections,
    basis: PGDBasisTower,
    accepted_indicator: float,
    indicator_values: Sequence[float],
    baseline_indicator_values: Sequence[float],
    trial_indicator_values: Sequence[float],
    saturation_values: Sequence[float],
    trial_basis_sizes: Sequence[int],
    trial_reduced_residual_values: Sequence[float],
    modes_added_values: Sequence[int],
    trial_kind_values: Sequence[str],
    commit_kind_values: Sequence[str],
    converged: bool,
    termination_reason: TowerLatinPGDTerminationReason,
    failure_reason: Optional[str],
    last_trial_a: Optional[TowerTrialEvaluation],
    last_trial_b: Optional[TowerTrialEvaluation],
    last_enrichment_result: Optional[TowerEnrichmentResult],
    tolerance: float,
    relaxation: float,
    saturation_enrichment_tolerance: float,
    saturation_stopping_tolerance: float,
    reduced_tolerance: float,
    mode_significance_tolerance: float,
    acceptance_tolerance: float,
) -> TowerLatinPGDResult:
    return TowerLatinPGDResult(
        state=state,
        local_state=local_state,
        directions=directions,
        basis=basis,
        accepted_indicator=float(accepted_indicator),
        indicator_history=np.asarray(indicator_values, dtype=np.float64),
        baseline_indicator_history=np.asarray(
            baseline_indicator_values,
            dtype=np.float64,
        ),
        trial_indicator_history=np.asarray(
            trial_indicator_values,
            dtype=np.float64,
        ),
        saturation_history=np.asarray(saturation_values, dtype=np.float64),
        trial_basis_size_history=np.asarray(
            trial_basis_sizes,
            dtype=np.int64,
        ),
        trial_reduced_residual_history=np.asarray(
            trial_reduced_residual_values,
            dtype=np.float64,
        ),
        modes_added_history=np.asarray(modes_added_values, dtype=np.int64),
        trial_kind_history=tuple(trial_kind_values),
        commit_kind_history=tuple(commit_kind_values),
        converged=bool(converged),
        iterations=len(indicator_values),
        trial_evaluations=len(trial_indicator_values),
        termination_reason=termination_reason,
        failure_reason=failure_reason,
        last_trial_a=last_trial_a,
        last_trial_b=last_trial_b,
        last_enrichment_result=last_enrichment_result,
        tolerance=float(tolerance),
        relaxation=float(relaxation),
        saturation_enrichment_tolerance=float(
            saturation_enrichment_tolerance
        ),
        saturation_stopping_tolerance=float(
            saturation_stopping_tolerance
        ),
        reduced_tolerance=float(reduced_tolerance),
        mode_significance_tolerance=float(mode_significance_tolerance),
        acceptance_tolerance=float(acceptance_tolerance),
    )


def solve_tower_latin_pgd(
    initial_state: LatinStateTower,
    materials: MaterialInput,
    metric: MaterialPointMetric,
    equilibrium_operator: TowerEquilibriumOperator,
    *,
    mode_significance_tolerance: float,
    acceptance_tolerance: float,
    initial_basis: Optional[PGDBasisTower] = None,
    initial_indicator: float = 1.0,
    tolerance: float = 1.0e-4,
    stagnation_indicator_threshold: float = 1.0e-3,
    stagnation_absolute_tolerance: float = 1.0e-6,
    stagnation_required_iterations: int = 3,
    max_iterations: int = 50,
    relaxation: float = 0.8,
    saturation_enrichment_tolerance: float = 0.1,
    saturation_stopping_tolerance: float = 1.0e-4,
    reduced_tolerance: float = 1.0e-4,
    search_regularization: float = 0.15,
    fixed_point_tolerance: float = 1.0e-6,
    max_fixed_point_iterations: int = 30,
    minimum_spatial_norm: float = 1.0e-14,
    minimum_spatial_novelty: float = 1.0e-12,
    basis_health_tolerance: float = 1.0e-8,
    field_invariance_tolerance: float = 1.0e-10,
    reorthogonalization_passes: int = 2,
    rcond: float = 1.0e-12,
    spatial_strategy: str = "paper_galerkin",
) -> TowerLatinPGDResult:
    """Run tower LATIN-PGD iterations with strict Trial-A/B transactions."""
    material_tuple = _validate_solver_inputs(
        initial_state=initial_state,
        materials=materials,
        metric=metric,
        equilibrium_operator=equilibrium_operator,
        initial_basis=initial_basis,
        initial_indicator=initial_indicator,
        tolerance=tolerance,
        stagnation_indicator_threshold=stagnation_indicator_threshold,
        stagnation_absolute_tolerance=stagnation_absolute_tolerance,
        stagnation_required_iterations=stagnation_required_iterations,
        max_iterations=max_iterations,
        relaxation=relaxation,
        saturation_enrichment_tolerance=saturation_enrichment_tolerance,
        saturation_stopping_tolerance=saturation_stopping_tolerance,
        reduced_tolerance=reduced_tolerance,
        search_regularization=search_regularization,
        mode_significance_tolerance=mode_significance_tolerance,
        acceptance_tolerance=acceptance_tolerance,
        fixed_point_tolerance=fixed_point_tolerance,
        max_fixed_point_iterations=max_fixed_point_iterations,
        minimum_spatial_norm=minimum_spatial_norm,
        minimum_spatial_novelty=minimum_spatial_novelty,
        basis_health_tolerance=basis_health_tolerance,
        field_invariance_tolerance=field_invariance_tolerance,
        reorthogonalization_passes=reorthogonalization_passes,
        rcond=rcond,
        spatial_strategy=spatial_strategy,
    )

    accepted_state = initial_state.copy()
    if initial_basis is None:
        accepted_basis = PGDBasisTower(
            n_material_points=initial_state.n_material_points,
            n_time=initial_state.n_time,
        )
    else:
        accepted_basis = initial_basis.copy()
    accepted_indicator = float(initial_indicator)

    indicator_values = []
    baseline_indicator_values = []
    trial_indicator_values = []
    saturation_values = []
    trial_basis_sizes = []
    trial_reduced_residual_values = []
    modes_added_values = []
    trial_kind_values = []
    commit_kind_values = []

    stagnation_count = 0
    last_local_state = accepted_state.copy()
    last_directions: Optional[DescentSearchDirections] = None
    last_trial_a: Optional[TowerTrialEvaluation] = None
    last_trial_b: Optional[TowerTrialEvaluation] = None
    last_enrichment_result: Optional[TowerEnrichmentResult] = None

    for latin_iteration in range(1, int(max_iterations) + 1):
        # Persistent transaction baseline for this complete iteration.
        baseline_state = accepted_state
        baseline_basis = accepted_basis
        baseline_indicator = accepted_indicator
        baseline_indicator_values.append(baseline_indicator)

        local_state = solve_tower_local_stage(
            global_state=baseline_state,
            materials=material_tuple,
        )
        directions = compute_tower_descent_search_directions(
            local_state=local_state,
            materials=material_tuple,
            regularization=search_regularization,
        )
        frozen_data = prepare_frozen_global_data(
            baseline_state=baseline_state,
            local_state=local_state,
            directions=directions,
            equilibrium_operator=equilibrium_operator,
        )
        last_local_state = local_state
        last_directions = directions
        last_trial_b = None
        last_enrichment_result = None

        fixed_a = update_tower_pgd_time_functions(
            basis=baseline_basis,
            time=baseline_state.time,
            forcing=frozen_data.full_plastic_forcing,
            H_sigma=directions.H_sigma,
            metric=metric,
            equilibrium_operator=equilibrium_operator,
            reduced_tolerance=reduced_tolerance,
            rcond=rcond,
        )
        candidate_a = build_unrelaxed_candidate(
            baseline_state=baseline_state,
            local_state=local_state,
            directions=directions,
            frozen_data=frozen_data,
            fixed_basis_result=fixed_a,
            materials=material_tuple,
        )
        trial_a = evaluate_tower_trial(
            baseline_state=baseline_state,
            local_state=local_state,
            unrelaxed_state=candidate_a.state,
            directions=directions,
            metric=metric,
            reference_modulus=equilibrium_operator.reference_modulus,
            previous_indicator=baseline_indicator,
            relaxation=relaxation,
            latin_tolerance=tolerance,
        )
        last_trial_a = trial_a
        _record_trial(
            trial=trial_a,
            fixed_result=fixed_a,
            kind="A",
            trial_indicator_values=trial_indicator_values,
            saturation_values=saturation_values,
            trial_basis_sizes=trial_basis_sizes,
            trial_reduced_residual_values=trial_reduced_residual_values,
            trial_kind_values=trial_kind_values,
        )

        # Absolute LATIN convergence has first priority.
        if trial_a.converged:
            accepted_state, accepted_basis, accepted_indicator = (
                trial_a.relaxed_state.copy(),
                fixed_a.basis.copy(),
                float(trial_a.indicator),
            )
            indicator_values.append(accepted_indicator)
            modes_added_values.append(0)
            commit_kind_values.append("A")
            return _build_result(
                state=accepted_state,
                local_state=last_local_state,
                directions=directions,
                basis=accepted_basis,
                accepted_indicator=accepted_indicator,
                indicator_values=indicator_values,
                baseline_indicator_values=baseline_indicator_values,
                trial_indicator_values=trial_indicator_values,
                saturation_values=saturation_values,
                trial_basis_sizes=trial_basis_sizes,
                trial_reduced_residual_values=trial_reduced_residual_values,
                modes_added_values=modes_added_values,
                trial_kind_values=trial_kind_values,
                commit_kind_values=commit_kind_values,
                converged=True,
                termination_reason=TowerLatinPGDTerminationReason.CONVERGED,
                failure_reason=None,
                last_trial_a=last_trial_a,
                last_trial_b=None,
                last_enrichment_result=None,
                tolerance=tolerance,
                relaxation=relaxation,
                saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                saturation_stopping_tolerance=saturation_stopping_tolerance,
                reduced_tolerance=reduced_tolerance,
                mode_significance_tolerance=mode_significance_tolerance,
                acceptance_tolerance=acceptance_tolerance,
            )

        saturation_decision = decide_pgd_saturation(
            previous_indicator=baseline_indicator,
            current_indicator=trial_a.indicator,
            enrichment_tolerance=saturation_enrichment_tolerance,
            stopping_tolerance=saturation_stopping_tolerance,
        )

        # Empty basis + unresolved reduced residual is the first-mode bootstrap.
        bootstrap_first_mode = bool(
            fixed_a.basis.n_modes == 0
            and fixed_a.relative_residual > reduced_tolerance
        )

        accept_trial_a = False
        if not bootstrap_first_mode:
            if saturation_decision.should_advance_latin:
                accept_trial_a = True
            elif fixed_a.relative_residual <= reduced_tolerance:
                # Very-small/negative saturation cannot manufacture a useful
                # residual-driven mode once the reduced problem is solved.
                accept_trial_a = True

        if accept_trial_a:
            previous_indicator = baseline_indicator
            accepted_state, accepted_basis, accepted_indicator = (
                trial_a.relaxed_state.copy(),
                fixed_a.basis.copy(),
                float(trial_a.indicator),
            )
            indicator_values.append(accepted_indicator)
            modes_added_values.append(0)
            commit_kind_values.append("A")
            stagnation_count = _updated_stagnation_count(
                previous_indicator=previous_indicator,
                current_indicator=accepted_indicator,
                current_count=stagnation_count,
                indicator_threshold=stagnation_indicator_threshold,
                absolute_tolerance=stagnation_absolute_tolerance,
            )
            if stagnation_count >= stagnation_required_iterations:
                return _build_result(
                    state=accepted_state,
                    local_state=last_local_state,
                    directions=directions,
                    basis=accepted_basis,
                    accepted_indicator=accepted_indicator,
                    indicator_values=indicator_values,
                    baseline_indicator_values=baseline_indicator_values,
                    trial_indicator_values=trial_indicator_values,
                    saturation_values=saturation_values,
                    trial_basis_sizes=trial_basis_sizes,
                    trial_reduced_residual_values=trial_reduced_residual_values,
                    modes_added_values=modes_added_values,
                    trial_kind_values=trial_kind_values,
                    commit_kind_values=commit_kind_values,
                    converged=True,
                    termination_reason=TowerLatinPGDTerminationReason.STAGNATED,
                    failure_reason=None,
                    last_trial_a=last_trial_a,
                    last_trial_b=None,
                    last_enrichment_result=None,
                    tolerance=tolerance,
                    relaxation=relaxation,
                    saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                    saturation_stopping_tolerance=saturation_stopping_tolerance,
                    reduced_tolerance=reduced_tolerance,
                    mode_significance_tolerance=mode_significance_tolerance,
                    acceptance_tolerance=acceptance_tolerance,
                )
            continue

        # At this point enrichment is genuinely required and the reduced
        # residual is unresolved.  The working baseline is B_m^A, whereas the
        # whole-iteration rollback target remains persistent B_m.
        enrichment = enrich_tower_pgd_basis_once(
            fixed_basis_result=fixed_a,
            time=baseline_state.time,
            full_forcing=frozen_data.full_plastic_forcing,
            shifted_defect=fixed_a.mechanical_residual,
            H_sigma=directions.H_sigma,
            metric=metric,
            equilibrium_operator=equilibrium_operator,
            mode_significance_tolerance=mode_significance_tolerance,
            acceptance_tolerance=acceptance_tolerance,
            iteration_added=latin_iteration,
            fixed_point_tolerance=fixed_point_tolerance,
            max_fixed_point_iterations=max_fixed_point_iterations,
            minimum_spatial_norm=minimum_spatial_norm,
            minimum_spatial_novelty=minimum_spatial_novelty,
            basis_health_tolerance=basis_health_tolerance,
            field_invariance_tolerance=field_invariance_tolerance,
            reorthogonalization_passes=reorthogonalization_passes,
            reduced_tolerance=reduced_tolerance,
            rcond=rcond,
            spatial_strategy=spatial_strategy,
        )
        last_enrichment_result = enrichment

        if not enrichment.accepted:
            return _build_result(
                state=baseline_state,
                local_state=last_local_state,
                directions=directions,
                basis=baseline_basis,
                accepted_indicator=baseline_indicator,
                indicator_values=indicator_values,
                baseline_indicator_values=baseline_indicator_values,
                trial_indicator_values=trial_indicator_values,
                saturation_values=saturation_values,
                trial_basis_sizes=trial_basis_sizes,
                trial_reduced_residual_values=trial_reduced_residual_values,
                modes_added_values=modes_added_values,
                trial_kind_values=trial_kind_values,
                commit_kind_values=commit_kind_values,
                converged=False,
                termination_reason=(
                    TowerLatinPGDTerminationReason.ENRICHMENT_FAILED
                ),
                failure_reason=enrichment.failure_reason,
                last_trial_a=last_trial_a,
                last_trial_b=None,
                last_enrichment_result=last_enrichment_result,
                tolerance=tolerance,
                relaxation=relaxation,
                saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                saturation_stopping_tolerance=saturation_stopping_tolerance,
                reduced_tolerance=reduced_tolerance,
                mode_significance_tolerance=mode_significance_tolerance,
                acceptance_tolerance=acceptance_tolerance,
            )

        fixed_b = enrichment.candidate_fixed_basis_result
        if fixed_b is None or fixed_b.basis.n_modes != fixed_a.basis.n_modes + 1:
            return _build_result(
                state=baseline_state,
                local_state=last_local_state,
                directions=directions,
                basis=baseline_basis,
                accepted_indicator=baseline_indicator,
                indicator_values=indicator_values,
                baseline_indicator_values=baseline_indicator_values,
                trial_indicator_values=trial_indicator_values,
                saturation_values=saturation_values,
                trial_basis_sizes=trial_basis_sizes,
                trial_reduced_residual_values=trial_reduced_residual_values,
                modes_added_values=modes_added_values,
                trial_kind_values=trial_kind_values,
                commit_kind_values=commit_kind_values,
                converged=False,
                termination_reason=(
                    TowerLatinPGDTerminationReason.ENRICHMENT_FAILED
                ),
                failure_reason="accepted_enrichment_invalid_basis_increment",
                last_trial_a=last_trial_a,
                last_trial_b=None,
                last_enrichment_result=last_enrichment_result,
                tolerance=tolerance,
                relaxation=relaxation,
                saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                saturation_stopping_tolerance=saturation_stopping_tolerance,
                reduced_tolerance=reduced_tolerance,
                mode_significance_tolerance=mode_significance_tolerance,
                acceptance_tolerance=acceptance_tolerance,
            )

        # Trial B is rebuilt from the same persistent state/local/search/frozen
        # inputs.  No Trial-A state chaining is permitted.
        try:
            candidate_b = build_unrelaxed_candidate(
                baseline_state=baseline_state,
                local_state=local_state,
                directions=directions,
                frozen_data=frozen_data,
                fixed_basis_result=fixed_b,
                materials=material_tuple,
            )
            trial_b = evaluate_tower_trial(
                baseline_state=baseline_state,
                local_state=local_state,
                unrelaxed_state=candidate_b.state,
                directions=directions,
                metric=metric,
                reference_modulus=equilibrium_operator.reference_modulus,
                previous_indicator=baseline_indicator,
                relaxation=relaxation,
                latin_tolerance=tolerance,
            )
        except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
            return _build_result(
                state=baseline_state,
                local_state=last_local_state,
                directions=directions,
                basis=baseline_basis,
                accepted_indicator=baseline_indicator,
                indicator_values=indicator_values,
                baseline_indicator_values=baseline_indicator_values,
                trial_indicator_values=trial_indicator_values,
                saturation_values=saturation_values,
                trial_basis_sizes=trial_basis_sizes,
                trial_reduced_residual_values=trial_reduced_residual_values,
                modes_added_values=modes_added_values,
                trial_kind_values=trial_kind_values,
                commit_kind_values=commit_kind_values,
                converged=False,
                termination_reason=TowerLatinPGDTerminationReason.TRIAL_B_FAILED,
                failure_reason=(
                    "trial_b_build_or_evaluation_failed:"
                    + type(error).__name__
                    + ":"
                    + str(error)
                ),
                last_trial_a=last_trial_a,
                last_trial_b=None,
                last_enrichment_result=last_enrichment_result,
                tolerance=tolerance,
                relaxation=relaxation,
                saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                saturation_stopping_tolerance=saturation_stopping_tolerance,
                reduced_tolerance=reduced_tolerance,
                mode_significance_tolerance=mode_significance_tolerance,
                acceptance_tolerance=acceptance_tolerance,
            )

        last_trial_b = trial_b
        _record_trial(
            trial=trial_b,
            fixed_result=fixed_b,
            kind="B",
            trial_indicator_values=trial_indicator_values,
            saturation_values=saturation_values,
            trial_basis_sizes=trial_basis_sizes,
            trial_reduced_residual_values=trial_reduced_residual_values,
            trial_kind_values=trial_kind_values,
        )

        # Valid Trial B commits atomically with no xi_B < xi_A requirement.
        previous_indicator = baseline_indicator
        accepted_state, accepted_basis, accepted_indicator = (
            trial_b.relaxed_state.copy(),
            fixed_b.basis.copy(),
            float(trial_b.indicator),
        )
        indicator_values.append(accepted_indicator)
        modes_added_values.append(1)
        commit_kind_values.append("B")

        if trial_b.converged:
            return _build_result(
                state=accepted_state,
                local_state=last_local_state,
                directions=directions,
                basis=accepted_basis,
                accepted_indicator=accepted_indicator,
                indicator_values=indicator_values,
                baseline_indicator_values=baseline_indicator_values,
                trial_indicator_values=trial_indicator_values,
                saturation_values=saturation_values,
                trial_basis_sizes=trial_basis_sizes,
                trial_reduced_residual_values=trial_reduced_residual_values,
                modes_added_values=modes_added_values,
                trial_kind_values=trial_kind_values,
                commit_kind_values=commit_kind_values,
                converged=True,
                termination_reason=TowerLatinPGDTerminationReason.CONVERGED,
                failure_reason=None,
                last_trial_a=last_trial_a,
                last_trial_b=last_trial_b,
                last_enrichment_result=last_enrichment_result,
                tolerance=tolerance,
                relaxation=relaxation,
                saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                saturation_stopping_tolerance=saturation_stopping_tolerance,
                reduced_tolerance=reduced_tolerance,
                mode_significance_tolerance=mode_significance_tolerance,
                acceptance_tolerance=acceptance_tolerance,
            )

        stagnation_count = _updated_stagnation_count(
            previous_indicator=previous_indicator,
            current_indicator=accepted_indicator,
            current_count=stagnation_count,
            indicator_threshold=stagnation_indicator_threshold,
            absolute_tolerance=stagnation_absolute_tolerance,
        )
        if stagnation_count >= stagnation_required_iterations:
            return _build_result(
                state=accepted_state,
                local_state=last_local_state,
                directions=directions,
                basis=accepted_basis,
                accepted_indicator=accepted_indicator,
                indicator_values=indicator_values,
                baseline_indicator_values=baseline_indicator_values,
                trial_indicator_values=trial_indicator_values,
                saturation_values=saturation_values,
                trial_basis_sizes=trial_basis_sizes,
                trial_reduced_residual_values=trial_reduced_residual_values,
                modes_added_values=modes_added_values,
                trial_kind_values=trial_kind_values,
                commit_kind_values=commit_kind_values,
                converged=True,
                termination_reason=TowerLatinPGDTerminationReason.STAGNATED,
                failure_reason=None,
                last_trial_a=last_trial_a,
                last_trial_b=last_trial_b,
                last_enrichment_result=last_enrichment_result,
                tolerance=tolerance,
                relaxation=relaxation,
                saturation_enrichment_tolerance=saturation_enrichment_tolerance,
                saturation_stopping_tolerance=saturation_stopping_tolerance,
                reduced_tolerance=reduced_tolerance,
                mode_significance_tolerance=mode_significance_tolerance,
                acceptance_tolerance=acceptance_tolerance,
            )

    if last_directions is None:
        raise RuntimeError(
            "The tower LATIN-PGD solver ended before opening an iteration."
        )

    return _build_result(
        state=accepted_state,
        local_state=last_local_state,
        directions=last_directions,
        basis=accepted_basis,
        accepted_indicator=accepted_indicator,
        indicator_values=indicator_values,
        baseline_indicator_values=baseline_indicator_values,
        trial_indicator_values=trial_indicator_values,
        saturation_values=saturation_values,
        trial_basis_sizes=trial_basis_sizes,
        trial_reduced_residual_values=trial_reduced_residual_values,
        modes_added_values=modes_added_values,
        trial_kind_values=trial_kind_values,
        commit_kind_values=commit_kind_values,
        converged=False,
        termination_reason=TowerLatinPGDTerminationReason.MAX_ITERATIONS,
        failure_reason=None,
        last_trial_a=last_trial_a,
        last_trial_b=last_trial_b,
        last_enrichment_result=last_enrichment_result,
        tolerance=tolerance,
        relaxation=relaxation,
        saturation_enrichment_tolerance=saturation_enrichment_tolerance,
        saturation_stopping_tolerance=saturation_stopping_tolerance,
        reduced_tolerance=reduced_tolerance,
        mode_significance_tolerance=mode_significance_tolerance,
        acceptance_tolerance=acceptance_tolerance,
    )

# -*- coding: utf-8 -*-
"""
Tower LATIN-PGD relaxation and iteration-level convergence diagnostics.

This module implements the control operations that act on one *complete*
unrelaxed tower global candidate:

    breve{s}_{i+1}
        -> uniform field-wise relaxation
        -> s_trial
        -> Eq. (76) LATIN indicator using Eq. (77)
        -> Eq. (60) raw saturation value.

The tower-v1 contract is intentionally narrow:

* all 13 stored material histories use the same relaxation parameter mu;
* the relaxed energy-release-rate field is the convex blend of the two stored
  fields and is NOT reprojected through Y(sigma, D);
* Eq. (77) reads only the seven mechanical groups
      sigma, beta, R_bar, eps_p_dot, eps_e, alpha_dot, r_bar_dot;
* spatial integration uses MaterialPointMetric weights v_q;
* time integration inherits the validated one-dimensional diagnostic convention
  and uses trapezoidal quadrature;
* Eq. (60) is returned only as data.  This module does not decide enrichment,
  modify a PGD basis, or commit persistent solver state.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple, Union

import numpy as np
from numpy.typing import NDArray

from latin.pgd_saturation import saturation_indicator
from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import MaterialPointMetric
from latin.tower_state import LatinStateTower


FloatArray = NDArray[np.float64]
ReferenceModulusInput = Union[float, FloatArray]


def _finite_scalar(value: float, name: str) -> float:
    """Return a finite real scalar, rejecting Boolean aliases."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be a real scalar.")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error
    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _relaxation_parameter(value: float) -> float:
    """Return a valid LATIN relaxation parameter."""
    result = _finite_scalar(value, "relaxation")
    if result <= 0.0 or result > 1.0:
        raise ValueError(
            "relaxation must satisfy 0 < relaxation <= 1."
        )
    return result


def _nonnegative_scalar(value: float, name: str) -> float:
    """Return a finite non-negative scalar."""
    result = _finite_scalar(value, name)
    if result < 0.0:
        raise ValueError(name + " must be non-negative.")
    return result


def _validate_state_pair(
    first_state: LatinStateTower,
    second_state: LatinStateTower,
) -> None:
    """Validate a common immutable tower time-material-point grid."""
    if not isinstance(first_state, LatinStateTower):
        raise TypeError("first_state must be a LatinStateTower.")
    if not isinstance(second_state, LatinStateTower):
        raise TypeError("second_state must be a LatinStateTower.")
    if first_state.field_shape != second_state.field_shape:
        raise ValueError(
            "The two tower LATIN states must have the same field shape."
        )
    if not np.array_equal(first_state.time, second_state.time):
        raise ValueError(
            "The two tower LATIN states must use the same time grid."
        )


def _validate_metric_and_directions(
    state: LatinStateTower,
    directions: DescentSearchDirections,
    metric: MaterialPointMetric,
) -> None:
    """Validate the Eq. (77) material-point metric and search directions."""
    if not isinstance(directions, DescentSearchDirections):
        raise TypeError(
            "directions must be a DescentSearchDirections object."
        )
    if not isinstance(metric, MaterialPointMetric):
        raise TypeError("metric must be a MaterialPointMetric.")
    if directions.field_shape != state.field_shape:
        raise ValueError(
            "Search directions must use the tower state field shape."
        )
    if metric.n_material_points != state.n_material_points:
        raise ValueError(
            "metric and tower state must use the same material-point count."
        )
    if np.any(directions.b_damage != 0.0):
        raise ValueError("tower v1 requires b_damage = 0.")


def _reference_modulus_vector(
    reference_modulus: ReferenceModulusInput,
    n_material_points: int,
) -> FloatArray:
    """Return one positive reference elastic modulus per material point."""
    if np.isscalar(reference_modulus):
        value = _finite_scalar(
            reference_modulus,
            "reference_modulus",
        )
        result = np.full(
            n_material_points,
            value,
            dtype=np.float64,
        )
    else:
        result = np.array(
            reference_modulus,
            dtype=np.float64,
            copy=True,
        )
        if result.shape != (n_material_points,):
            raise ValueError(
                "reference_modulus must be a scalar or have shape "
                + str((n_material_points,))
                + "."
            )
        if np.any(~np.isfinite(result)):
            raise ValueError(
                "reference_modulus must contain only finite values."
            )

    if np.any(result <= 0.0):
        raise ValueError(
            "reference_modulus must be strictly positive."
        )
    result.setflags(write=False)
    return result


def relax_tower_global_state(
    previous_state: LatinStateTower,
    candidate_state: LatinStateTower,
    *,
    relaxation: float = 0.8,
) -> LatinStateTower:
    """
    Apply the LATIN relaxation step to every stored tower material history.

    The same coefficient is used for the nine primary fields and the four
    integrated support histories.  In particular, energy_release_rate is
    blended exactly like the other stored fields; it is not recomputed from the
    relaxed stress and damage.
    """
    _validate_state_pair(
        first_state=previous_state,
        second_state=candidate_state,
    )
    mu = _relaxation_parameter(relaxation)
    one_minus_mu = 1.0 - mu

    blended_fields = {}
    for field_name in LatinStateTower.MATERIAL_FIELD_NAMES:
        blended_fields[field_name] = (
            one_minus_mu * getattr(previous_state, field_name)
            + mu * getattr(candidate_state, field_name)
        )

    return LatinStateTower(
        time=previous_state.time,
        **blended_fields,
    )


def _mechanical_fields(
    state: LatinStateTower,
) -> Tuple[
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
]:
    """Return the seven tower fields entering the scalar-fiber Eq. (77)."""
    return (
        state.stress,
        state.beta,
        state.R_bar,
        state.plastic_strain_rate,
        state.elastic_strain,
        state.alpha_rate,
        state.r_bar_rate,
    )


def _integrate_norm_density(
    density: FloatArray,
    time: FloatArray,
    metric: MaterialPointMetric,
) -> float:
    """Integrate one finite non-negative Eq. (77) density over q and time."""
    if density.shape != (time.size, metric.n_material_points):
        raise ValueError(
            "density must have shape "
            + str((time.size, metric.n_material_points))
            + "."
        )
    if np.any(~np.isfinite(density)):
        raise ValueError("Norm density contains non-finite values.")

    space_integral = np.asarray(
        metric.integrate(density),
        dtype=np.float64,
    )
    value = float(np.trapz(space_integral, x=time))

    scale = max(
        1.0,
        float(np.max(np.abs(space_integral)))
        if space_integral.size
        else 1.0,
    )
    if value < -1.0e-12 * scale:
        raise FloatingPointError(
            "The squared tower LATIN norm became negative."
        )
    return max(value, 0.0)


def _mechanical_norm_squared_from_fields(
    fields: Tuple[
        FloatArray,
        FloatArray,
        FloatArray,
        FloatArray,
        FloatArray,
        FloatArray,
        FloatArray,
    ],
    time: FloatArray,
    directions: DescentSearchDirections,
    metric: MaterialPointMetric,
    reference_modulus: FloatArray,
) -> float:
    """Evaluate the scalar fiber/material-point form of paper Eq. (77)."""
    (
        stress,
        beta,
        R_bar,
        plastic_strain_rate,
        elastic_strain,
        alpha_rate,
        r_bar_rate,
    ) = fields

    expected_shape = directions.field_shape
    for field in fields:
        if field.shape != expected_shape:
            raise ValueError(
                "Every Eq. (77) field must have the search-direction shape."
            )
        if np.any(~np.isfinite(field)):
            raise ValueError(
                "An Eq. (77) field contains non-finite values."
            )

    density = (
        stress**2 * directions.H_sigma
        + beta**2 * directions.H_beta
        + R_bar**2 * directions.H_R_bar
        + plastic_strain_rate**2 / directions.H_sigma
        + elastic_strain**2 * reference_modulus[np.newaxis, :]
        + alpha_rate**2 / directions.H_beta
        + r_bar_rate**2 / directions.H_R_bar
    )

    return _integrate_norm_density(
        density=density,
        time=time,
        metric=metric,
    )


def tower_plastic_state_norm_squared(
    state: LatinStateTower,
    directions: DescentSearchDirections,
    metric: MaterialPointMetric,
    reference_modulus: ReferenceModulusInput,
) -> float:
    """Return the tower mechanical norm squared from paper Eq. (77)."""
    if not isinstance(state, LatinStateTower):
        raise TypeError("state must be a LatinStateTower.")
    _validate_metric_and_directions(
        state=state,
        directions=directions,
        metric=metric,
    )
    C0 = _reference_modulus_vector(
        reference_modulus=reference_modulus,
        n_material_points=state.n_material_points,
    )
    return _mechanical_norm_squared_from_fields(
        fields=_mechanical_fields(state),
        time=state.time,
        directions=directions,
        metric=metric,
        reference_modulus=C0,
    )


def _relative_indicator_components(
    local_state: LatinStateTower,
    global_state: LatinStateTower,
    directions: DescentSearchDirections,
    metric: MaterialPointMetric,
    reference_modulus: ReferenceModulusInput,
) -> Tuple[float, float, float, float]:
    """Return xi together with difference, local, and global norms."""
    _validate_state_pair(
        first_state=local_state,
        second_state=global_state,
    )
    _validate_metric_and_directions(
        state=local_state,
        directions=directions,
        metric=metric,
    )
    C0 = _reference_modulus_vector(
        reference_modulus=reference_modulus,
        n_material_points=local_state.n_material_points,
    )

    local_fields = _mechanical_fields(local_state)
    global_fields = _mechanical_fields(global_state)
    difference_fields = tuple(
        local_field - global_field
        for local_field, global_field in zip(
            local_fields,
            global_fields,
        )
    )

    difference_norm = float(
        np.sqrt(
            _mechanical_norm_squared_from_fields(
                fields=difference_fields,
                time=local_state.time,
                directions=directions,
                metric=metric,
                reference_modulus=C0,
            )
        )
    )
    local_norm = float(
        np.sqrt(
            _mechanical_norm_squared_from_fields(
                fields=local_fields,
                time=local_state.time,
                directions=directions,
                metric=metric,
                reference_modulus=C0,
            )
        )
    )
    global_norm = float(
        np.sqrt(
            _mechanical_norm_squared_from_fields(
                fields=global_fields,
                time=global_state.time,
                directions=directions,
                metric=metric,
                reference_modulus=C0,
            )
        )
    )

    denominator = local_norm + global_norm
    numerical_zero = float(np.finfo(np.float64).eps)

    if denominator <= numerical_zero:
        if difference_norm <= numerical_zero:
            indicator = 0.0
        else:
            raise FloatingPointError(
                "Degenerate Eq. (76) denominator with a nonzero difference."
            )
    else:
        indicator = float(difference_norm / denominator)

    if not np.isfinite(indicator):
        raise FloatingPointError(
            "The tower relative LATIN indicator is non-finite."
        )
    if indicator < -1.0e-14 or indicator > 1.0 + 1.0e-12:
        raise FloatingPointError(
            "The tower relative LATIN indicator left its metric bound."
        )

    indicator = min(max(indicator, 0.0), 1.0)
    return (
        indicator,
        difference_norm,
        local_norm,
        global_norm,
    )


def tower_relative_latin_indicator(
    local_state: LatinStateTower,
    global_state: LatinStateTower,
    directions: DescentSearchDirections,
    metric: MaterialPointMetric,
    reference_modulus: ReferenceModulusInput,
) -> float:
    """Evaluate the tower relative LATIN indicator from Eqs. (76)-(77)."""
    return _relative_indicator_components(
        local_state=local_state,
        global_state=global_state,
        directions=directions,
        metric=metric,
        reference_modulus=reference_modulus,
    )[0]


@dataclass(frozen=True)
class TowerTrialEvaluation:
    """
    Complete iteration-control data for one tower trial.

    `saturation` is only the raw Eq. (60) value relative to the persistent
    previous_indicator supplied to `evaluate_tower_trial`.  Enrichment,
    rollback, and commit semantics belong to the future solver orchestrator.
    """

    unrelaxed_state: LatinStateTower
    relaxed_state: LatinStateTower
    indicator: float
    saturation: float
    converged: bool
    difference_norm: float
    local_norm: float
    global_norm: float
    relaxation: float
    latin_tolerance: float
    previous_indicator: float

    def __post_init__(self) -> None:
        if not isinstance(self.unrelaxed_state, LatinStateTower):
            raise TypeError(
                "unrelaxed_state must be a LatinStateTower."
            )
        if not isinstance(self.relaxed_state, LatinStateTower):
            raise TypeError(
                "relaxed_state must be a LatinStateTower."
            )
        _validate_state_pair(
            first_state=self.unrelaxed_state,
            second_state=self.relaxed_state,
        )

        for name in (
            "indicator",
            "difference_norm",
            "local_norm",
            "global_norm",
            "latin_tolerance",
            "previous_indicator",
        ):
            value = _nonnegative_scalar(getattr(self, name), name)
            object.__setattr__(self, name, value)

        saturation = _finite_scalar(
            self.saturation,
            "saturation",
        )
        if saturation < -1.0 - 1.0e-12 or saturation > 1.0 + 1.0e-12:
            raise ValueError(
                "saturation must lie within the Eq. (60) range [-1, 1]."
            )
        object.__setattr__(
            self,
            "saturation",
            min(max(saturation, -1.0), 1.0),
        )

        relaxation = _relaxation_parameter(self.relaxation)
        object.__setattr__(self, "relaxation", relaxation)
        object.__setattr__(self, "converged", bool(self.converged))

        if self.indicator > 1.0 + 1.0e-12:
            raise ValueError(
                "indicator must not exceed one for the Eq. (77) norm."
            )
        expected_converged = self.indicator <= self.latin_tolerance
        if self.converged != expected_converged:
            raise ValueError(
                "converged must equal indicator <= latin_tolerance."
            )

    @property
    def finite(self) -> bool:
        """Return True because all state and scalar invariants were validated."""
        scalars = (
            self.indicator,
            self.saturation,
            self.difference_norm,
            self.local_norm,
            self.global_norm,
            self.relaxation,
            self.latin_tolerance,
            self.previous_indicator,
        )
        return bool(all(np.isfinite(value) for value in scalars))


def evaluate_tower_trial(
    baseline_state: LatinStateTower,
    local_state: LatinStateTower,
    unrelaxed_state: LatinStateTower,
    directions: DescentSearchDirections,
    metric: MaterialPointMetric,
    reference_modulus: ReferenceModulusInput,
    previous_indicator: float,
    *,
    relaxation: float = 0.8,
    latin_tolerance: float = 1.0e-4,
) -> TowerTrialEvaluation:
    """
    Relax and evaluate one complete tower global candidate.

    The supplied `previous_indicator` is the persistent xi_i baseline.  The
    function never changes it and returns the raw Eq. (60) saturation value
    computed against the current trial indicator.
    """
    _validate_state_pair(
        first_state=baseline_state,
        second_state=local_state,
    )
    _validate_state_pair(
        first_state=baseline_state,
        second_state=unrelaxed_state,
    )
    _validate_metric_and_directions(
        state=baseline_state,
        directions=directions,
        metric=metric,
    )

    previous_xi = _nonnegative_scalar(
        previous_indicator,
        "previous_indicator",
    )
    tolerance = _nonnegative_scalar(
        latin_tolerance,
        "latin_tolerance",
    )
    mu = _relaxation_parameter(relaxation)

    relaxed_state = relax_tower_global_state(
        previous_state=baseline_state,
        candidate_state=unrelaxed_state,
        relaxation=mu,
    )
    (
        indicator,
        difference_norm,
        local_norm,
        global_norm,
    ) = _relative_indicator_components(
        local_state=local_state,
        global_state=relaxed_state,
        directions=directions,
        metric=metric,
        reference_modulus=reference_modulus,
    )

    zeta = saturation_indicator(
        previous_indicator=previous_xi,
        current_indicator=indicator,
    )

    return TowerTrialEvaluation(
        unrelaxed_state=unrelaxed_state,
        relaxed_state=relaxed_state,
        indicator=indicator,
        saturation=zeta,
        converged=indicator <= tolerance,
        difference_norm=difference_norm,
        local_norm=local_norm,
        global_norm=global_norm,
        relaxation=mu,
        latin_tolerance=tolerance,
        previous_indicator=previous_xi,
    )

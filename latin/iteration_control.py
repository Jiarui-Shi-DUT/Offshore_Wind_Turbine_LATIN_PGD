# -*- coding: utf-8 -*-
"""
Relaxation and convergence control for the one-dimensional LATIN solver.

The reference formulation applies the relaxed global update

    s_(i+1) = mu * s_breve_(i+1) + (1 - mu) * s_i,

with mu = 0.8.

The relative LATIN indicator is

    xi = ||s_hat_(i+1/2)^p - s_(i+1)^p||
         / (||s_hat_(i+1/2)^p|| + ||s_(i+1)^p||),

where the one-dimensional form of the mechanical norm is

    ||s^p||^2 = integral [
        sigma^2 H_sigma
        + beta^2 H_beta
        + R_bar^2 H_R_bar
        + eps_p_dot^2 / H_sigma
        + E eps_e^2
        + alpha_dot^2 / H_beta
        + r_bar_dot^2 / H_R_bar
    ] dOmega dt.

Damage variables do not appear explicitly in this norm, but they affect
stress, elastic strain and the search-direction operators.
"""

from __future__ import annotations

from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from latin.search_directions import DescentSearchDirections
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


def _validate_common_discretisation(
    first_state: LatinState,
    second_state: LatinState,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
) -> None:
    """Validate dimensions shared by relaxation and indicator operations."""
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if first_state.field_shape != second_state.field_shape:
        raise ValueError(
            "The two LATIN states must have the same field shape."
        )
    if first_state.field_shape != directions.field_shape:
        raise ValueError(
            "Search directions must have the LATIN state field shape."
        )
    if first_state.n_elements != mesh.n_elements:
        raise ValueError(
            "The LATIN states and finite-element mesh are inconsistent."
        )
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if not np.array_equal(first_state.time, second_state.time):
        raise ValueError(
            "The two LATIN states must use the same time grid."
        )


def relax_global_state(
    previous_state: LatinState,
    candidate_state: LatinState,
    *,
    relaxation: float = 0.8,
) -> LatinState:
    """
    Apply the LATIN relaxation step field by field.

    Parameters
    ----------
    previous_state:
        Global state s_i from the previous LATIN iteration.
    candidate_state:
        Unrelaxed global-stage state s_breve_(i+1).
    relaxation:
        Relaxation parameter mu. The paper uses mu = 0.8.
    """
    if previous_state.field_shape != candidate_state.field_shape:
        raise ValueError(
            "previous_state and candidate_state must have the same shape."
        )
    if not np.array_equal(previous_state.time, candidate_state.time):
        raise ValueError(
            "previous_state and candidate_state must use the same time grid."
        )
    if not 0.0 < relaxation <= 1.0:
        raise ValueError("relaxation must satisfy 0 < relaxation <= 1.")

    mu = float(relaxation)
    one_minus_mu = 1.0 - mu

    def blend(
        previous_field: FloatArray,
        candidate_field: FloatArray,
    ) -> FloatArray:
        return (
            one_minus_mu * previous_field
            + mu * candidate_field
        )

    return LatinState(
        time=previous_state.time.copy(),
        plastic_strain_rate=blend(
            previous_state.plastic_strain_rate,
            candidate_state.plastic_strain_rate,
        ),
        elastic_strain=blend(
            previous_state.elastic_strain,
            candidate_state.elastic_strain,
        ),
        alpha_rate=blend(
            previous_state.alpha_rate,
            candidate_state.alpha_rate,
        ),
        r_bar_rate=blend(
            previous_state.r_bar_rate,
            candidate_state.r_bar_rate,
        ),
        damage_rate=blend(
            previous_state.damage_rate,
            candidate_state.damage_rate,
        ),
        stress=blend(
            previous_state.stress,
            candidate_state.stress,
        ),
        beta=blend(
            previous_state.beta,
            candidate_state.beta,
        ),
        R_bar=blend(
            previous_state.R_bar,
            candidate_state.R_bar,
        ),
        energy_release_rate=blend(
            previous_state.energy_release_rate,
            candidate_state.energy_release_rate,
        ),
        plastic_strain=blend(
            previous_state.plastic_strain,
            candidate_state.plastic_strain,
        ),
        alpha=blend(
            previous_state.alpha,
            candidate_state.alpha,
        ),
        r_bar=blend(
            previous_state.r_bar,
            candidate_state.r_bar,
        ),
        damage=blend(
            previous_state.damage,
            candidate_state.damage,
        ),
    )


def _plastic_fields(
    state: LatinState,
) -> Tuple[
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
]:
    """Return the seven fields entering the paper's mechanical norm."""
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
    mesh: BarMesh1D,
    area: float,
) -> float:
    """Integrate a time-by-element norm density over space and time."""
    if density.shape != (time.size, mesh.n_elements):
        raise ValueError(
            "density must have shape (n_time, n_elements)."
        )
    if not np.all(np.isfinite(density)):
        raise ValueError("Norm density contains non-finite values.")

    element_volumes = area * mesh.element_lengths
    space_integral = density @ element_volumes
    value = float(np.trapz(space_integral, x=time))

    if value < -1.0e-10:
        raise FloatingPointError(
            "The squared LATIN norm became negative."
        )

    return max(value, 0.0)


def _plastic_norm_squared_from_fields(
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
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
) -> float:
    """Evaluate the discrete one-dimensional form of Eq. (77)."""
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
                "Every norm field must have the search-direction shape."
            )
        if not np.all(np.isfinite(field)):
            raise ValueError("A norm field contains non-finite values.")

    elastic_moduli = np.array(
        [material.E for material in materials],
        dtype=np.float64,
    )
    if np.any(elastic_moduli <= 0.0):
        raise ValueError("All Young's moduli must be positive.")

    density = (
        stress**2 * directions.H_sigma
        + beta**2 * directions.H_beta
        + R_bar**2 * directions.H_R_bar
        + plastic_strain_rate**2 / directions.H_sigma
        + elastic_strain**2
        * elastic_moduli[np.newaxis, :]
        + alpha_rate**2 / directions.H_beta
        + r_bar_rate**2 / directions.H_R_bar
    )

    return _integrate_norm_density(
        density=density,
        time=time,
        mesh=mesh,
        area=area,
    )


def plastic_state_norm_squared(
    state: LatinState,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
) -> float:
    """Return the squared mechanical norm ||s^p||^2 from Eq. (77)."""
    if state.field_shape != directions.field_shape:
        raise ValueError(
            "state and directions must have the same field shape."
        )
    if state.n_elements != mesh.n_elements:
        raise ValueError(
            "state and mesh must contain the same number of elements."
        )
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if area <= 0.0:
        raise ValueError("area must be positive.")

    return _plastic_norm_squared_from_fields(
        fields=_plastic_fields(state),
        time=state.time,
        directions=directions,
        mesh=mesh,
        area=area,
        materials=materials,
    )


def relative_latin_indicator(
    local_state: LatinState,
    global_state: LatinState,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
) -> float:
    """
    Evaluate the relative LATIN indicator xi from Eqs. (76)-(77).
    """
    _validate_common_discretisation(
        first_state=local_state,
        second_state=global_state,
        directions=directions,
        mesh=mesh,
        area=area,
        materials=materials,
    )

    local_fields = _plastic_fields(local_state)
    global_fields = _plastic_fields(global_state)
    difference_fields = (
        local_fields[0] - global_fields[0],
        local_fields[1] - global_fields[1],
        local_fields[2] - global_fields[2],
        local_fields[3] - global_fields[3],
        local_fields[4] - global_fields[4],
        local_fields[5] - global_fields[5],
        local_fields[6] - global_fields[6],
    )

    difference_norm = np.sqrt(
        _plastic_norm_squared_from_fields(
            fields=difference_fields,
            time=local_state.time,
            directions=directions,
            mesh=mesh,
            area=area,
            materials=materials,
        )
    )
    local_norm = np.sqrt(
        _plastic_norm_squared_from_fields(
            fields=local_fields,
            time=local_state.time,
            directions=directions,
            mesh=mesh,
            area=area,
            materials=materials,
        )
    )
    global_norm = np.sqrt(
        _plastic_norm_squared_from_fields(
            fields=global_fields,
            time=global_state.time,
            directions=directions,
            mesh=mesh,
            area=area,
            materials=materials,
        )
    )

    denominator = local_norm + global_norm

    if denominator <= np.finfo(float).eps:
        return 0.0 if difference_norm <= np.finfo(float).eps else np.inf

    indicator = float(difference_norm / denominator)

    if not np.isfinite(indicator):
        raise FloatingPointError(
            "The relative LATIN indicator is non-finite."
        )

    return indicator

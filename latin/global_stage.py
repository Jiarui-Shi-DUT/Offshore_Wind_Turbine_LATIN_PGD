# -*- coding: utf-8 -*-
"""
Full-order linear global stage for the one-dimensional LATIN solver.

This module implements the global correction before introducing PGD.  It is
used as a reference for the reduced global stage.  The stress and total-strain
corrections are split into a plastic part and a damage-residual part, following
Eqs. (32)-(57) of the reference formulation.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from latin.equilibrium_operator import (
    EquilibriumProjection,
    apply_equilibrium_operator,
)
from latin.search_directions import DescentSearchDirections
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class GlobalStageResult:
    """State and correction fields produced by one global stage."""

    state: LatinState
    plastic_strain_correction: FloatArray
    plastic_strain_rate_correction: FloatArray
    plastic_projection: EquilibriumProjection
    residual_strain: FloatArray
    damage_projection: EquilibriumProjection
    displacement_correction: FloatArray


def _validate_global_stage_inputs(
    global_state: LatinState,
    local_state: LatinState,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
) -> None:
    """Check dimensions and common time-space discretisation."""
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if global_state.field_shape != local_state.field_shape:
        raise ValueError(
            "global_state and local_state must have the same field shape."
        )
    if global_state.field_shape != directions.field_shape:
        raise ValueError(
            "Search directions must have the LATIN state field shape."
        )
    if global_state.n_elements != mesh.n_elements:
        raise ValueError(
            "The LATIN state and finite-element mesh are inconsistent."
        )
    if not np.array_equal(global_state.time, local_state.time):
        raise ValueError(
            "global_state and local_state must use the same time grid."
        )


def _compute_damage_residual_strain(
    global_state: LatinState,
    local_state: LatinState,
    materials: Sequence[MaterialParameters],
) -> FloatArray:
    """Compute Delta epsilon_R = C^(-1) Delta R from Eqs. (36)-(38)."""
    elastic_moduli = np.array(
        [material.E for material in materials],
        dtype=np.float64,
    )

    residual_stress = (
        global_state.stress
        - local_state.stress
        - elastic_moduli[np.newaxis, :]
        * (
            global_state.elastic_strain
            - local_state.elastic_strain
        )
    )
    return residual_stress / elastic_moduli[np.newaxis, :]


def _solve_plastic_correction(
    global_state: LatinState,
    local_state: LatinState,
    directions: DescentSearchDirections,
    damage_stress_correction: FloatArray,
    mesh: BarMesh1D,
    materials: Sequence[MaterialParameters],
) -> Tuple[FloatArray, FloatArray]:
    """
    Solve the full-order plastic correction with backward Euler in time.

    The correction equation is

        Delta eps_p_dot - H_sigma Delta sigma_plastic = forcing,

    while Delta sigma_plastic = C(E-I) Delta eps_p.  For a one-dimensional
    bar, the correction stress is constant in space, so every time step reduces
    to one scalar equilibrium equation.
    """
    shape = global_state.field_shape
    correction = np.zeros(shape, dtype=np.float64)
    rate_correction = np.zeros(shape, dtype=np.float64)

    forcing = (
        local_state.plastic_strain_rate
        - global_state.plastic_strain_rate
        - directions.H_sigma
        * (local_state.stress - global_state.stress)
        + directions.H_sigma * damage_stress_correction
    )

    element_lengths = mesh.element_lengths
    elastic_moduli = np.array(
        [material.E for material in materials],
        dtype=np.float64,
    )
    compliance_length = float(
        np.sum(element_lengths / elastic_moduli)
    )

    rate_correction[0, :] = forcing[0, :]

    for step in range(1, global_state.n_time):
        time_step = float(
            global_state.time[step] - global_state.time[step - 1]
        )
        known_part = (
            correction[step - 1, :]
            + time_step * forcing[step, :]
        )

        denominator = (
            compliance_length
            + time_step
            * float(
                np.dot(
                    element_lengths,
                    directions.H_sigma[step, :],
                )
            )
        )
        if denominator <= 0.0 or not np.isfinite(denominator):
            raise FloatingPointError(
                "Invalid denominator in the plastic global correction."
            )

        correction_stress = -float(
            np.dot(element_lengths, known_part)
        ) / denominator

        correction[step, :] = (
            known_part
            + time_step
            * directions.H_sigma[step, :]
            * correction_stress
        )
        rate_correction[step, :] = (
            correction[step, :] - correction[step - 1, :]
        ) / time_step

    if not np.all(np.isfinite(correction)):
        raise FloatingPointError(
            "Non-finite plastic-strain correction in the global stage."
        )
    if not np.all(np.isfinite(rate_correction)):
        raise FloatingPointError(
            "Non-finite plastic-strain-rate correction in the global stage."
        )

    return correction, rate_correction


def _update_hardening_variables(
    new_state: LatinState,
    local_state: LatinState,
    directions: DescentSearchDirections,
    materials: Sequence[MaterialParameters],
) -> None:
    """Solve the local linear hardening ODEs from the descent direction."""
    for element, material in enumerate(materials):
        new_state.alpha[0, element] = local_state.alpha[0, element]
        new_state.beta[0, element] = (
            material.C * new_state.alpha[0, element]
        )
        new_state.alpha_rate[0, element] = (
            local_state.alpha_rate[0, element]
            - directions.H_beta[0, element]
            * (
                new_state.beta[0, element]
                - local_state.beta[0, element]
            )
        )

        new_state.r_bar[0, element] = local_state.r_bar[0, element]
        new_state.R_bar[0, element] = (
            material.R_inf * new_state.r_bar[0, element]
        )
        new_state.r_bar_rate[0, element] = (
            local_state.r_bar_rate[0, element]
            - directions.H_R_bar[0, element]
            * (
                new_state.R_bar[0, element]
                - local_state.R_bar[0, element]
            )
        )

        for step in range(1, new_state.n_time):
            time_step = float(
                new_state.time[step] - new_state.time[step - 1]
            )

            H_beta = directions.H_beta[step, element]
            alpha_rhs = (
                local_state.alpha_rate[step, element]
                + H_beta * local_state.beta[step, element]
            )
            new_state.alpha[step, element] = (
                new_state.alpha[step - 1, element]
                + time_step * alpha_rhs
            ) / (1.0 + time_step * H_beta * material.C)
            new_state.alpha_rate[step, element] = (
                new_state.alpha[step, element]
                - new_state.alpha[step - 1, element]
            ) / time_step
            new_state.beta[step, element] = (
                material.C * new_state.alpha[step, element]
            )

            H_R_bar = directions.H_R_bar[step, element]
            r_bar_rhs = (
                local_state.r_bar_rate[step, element]
                + H_R_bar * local_state.R_bar[step, element]
            )
            new_state.r_bar[step, element] = (
                new_state.r_bar[step - 1, element]
                + time_step * r_bar_rhs
            ) / (
                1.0
                + time_step * H_R_bar * material.R_inf
            )
            new_state.r_bar_rate[step, element] = (
                new_state.r_bar[step, element]
                - new_state.r_bar[step - 1, element]
            ) / time_step
            new_state.R_bar[step, element] = (
                material.R_inf * new_state.r_bar[step, element]
            )


def _update_damage_and_energy(
    new_state: LatinState,
    local_state: LatinState,
    materials: Sequence[MaterialParameters],
) -> None:
    """Apply b_minus = 0 and post-process the damage state law Y."""
    new_state.damage_rate[:, :] = local_state.damage_rate

    for element, material in enumerate(materials):
        new_state.damage[0, element] = local_state.damage[0, element]

        for step in range(1, new_state.n_time):
            time_step = float(
                new_state.time[step] - new_state.time[step - 1]
            )
            new_state.damage[step, element] = np.clip(
                new_state.damage[step - 1, element]
                + time_step * new_state.damage_rate[step, element],
                0.0,
                material.damage_upper_bound,
            )

        stress = new_state.stress[:, element]
        damage = new_state.damage[:, element]
        tensile = stress >= 0.0

        new_state.energy_release_rate[:, element] = np.where(
            tensile,
            stress**2
            / (2.0 * material.E * (1.0 - damage) ** 2),
            material.h
            * stress**2
            / (
                2.0
                * material.E
                * (1.0 - material.h * damage) ** 2
            ),
        )


def solve_full_order_global_stage(
    global_state: LatinState,
    local_state: LatinState,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
) -> GlobalStageResult:
    """
    Perform one full-order LATIN global stage over the whole time domain.

    This routine is deliberately unreduced.  It provides the reference global
    correction that the subsequent PGD implementation must reproduce.
    """
    _validate_global_stage_inputs(
        global_state=global_state,
        local_state=local_state,
        directions=directions,
        mesh=mesh,
        area=area,
        materials=materials,
    )

    residual_strain = _compute_damage_residual_strain(
        global_state=global_state,
        local_state=local_state,
        materials=materials,
    )
    damage_projection = apply_equilibrium_operator(
        mesh=mesh,
        area=area,
        materials=materials,
        source_strain=residual_strain,
    )

    (
        plastic_strain_correction,
        plastic_strain_rate_correction,
    ) = _solve_plastic_correction(
        global_state=global_state,
        local_state=local_state,
        directions=directions,
        damage_stress_correction=damage_projection.stress,
        mesh=mesh,
        materials=materials,
    )
    plastic_projection = apply_equilibrium_operator(
        mesh=mesh,
        area=area,
        materials=materials,
        source_strain=plastic_strain_correction,
    )

    new_state = global_state.copy()
    new_state.plastic_strain[:, :] = (
        global_state.plastic_strain
        + plastic_strain_correction
    )
    new_state.plastic_strain_rate[:, :] = (
        global_state.plastic_strain_rate
        + plastic_strain_rate_correction
    )
    new_state.stress[:, :] = (
        global_state.stress
        + plastic_projection.stress
        + damage_projection.stress
    )
    new_state.elastic_strain[:, :] = (
        global_state.elastic_strain
        + plastic_projection.compatible_strain
        - plastic_strain_correction
        + damage_projection.compatible_strain
    )

    _update_hardening_variables(
        new_state=new_state,
        local_state=local_state,
        directions=directions,
        materials=materials,
    )
    _update_damage_and_energy(
        new_state=new_state,
        local_state=local_state,
        materials=materials,
    )

    displacement_correction = (
        plastic_projection.displacement
        + damage_projection.displacement
    )

    return GlobalStageResult(
        state=new_state,
        plastic_strain_correction=plastic_strain_correction,
        plastic_strain_rate_correction=(
            plastic_strain_rate_correction
        ),
        plastic_projection=plastic_projection,
        residual_strain=residual_strain,
        damage_projection=damage_projection,
        displacement_correction=displacement_correction,
    )

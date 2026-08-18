# -*- coding: utf-8 -*-
"""
Tower-compatible LATIN descent search directions on material points q.

This module reuses the current scalar-fiber tangent approximation without
changing its constitutive meaning.  The one-dimensional element index is
replaced by the canonical tower material-point index q.

For every (time, q) entry,

    H_minus = diag(H_sigma, H_beta, H_R_bar),

with the same regularisation

    H_minus <- H_minus + zeta * M_material^(-1),

where M_material = diag(E, C, R_inf), and tower v1 keeps

    b_damage = 0.

No FEM topology, PGD basis, state relaxation, or persistent solver transaction
is owned here.
"""

from __future__ import annotations

from typing import Sequence, Tuple, Union

import numpy as np

from latin.local_stage import (
    isotropic_force_from_transformed_force,
)
from latin.search_directions import DescentSearchDirections
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


MaterialInput = Union[
    MaterialParameters,
    Sequence[MaterialParameters],
]


def _material_sequence(
    materials: MaterialInput,
    n_material_points: int,
) -> Tuple[MaterialParameters, ...]:
    """Broadcast one material or validate one material object per q-point."""
    if isinstance(materials, MaterialParameters):
        material_tuple = (materials,) * n_material_points
    else:
        material_tuple = tuple(materials)
        if len(material_tuple) == 1:
            material_tuple = material_tuple * n_material_points

    if len(material_tuple) != n_material_points:
        raise ValueError(
            "materials must contain one MaterialParameters object or one "
            "object per material point."
        )
    if any(
        not isinstance(material, MaterialParameters)
        for material in material_tuple
    ):
        raise TypeError(
            "materials must contain only MaterialParameters objects."
        )

    return material_tuple


def _positive_part(values: np.ndarray) -> np.ndarray:
    """Return the componentwise positive part."""
    return np.maximum(values, 0.0)


def compute_tower_descent_search_directions(
    local_state: LatinStateTower,
    materials: MaterialInput,
    *,
    regularization: float = 0.15,
) -> DescentSearchDirections:
    """
    Compute regularised diagonal descent operators on the tower (time, q) grid.
    """
    if not isinstance(local_state, LatinStateTower):
        raise TypeError(
            "local_state must be a LatinStateTower."
        )

    regularization_value = float(regularization)
    if (
        not np.isfinite(regularization_value)
        or regularization_value <= 0.0
    ):
        raise ValueError(
            "regularization must be positive and finite."
        )

    material_tuple = _material_sequence(
        materials=materials,
        n_material_points=local_state.n_material_points,
    )

    shape = local_state.field_shape
    H_sigma = np.zeros(shape, dtype=np.float64)
    H_beta = np.zeros(shape, dtype=np.float64)
    H_R_bar = np.zeros(shape, dtype=np.float64)
    b_damage = np.zeros(shape, dtype=np.float64)

    for q, material in enumerate(material_tuple):
        if material.E <= 0.0:
            raise ValueError(
                "Every material-point Young's modulus must be positive."
            )
        if material.C <= 0.0:
            raise ValueError(
                "Every material-point kinematic modulus C must be positive."
            )
        if material.R_inf <= 0.0:
            raise ValueError(
                "Every material-point R_inf must be positive."
            )
        if material.gamma <= 0.0:
            raise ValueError(
                "Every material-point gamma must be positive."
            )

        stress = local_state.stress[:, q]
        beta = local_state.beta[:, q]
        transformed_force = local_state.R_bar[:, q]
        damage = np.clip(
            local_state.damage[:, q],
            0.0,
            material.damage_upper_bound,
        )

        relative_effective_stress = (
            stress / (1.0 - damage) - beta
        )
        flow_direction = np.sign(relative_effective_stress)

        isotropic_force = np.empty_like(transformed_force)
        for step, value in enumerate(transformed_force):
            isotropic_force[step] = (
                isotropic_force_from_transformed_force(
                    transformed_force=float(value),
                    material=material,
                )
            )

        yield_function = (
            np.abs(relative_effective_stress)
            + material.a * beta**2 / (2.0 * material.C)
            - isotropic_force
            - material.sigma_y
        )
        positive_yield = _positive_part(yield_function)

        common_tangent_factor = (
            material.k_viscoplastic
            * material.n
            * positive_yield ** (material.n - 1.0)
        )
        plastic_multiplier = (
            material.k_viscoplastic
            * positive_yield**material.n
        )

        H_sigma[:, q] = (
            common_tangent_factor / (1.0 - damage) ** 2
            + regularization_value / material.E
        )

        beta_gradient = (
            -flow_direction
            + material.a * beta / material.C
        )
        H_beta[:, q] = (
            common_tangent_factor * beta_gradient**2
            + plastic_multiplier * material.a / material.C
            + regularization_value / material.C
        )

        transformed_gradient = (
            1.0
            - transformed_force
            * np.sqrt(material.gamma)
            / (2.0 * material.R_inf)
        )
        H_R_bar[:, q] = (
            common_tangent_factor
            * material.gamma
            * transformed_gradient**2
            + plastic_multiplier
            * material.gamma
            / (2.0 * material.R_inf)
            + regularization_value / material.R_inf
        )

    return DescentSearchDirections(
        H_sigma=H_sigma,
        H_beta=H_beta,
        H_R_bar=H_R_bar,
        b_damage=b_damage,
        regularization=regularization_value,
    )

# -*- coding: utf-8 -*-
"""
Descent search-direction operators for the one-dimensional LATIN solver.

The global stage uses the diagonal approximation

    H_minus = diag(H_sigma, H_beta, H_R_bar)

of the tangent of the viscoplastic dissipation potential. Following the
reference paper, the operator is regularized by

    H_minus <- H_minus + zeta * M^(-1),

with

    M = diag(E, C, R_inf)

in the present one-dimensional formulation and zeta = 0.15 by default.

The damage descent operator is

    b_minus = 0,

so damage is frozen during the linear global stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from numpy.typing import NDArray

from latin.local_stage import (
    isotropic_force_from_transformed_force,
)
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class DescentSearchDirections:
    """Diagonal LATIN descent operators on the full time-space grid."""

    H_sigma: FloatArray
    H_beta: FloatArray
    H_R_bar: FloatArray
    b_damage: FloatArray
    regularization: float

    def __post_init__(self) -> None:
        arrays = (
            np.asarray(self.H_sigma, dtype=np.float64),
            np.asarray(self.H_beta, dtype=np.float64),
            np.asarray(self.H_R_bar, dtype=np.float64),
            np.asarray(self.b_damage, dtype=np.float64),
        )

        reference_shape = arrays[0].shape

        if len(reference_shape) != 2:
            raise ValueError(
                "Search-direction fields must be two-dimensional."
            )

        for array in arrays:
            if array.shape != reference_shape:
                raise ValueError(
                    "All search-direction fields must have the same shape."
                )
            if not np.all(np.isfinite(array)):
                raise ValueError(
                    "Search-direction fields contain non-finite values."
                )

        if np.any(arrays[0] <= 0.0):
            raise ValueError("H_sigma must be strictly positive.")
        if np.any(arrays[1] <= 0.0):
            raise ValueError("H_beta must be strictly positive.")
        if np.any(arrays[2] <= 0.0):
            raise ValueError("H_R_bar must be strictly positive.")
        if np.any(arrays[3] != 0.0):
            raise ValueError("b_damage must be zero in this formulation.")
        if self.regularization <= 0.0:
            raise ValueError("regularization must be positive.")

        object.__setattr__(self, "H_sigma", arrays[0])
        object.__setattr__(self, "H_beta", arrays[1])
        object.__setattr__(self, "H_R_bar", arrays[2])
        object.__setattr__(self, "b_damage", arrays[3])

    @property
    def field_shape(self) -> tuple:
        """Common shape (n_time, n_elements)."""
        return self.H_sigma.shape


def _positive_part(value: FloatArray) -> FloatArray:
    """Return the componentwise positive part."""
    return np.maximum(value, 0.0)


def _yield_function(
    stress: FloatArray,
    beta: FloatArray,
    transformed_force: FloatArray,
    damage: FloatArray,
    material: MaterialParameters,
) -> FloatArray:
    """Evaluate the one-dimensional transformed yield function."""
    damage_safe = np.clip(
        damage,
        0.0,
        material.damage_upper_bound,
    )
    relative_effective_stress = (
        stress / (1.0 - damage_safe) - beta
    )

    isotropic_force = np.empty_like(transformed_force)

    for index, value in np.ndenumerate(transformed_force):
        isotropic_force[index] = (
            isotropic_force_from_transformed_force(
                transformed_force=float(value),
                material=material,
            )
        )

    return (
        np.abs(relative_effective_stress)
        + material.a * beta**2 / (2.0 * material.C)
        - isotropic_force
        - material.sigma_y
    )


def compute_descent_search_directions(
    local_state: LatinState,
    materials: Sequence[MaterialParameters],
    *,
    regularization: float = 0.15,
) -> DescentSearchDirections:
    """
    Compute the regularized diagonal descent operators.

    Parameters
    ----------
    local_state:
        State obtained from the nonlinear local stage.
    materials:
        One material parameter object per finite element.
    regularization:
        Dimensionless coefficient zeta in H <- H + zeta M^(-1).

    Returns
    -------
    DescentSearchDirections
        Full time-space arrays H_sigma, H_beta, H_R_bar and b_damage.
    """
    if len(materials) != local_state.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if regularization <= 0.0:
        raise ValueError("regularization must be positive.")

    shape = local_state.field_shape
    H_sigma = np.zeros(shape, dtype=np.float64)
    H_beta = np.zeros(shape, dtype=np.float64)
    H_R_bar = np.zeros(shape, dtype=np.float64)
    b_damage = np.zeros(shape, dtype=np.float64)

    for element, material in enumerate(materials):
        if material.E <= 0.0:
            raise ValueError("Young's modulus must be positive.")
        if material.C <= 0.0:
            raise ValueError(
                "Kinematic hardening modulus C must be positive."
            )
        if material.R_inf <= 0.0:
            raise ValueError("R_inf must be positive.")
        if material.gamma <= 0.0:
            raise ValueError("gamma must be positive.")

        stress = local_state.stress[:, element]
        beta = local_state.beta[:, element]
        transformed_force = local_state.R_bar[:, element]
        damage = np.clip(
            local_state.damage[:, element],
            0.0,
            material.damage_upper_bound,
        )

        relative_effective_stress = (
            stress / (1.0 - damage) - beta
        )
        flow_direction = np.sign(relative_effective_stress)

        yield_function = _yield_function(
            stress=stress,
            beta=beta,
            transformed_force=transformed_force,
            damage=damage,
            material=material,
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

        H_sigma[:, element] = (
            common_tangent_factor / (1.0 - damage) ** 2
            + regularization / material.E
        )

        beta_gradient = (
            -flow_direction
            + material.a * beta / material.C
        )
        H_beta[:, element] = (
            common_tangent_factor * beta_gradient**2
            + plastic_multiplier * material.a / material.C
            + regularization / material.C
        )

        transformed_gradient = (
            1.0
            - transformed_force
            * np.sqrt(material.gamma)
            / (2.0 * material.R_inf)
        )
        H_R_bar[:, element] = (
            common_tangent_factor
            * material.gamma
            * transformed_gradient**2
            + plastic_multiplier
            * material.gamma
            / (2.0 * material.R_inf)
            + regularization / material.R_inf
        )

    return DescentSearchDirections(
        H_sigma=H_sigma,
        H_beta=H_beta,
        H_R_bar=H_R_bar,
        b_damage=b_damage,
        regularization=float(regularization),
    )

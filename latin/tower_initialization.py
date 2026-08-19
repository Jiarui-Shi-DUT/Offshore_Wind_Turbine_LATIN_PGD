# -*- coding: utf-8 -*-
"""Elastic whole-history initialization for the tower LATIN-PGD solver."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from latin.tower_equilibrium_operator import TowerEquilibriumOperator
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import (
    MaterialParameters,
    damage_energy_release_rate,
)

FloatArray = NDArray[np.float64]
MaterialInput = Union[MaterialParameters, Sequence[MaterialParameters]]


def _readonly_float_array(values: FloatArray, name: str, ndim: int) -> FloatArray:
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != ndim:
        raise ValueError(name + " must have " + str(ndim) + " dimension(s).")
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain only finite values.")
    array.setflags(write=False)
    return array


def _positive_finite(value: float, name: str) -> float:
    result = float(value)
    if not np.isfinite(result) or result <= 0.0:
        raise ValueError(name + " must be positive and finite.")
    return result


def _material_sequence(
    materials: MaterialInput,
    n_material_points: int,
) -> Tuple[MaterialParameters, ...]:
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


@dataclass(frozen=True)
class TowerElasticInitialization:
    """Elastic LATIN state plus structural initialization diagnostics."""

    state: LatinStateTower
    displacements: FloatArray
    load_vectors: FloatArray
    free_equilibrium_residual: FloatArray
    stress_to_force_factor: float

    def __post_init__(self) -> None:
        if not isinstance(self.state, LatinStateTower):
            raise TypeError("state must be a LatinStateTower.")
        displacements = _readonly_float_array(
            self.displacements, "displacements", ndim=2
        )
        load_vectors = _readonly_float_array(
            self.load_vectors, "load_vectors", ndim=2
        )
        residual = _readonly_float_array(
            self.free_equilibrium_residual,
            "free_equilibrium_residual",
            ndim=2,
        )
        factor = _positive_finite(
            self.stress_to_force_factor, "stress_to_force_factor"
        )

        if displacements.shape != load_vectors.shape:
            raise ValueError(
                "displacements and load_vectors must have the same shape."
            )
        if displacements.shape[0] != self.state.n_time:
            raise ValueError(
                "Structural histories must contain one row per LATIN time."
            )
        if residual.shape[0] != self.state.n_time:
            raise ValueError(
                "free_equilibrium_residual must contain one row per time."
            )
        if residual.shape[1] > displacements.shape[1]:
            raise ValueError(
                "free_equilibrium_residual cannot contain more DOFs than "
                "the full structural history."
            )

        object.__setattr__(self, "state", self.state.copy())
        object.__setattr__(self, "displacements", displacements)
        object.__setattr__(self, "load_vectors", load_vectors)
        object.__setattr__(self, "free_equilibrium_residual", residual)
        object.__setattr__(self, "stress_to_force_factor", factor)

    @property
    def maximum_free_equilibrium_residual(self) -> float:
        if self.free_equilibrium_residual.size == 0:
            return 0.0
        return float(np.max(np.abs(self.free_equilibrium_residual)))


def compute_tower_elastic_initialization(
    time: FloatArray,
    load_vectors: FloatArray,
    materials: MaterialInput,
    equilibrium_operator: TowerEquilibriumOperator,
    stress_to_force_factor: float,
) -> TowerElasticInitialization:
    """Construct the globally admissible elastic whole-history state s0."""
    if not isinstance(equilibrium_operator, TowerEquilibriumOperator):
        raise TypeError(
            "equilibrium_operator must be a TowerEquilibriumOperator."
        )

    time_array = np.asarray(time, dtype=np.float64)
    if time_array.ndim != 1:
        raise ValueError("time must be one-dimensional.")
    if time_array.size < 2:
        raise ValueError("At least two time points are required.")
    if np.any(~np.isfinite(time_array)):
        raise ValueError("time must contain only finite values.")
    if np.any(np.diff(time_array) <= 0.0):
        raise ValueError("time must be strictly increasing.")

    loads = np.asarray(load_vectors, dtype=np.float64)
    expected_load_shape = (time_array.size, equilibrium_operator.n_dof)
    if loads.shape != expected_load_shape:
        raise ValueError(
            "load_vectors must have shape " + str(expected_load_shape) + "."
        )
    if np.any(~np.isfinite(loads)):
        raise ValueError("load_vectors must contain only finite values.")

    factor = _positive_finite(
        stress_to_force_factor, "stress_to_force_factor"
    )
    material_tuple = _material_sequence(
        materials, equilibrium_operator.n_material_points
    )

    material_modulus = np.array(
        [material.E for material in material_tuple], dtype=np.float64
    )
    if np.any(~np.isfinite(material_modulus)) or np.any(material_modulus <= 0.0):
        raise ValueError(
            "Every material-point elastic modulus must be positive and finite."
        )

    reference_modulus = equilibrium_operator.reference_modulus
    modulus_scale = max(
        1.0,
        float(np.max(np.abs(reference_modulus))),
        float(np.max(np.abs(material_modulus))),
    )
    if not np.allclose(
        material_modulus,
        reference_modulus,
        rtol=1.0e-12,
        atol=1.0e-12 * modulus_scale,
    ):
        raise ValueError(
            "materials and equilibrium_operator must use the same "
            "reference elastic modulus at every material point."
        )

    free_dofs = equilibrium_operator.free_dofs
    rhs = loads[:, free_dofs] / factor

    try:
        displacement_free = np.linalg.solve(
            equilibrium_operator.reduced_stiffness,
            rhs.T,
        ).T
    except np.linalg.LinAlgError as error:
        raise np.linalg.LinAlgError(
            "The reference tower stiffness could not be solved during "
            "elastic initialization."
        ) from error

    displacements = np.zeros(
        (time_array.size, equilibrium_operator.n_dof),
        dtype=np.float64,
    )
    displacements[:, free_dofs] = displacement_free

    elastic_strain = (
        displacement_free @ equilibrium_operator.compatibility_matrix.T
    )
    stress = (
        elastic_strain * reference_modulus[np.newaxis, :]
    )

    energy_release_rate = np.empty_like(stress)
    for q, material in enumerate(material_tuple):
        for time_index in range(time_array.size):
            energy_release_rate[time_index, q] = damage_energy_release_rate(
                stress=float(stress[time_index, q]),
                damage=0.0,
                material=material,
            )

    shape = (
        time_array.size,
        equilibrium_operator.n_material_points,
    )

    def zero_field() -> FloatArray:
        return np.zeros(shape, dtype=np.float64)

    state = LatinStateTower(
        time=time_array,
        plastic_strain_rate=zero_field(),
        elastic_strain=elastic_strain,
        alpha_rate=zero_field(),
        r_bar_rate=zero_field(),
        damage_rate=zero_field(),
        stress=stress,
        beta=zero_field(),
        R_bar=zero_field(),
        energy_release_rate=energy_release_rate,
        plastic_strain=zero_field(),
        alpha=zero_field(),
        r_bar=zero_field(),
        damage=zero_field(),
    )

    free_equilibrium_residual = (
        factor * equilibrium_operator.equilibrium_residual(state.stress)
        - loads[:, free_dofs]
    )

    return TowerElasticInitialization(
        state=state,
        displacements=displacements,
        load_vectors=loads,
        free_equilibrium_residual=free_equilibrium_residual,
        stress_to_force_factor=factor,
    )

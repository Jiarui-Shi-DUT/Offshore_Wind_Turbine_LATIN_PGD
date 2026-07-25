# -*- coding: utf-8 -*-
"""
Elastic initialization for the one-dimensional LATIN formulation.

The initial state satisfies:
    - the prescribed displacement boundary conditions,
    - global equilibrium,
    - the linear elastic constitutive law,
    - zero plastic, hardening and damage variables.

For a one-dimensional bar without distributed loads, the axial force is
constant along the bar. The elastic solution can therefore be obtained
directly from the series compliance of all elements.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class ElasticInitialization:
    """Complete elastic space-time initialization of the bar."""

    time: FloatArray
    displacement: FloatArray
    strain: FloatArray
    stress: FloatArray
    state: FloatArray
    reaction_left: FloatArray
    reaction_right: FloatArray


def _validate_inputs(
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    time: FloatArray,
    right_displacement: FloatArray,
) -> Tuple[FloatArray, FloatArray]:
    """Validate and normalize the elastic-initialization inputs."""
    if area <= 0.0:
        raise ValueError("area must be positive.")

    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )

    time_array = np.asarray(time, dtype=np.float64)
    displacement_array = np.asarray(
        right_displacement,
        dtype=np.float64,
    )

    if time_array.ndim != 1:
        raise ValueError("time must be a one-dimensional array.")
    if displacement_array.ndim != 1:
        raise ValueError(
            "right_displacement must be a one-dimensional array."
        )
    if time_array.size < 2:
        raise ValueError("At least two time points are required.")
    if displacement_array.shape != time_array.shape:
        raise ValueError(
            "right_displacement must have the same shape as time."
        )
    if np.any(np.diff(time_array) <= 0.0):
        raise ValueError("time must be strictly increasing.")
    if not np.all(np.isfinite(time_array)):
        raise ValueError("time contains non-finite values.")
    if not np.all(np.isfinite(displacement_array)):
        raise ValueError(
            "right_displacement contains non-finite values."
        )

    for element, material in enumerate(materials):
        if material.E <= 0.0:
            raise ValueError(
                f"Young's modulus must be positive in element {element}."
            )

    return time_array, displacement_array


def compute_elastic_initialization(
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    time: FloatArray,
    right_displacement: FloatArray,
) -> ElasticInitialization:
    """
    Compute the globally admissible elastic initialization.

    Boundary conditions:
        u(0, t) = 0
        u(L, t) = right_displacement(t)

    Assumptions:
        - quasi-static response,
        - no body force,
        - no distributed or interior nodal force,
        - one constant cross-sectional area,
        - one-dimensional two-node bar elements.

    The internal-variable state ordering is:
        [plastic_strain, alpha, r_bar, damage]
    and all four variables are initialized to zero.
    """
    time_array, prescribed_displacement = _validate_inputs(
        mesh=mesh,
        area=area,
        materials=materials,
        time=time,
        right_displacement=right_displacement,
    )

    element_lengths = mesh.element_lengths
    elastic_moduli = np.array(
        [material.E for material in materials],
        dtype=np.float64,
    )

    element_compliances = element_lengths / (
        area * elastic_moduli
    )
    total_compliance = float(np.sum(element_compliances))

    if not np.isfinite(total_compliance) or total_compliance <= 0.0:
        raise FloatingPointError(
            "The assembled elastic compliance is invalid."
        )

    axial_force = prescribed_displacement / total_compliance
    stress_scalar = axial_force / area

    stress = stress_scalar[:, np.newaxis] * np.ones(
        (1, mesh.n_elements),
        dtype=np.float64,
    )
    strain = stress / elastic_moduli[np.newaxis, :]

    displacement = np.zeros(
        (time_array.size, mesh.n_nodes),
        dtype=np.float64,
    )
    displacement[:, 1:] = np.cumsum(
        strain * element_lengths[np.newaxis, :],
        axis=1,
    )

    displacement[:, 0] = 0.0
    displacement[:, -1] = prescribed_displacement

    state = np.zeros(
        (time_array.size, mesh.n_elements, 4),
        dtype=np.float64,
    )

    reaction_left = -axial_force
    reaction_right = axial_force

    if not np.all(np.isfinite(displacement)):
        raise FloatingPointError(
            "Non-finite displacement in elastic initialization."
        )
    if not np.all(np.isfinite(strain)):
        raise FloatingPointError(
            "Non-finite strain in elastic initialization."
        )
    if not np.all(np.isfinite(stress)):
        raise FloatingPointError(
            "Non-finite stress in elastic initialization."
        )

    return ElasticInitialization(
        time=time_array.copy(),
        displacement=displacement,
        strain=strain,
        stress=stress,
        state=state,
        reaction_left=reaction_left,
        reaction_right=reaction_right,
    )

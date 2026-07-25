# -*- coding: utf-8 -*-
"""
One-dimensional equilibrium operator used by the LATIN global stage.

For a prescribed element field q(x, t), the operator E returns a compatible
strain field

    epsilon = E q

such that the associated correction stress

    sigma = C (epsilon - q)

is in equilibrium and the correction displacement is zero at both ends.

This is the discrete one-dimensional counterpart of Eqs. (50)-(51) in the
reference LATIN-PGD formulation. The same operator is later used for:

    - the plastic correction q = Delta epsilon_p,
    - the damage residual correction q = Delta epsilon_R.
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
class EquilibriumProjection:
    """Result of applying the homogeneous equilibrium operator."""

    source_strain: FloatArray
    compatible_strain: FloatArray
    stress: FloatArray
    displacement: FloatArray
    reaction_left: FloatArray
    reaction_right: FloatArray

    @property
    def n_time(self) -> int:
        return int(self.compatible_strain.shape[0])

    @property
    def n_elements(self) -> int:
        return int(self.compatible_strain.shape[1])


def _normalize_source_strain(
    source_strain: FloatArray,
    n_elements: int,
) -> Tuple[FloatArray, bool]:
    """
    Convert a one-dimensional spatial field or a time-space field to 2D.

    Returns
    -------
    normalized, was_spatial_only
    """
    source = np.asarray(source_strain, dtype=np.float64)

    if source.ndim == 1:
        if source.shape != (n_elements,):
            raise ValueError(
                "A spatial source_strain must have shape (n_elements,)."
            )
        source = source[np.newaxis, :]
        spatial_only = True
    elif source.ndim == 2:
        if source.shape[1] != n_elements:
            raise ValueError(
                "A time-space source_strain must have "
                "shape (n_time, n_elements)."
            )
        spatial_only = False
    else:
        raise ValueError(
            "source_strain must be a one- or two-dimensional array."
        )

    if source.shape[0] < 1:
        raise ValueError("source_strain must contain at least one field.")
    if not np.all(np.isfinite(source)):
        raise ValueError("source_strain contains non-finite values.")

    return source, spatial_only


def apply_equilibrium_operator(
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    source_strain: FloatArray,
) -> EquilibriumProjection:
    """
    Apply the one-dimensional homogeneous equilibrium operator E.

    The correction problem satisfies

        delta_u(0, t) = delta_u(L, t) = 0,

    and has no body force or interior nodal force. Therefore, the axial force
    is constant along the bar. For each time point,

        epsilon_e = q_e + sigma / E_e,

    while compatibility requires

        sum_e L_e epsilon_e = 0.

    These two relations determine the constant stress exactly.
    """
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )

    source, _ = _normalize_source_strain(
        source_strain=source_strain,
        n_elements=mesh.n_elements,
    )

    element_lengths = mesh.element_lengths
    elastic_moduli = np.array(
        [material.E for material in materials],
        dtype=np.float64,
    )

    if np.any(elastic_moduli <= 0.0):
        raise ValueError("All Young's moduli must be positive.")

    total_compliance = float(
        np.sum(element_lengths / (area * elastic_moduli))
    )
    if not np.isfinite(total_compliance) or total_compliance <= 0.0:
        raise FloatingPointError(
            "The assembled bar compliance is invalid."
        )

    source_extension = source @ element_lengths
    axial_force = -source_extension / total_compliance
    stress_scalar = axial_force / area

    stress = stress_scalar[:, np.newaxis] * np.ones(
        (1, mesh.n_elements),
        dtype=np.float64,
    )
    compatible_strain = (
        source + stress / elastic_moduli[np.newaxis, :]
    )

    displacement = np.zeros(
        (source.shape[0], mesh.n_nodes),
        dtype=np.float64,
    )
    displacement[:, 1:] = np.cumsum(
        compatible_strain * element_lengths[np.newaxis, :],
        axis=1,
    )

    # Remove only round-off at the homogeneous right boundary.
    right_boundary_error = displacement[:, -1].copy()
    if np.max(np.abs(right_boundary_error)) > 1.0e-10:
        raise FloatingPointError(
            "The equilibrium projection violates the homogeneous "
            "right-end displacement condition."
        )
    displacement[:, -1] = 0.0

    reaction_left = -axial_force
    reaction_right = axial_force

    if not np.all(np.isfinite(compatible_strain)):
        raise FloatingPointError(
            "Non-finite compatible strain in equilibrium projection."
        )
    if not np.all(np.isfinite(stress)):
        raise FloatingPointError(
            "Non-finite stress in equilibrium projection."
        )
    if not np.all(np.isfinite(displacement)):
        raise FloatingPointError(
            "Non-finite displacement in equilibrium projection."
        )

    return EquilibriumProjection(
        source_strain=source.copy(),
        compatible_strain=compatible_strain,
        stress=stress,
        displacement=displacement,
        reaction_left=reaction_left,
        reaction_right=reaction_right,
    )

# -*- coding: utf-8 -*-
"""
Tower LATIN-PGD global-stage preparation and unrelaxed finishing.

This module closes the global stage after a tower PGD plastic correction has
been obtained.  Its responsibilities are deliberately limited to two pure
operations:

1. prepare the same-iteration data that are frozen between Trial A and Trial B;
2. build one complete *unrelaxed* global candidate from a FixedBasisPGDResult.

The tower-v1 formulation follows the project freeze:

* the plastic reduced branch uses the paper-separated forcing

      f = eps_p_dot_hat - eps_p_dot_i
          - H_sigma (sigma_hat - sigma_i),

  i.e. the damage stress correction is not folded into the PGD forcing;
* the nonlinear damaged-elastic mismatch is handled by a separate full-order
  reference-equilibrium projection;
* hardening is completed pointwise with Eq. (73)-(74) and backward Euler;
* b_minus = 0 gives D_dot_candidate = D_dot_hat, and tower v1 also inherits
  the already-integrated local damage history D_candidate = D_hat;
* the final energy-release rate is recomputed only after the candidate stress
  and damage histories are known.

No relaxation, LATIN indicator, saturation indicator, enrichment decision, or
persistent solver commit is performed here.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple, Union

import numpy as np
from numpy.typing import NDArray

from latin.search_directions import DescentSearchDirections
from latin.tower_equilibrium_operator import (
    EquilibriumProjectionTower,
    TowerEquilibriumOperator,
)
from latin.tower_pgd_time_update import FixedBasisPGDResult
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]
MaterialInput = Union[
    MaterialParameters,
    Sequence[MaterialParameters],
]


def _readonly_float_array(
    values: FloatArray,
    name: str,
    ndim: int,
) -> FloatArray:
    """Return a detached finite read-only float64 array."""
    array = np.array(values, dtype=np.float64, copy=True)
    if array.ndim != ndim:
        raise ValueError(
            name
            + " must have "
            + str(ndim)
            + " dimension(s)."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain only finite values.")
    array.setflags(write=False)
    return array


def _state_pair_validation(
    baseline_state: LatinStateTower,
    local_state: LatinStateTower,
) -> None:
    """Validate a common immutable tower time-material grid."""
    if not isinstance(baseline_state, LatinStateTower):
        raise TypeError(
            "baseline_state must be a LatinStateTower."
        )
    if not isinstance(local_state, LatinStateTower):
        raise TypeError(
            "local_state must be a LatinStateTower."
        )
    if baseline_state.field_shape != local_state.field_shape:
        raise ValueError(
            "baseline_state and local_state must have the same field shape."
        )
    if not np.array_equal(
        baseline_state.time,
        local_state.time,
    ):
        raise ValueError(
            "baseline_state and local_state must use the same time grid."
        )


def _search_direction_validation(
    directions: DescentSearchDirections,
    field_shape: Tuple[int, int],
) -> None:
    """Validate the reused diagonal descent-direction container."""
    if not isinstance(directions, DescentSearchDirections):
        raise TypeError(
            "directions must be a DescentSearchDirections object."
        )
    if directions.field_shape != field_shape:
        raise ValueError(
            "Search directions must use the tower state field shape."
        )
    if np.any(directions.b_damage != 0.0):
        raise ValueError(
            "tower v1 requires b_damage = 0."
        )


def _material_parameter_arrays(
    materials: MaterialInput,
    n_material_points: int,
) -> Tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]:
    """Broadcast one material or validate one material object per q-point."""
    if isinstance(materials, MaterialParameters):
        material_list = (materials,) * n_material_points
    else:
        material_list = tuple(materials)
        if len(material_list) == 1:
            material_list = material_list * n_material_points

    if len(material_list) != n_material_points:
        raise ValueError(
            "materials must contain one MaterialParameters object or one "
            "object per material point."
        )
    if any(
        not isinstance(material, MaterialParameters)
        for material in material_list
    ):
        raise TypeError(
            "materials must contain only MaterialParameters objects."
        )

    elastic_modulus = np.array(
        [material.E for material in material_list],
        dtype=np.float64,
    )
    kinematic_modulus = np.array(
        [material.C for material in material_list],
        dtype=np.float64,
    )
    isotropic_modulus = np.array(
        [material.R_inf for material in material_list],
        dtype=np.float64,
    )
    closure_parameter = np.array(
        [material.h for material in material_list],
        dtype=np.float64,
    )
    damage_upper_bound = np.array(
        [material.damage_upper_bound for material in material_list],
        dtype=np.float64,
    )

    for name, array in (
        ("elastic modulus", elastic_modulus),
        ("kinematic hardening modulus", kinematic_modulus),
        ("isotropic hardening modulus", isotropic_modulus),
        ("closure parameter", closure_parameter),
        ("damage upper bound", damage_upper_bound),
    ):
        if np.any(~np.isfinite(array)):
            raise ValueError(name + " must be finite at every material point.")

    if np.any(elastic_modulus <= 0.0):
        raise ValueError(
            "Every material-point elastic modulus must be positive."
        )
    if np.any(kinematic_modulus <= 0.0):
        raise ValueError(
            "Every material-point kinematic hardening modulus must be positive."
        )
    if np.any(isotropic_modulus <= 0.0):
        raise ValueError(
            "Every material-point isotropic hardening modulus must be positive."
        )
    if np.any(closure_parameter < 0.0) or np.any(closure_parameter > 1.0):
        raise ValueError(
            "Every material-point closure parameter must satisfy 0 <= h <= 1."
        )
    if (
        np.any(damage_upper_bound <= 0.0)
        or np.any(damage_upper_bound >= 1.0)
    ):
        raise ValueError(
            "Every damage upper bound must satisfy 0 < bound < 1."
        )

    return (
        elastic_modulus,
        kinematic_modulus,
        isotropic_modulus,
        closure_parameter,
        damage_upper_bound,
    )


@dataclass(frozen=True)
class FrozenGlobalData:
    """Same-iteration tower data shared by fixed-basis and enriched trials."""

    damage_residual_strain: FloatArray
    damage_projection: EquilibriumProjectionTower
    full_plastic_forcing: FloatArray
    reference_modulus: FloatArray

    def __post_init__(self) -> None:
        residual = _readonly_float_array(
            self.damage_residual_strain,
            "damage_residual_strain",
            ndim=2,
        )
        forcing = _readonly_float_array(
            self.full_plastic_forcing,
            "full_plastic_forcing",
            ndim=2,
        )
        reference_modulus = _readonly_float_array(
            self.reference_modulus,
            "reference_modulus",
            ndim=1,
        )

        if residual.shape != forcing.shape:
            raise ValueError(
                "damage_residual_strain and full_plastic_forcing must have "
                "the same shape."
            )
        if reference_modulus.shape != (residual.shape[1],):
            raise ValueError(
                "reference_modulus must contain one value per material point."
            )
        if np.any(reference_modulus <= 0.0):
            raise ValueError(
                "reference_modulus must be strictly positive."
            )
        if not isinstance(
            self.damage_projection,
            EquilibriumProjectionTower,
        ):
            raise TypeError(
                "damage_projection must be an EquilibriumProjectionTower."
            )
        if not self.damage_projection.is_history:
            raise ValueError(
                "damage_projection must be a time-history projection."
            )
        if self.damage_projection.source_strain.shape != residual.shape:
            raise ValueError(
                "damage_projection must use the frozen field shape."
            )

        scale = max(
            1.0,
            float(np.linalg.norm(residual)),
        )
        if not np.allclose(
            self.damage_projection.source_strain,
            residual,
            rtol=0.0,
            atol=1.0e-12 * scale,
        ):
            raise ValueError(
                "damage_projection source strain must equal "
                "damage_residual_strain."
            )

        object.__setattr__(
            self,
            "damage_residual_strain",
            residual,
        )
        object.__setattr__(
            self,
            "full_plastic_forcing",
            forcing,
        )
        object.__setattr__(
            self,
            "reference_modulus",
            reference_modulus,
        )

    @property
    def field_shape(self) -> Tuple[int, int]:
        """Return common (n_time, n_material_points) shape."""
        return (
            int(self.damage_residual_strain.shape[0]),
            int(self.damage_residual_strain.shape[1]),
        )

    @property
    def damage_displacement_correction(self) -> FloatArray:
        """Return a defensive copy of the frozen full-order damage displacement."""
        return self.damage_projection.displacement_free.copy()


@dataclass(frozen=True)
class TowerGlobalCandidate:
    """Complete unrelaxed tower global-stage candidate."""

    state: LatinStateTower
    total_displacement_correction: FloatArray
    plastic_displacement_correction: FloatArray
    damage_displacement_correction: FloatArray

    def __post_init__(self) -> None:
        if not isinstance(self.state, LatinStateTower):
            raise TypeError(
                "state must be a LatinStateTower."
            )

        total = _readonly_float_array(
            self.total_displacement_correction,
            "total_displacement_correction",
            ndim=2,
        )
        plastic = _readonly_float_array(
            self.plastic_displacement_correction,
            "plastic_displacement_correction",
            ndim=2,
        )
        damage = _readonly_float_array(
            self.damage_displacement_correction,
            "damage_displacement_correction",
            ndim=2,
        )

        if plastic.shape != damage.shape or total.shape != plastic.shape:
            raise ValueError(
                "All displacement-correction arrays must have the same shape."
            )
        if total.shape[0] != self.state.n_time:
            raise ValueError(
                "Displacement corrections must contain one row per time point."
            )
        if not np.allclose(
            total,
            plastic + damage,
            rtol=0.0,
            atol=1.0e-12
            * max(
                1.0,
                float(np.linalg.norm(total)),
            ),
        ):
            raise ValueError(
                "total_displacement_correction must equal plastic + damage."
            )

        object.__setattr__(
            self,
            "total_displacement_correction",
            total,
        )
        object.__setattr__(
            self,
            "plastic_displacement_correction",
            plastic,
        )
        object.__setattr__(
            self,
            "damage_displacement_correction",
            damage,
        )


def prepare_frozen_global_data(
    baseline_state: LatinStateTower,
    local_state: LatinStateTower,
    directions: DescentSearchDirections,
    equilibrium_operator: TowerEquilibriumOperator,
) -> FrozenGlobalData:
    """
    Prepare damage projection and paper-separated plastic forcing once.

    The residual-stress branch is

        Delta R_res = (sigma_i - sigma_hat)
                      - C0 (eps_e_i - eps_e_hat),

    followed by Delta eps_R = C0^(-1) Delta R_res and the full tower reference
    projection.  The PGD forcing intentionally remains paper-separated:

        f = eps_p_dot_hat - eps_p_dot_i
            - H_sigma (sigma_hat - sigma_i).
    """
    _state_pair_validation(
        baseline_state=baseline_state,
        local_state=local_state,
    )
    _search_direction_validation(
        directions=directions,
        field_shape=baseline_state.field_shape,
    )
    if not isinstance(
        equilibrium_operator,
        TowerEquilibriumOperator,
    ):
        raise TypeError(
            "equilibrium_operator must be a TowerEquilibriumOperator."
        )
    if (
        equilibrium_operator.n_material_points
        != baseline_state.n_material_points
    ):
        raise ValueError(
            "equilibrium_operator and tower states must use the same "
            "material-point count."
        )

    reference_modulus = (
        equilibrium_operator.reference_modulus[np.newaxis, :]
    )
    residual_stress = (
        baseline_state.stress
        - local_state.stress
        - reference_modulus
        * (
            baseline_state.elastic_strain
            - local_state.elastic_strain
        )
    )
    damage_residual_strain = (
        residual_stress / reference_modulus
    )
    damage_projection = equilibrium_operator.apply_history(
        damage_residual_strain
    )

    full_plastic_forcing = (
        local_state.plastic_strain_rate
        - baseline_state.plastic_strain_rate
        - directions.H_sigma
        * (
            local_state.stress
            - baseline_state.stress
        )
    )

    return FrozenGlobalData(
        damage_residual_strain=damage_residual_strain,
        damage_projection=damage_projection,
        full_plastic_forcing=full_plastic_forcing,
        reference_modulus=equilibrium_operator.reference_modulus,
    )


def _validate_fixed_basis_result(
    fixed_basis_result: FixedBasisPGDResult,
    frozen_data: FrozenGlobalData,
    directions: DescentSearchDirections,
    field_shape: Tuple[int, int],
) -> None:
    """Ensure that one reduced result belongs to the same frozen trial data."""
    if not isinstance(
        fixed_basis_result,
        FixedBasisPGDResult,
    ):
        raise TypeError(
            "fixed_basis_result must be a FixedBasisPGDResult."
        )
    if fixed_basis_result.basis.field_shape != field_shape:
        raise ValueError(
            "fixed_basis_result must use the tower state field shape."
        )
    if frozen_data.field_shape != field_shape:
        raise ValueError(
            "frozen_data must use the tower state field shape."
        )

    reconstructed_residual = (
        fixed_basis_result.plastic_strain_rate_correction
        - directions.H_sigma
        * fixed_basis_result.plastic_projection.stress
        - frozen_data.full_plastic_forcing
    )
    scale = max(
        1.0,
        float(np.linalg.norm(reconstructed_residual)),
        float(np.linalg.norm(fixed_basis_result.mechanical_residual)),
    )
    if not np.allclose(
        reconstructed_residual,
        fixed_basis_result.mechanical_residual,
        rtol=0.0,
        atol=1.0e-10 * scale,
    ):
        raise ValueError(
            "fixed_basis_result was not generated from the supplied "
            "FrozenGlobalData full plastic forcing."
        )


def _finish_hardening(
    local_state: LatinStateTower,
    directions: DescentSearchDirections,
    kinematic_modulus: FloatArray,
    isotropic_modulus: FloatArray,
) -> Tuple[
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
]:
    """Complete Eq. (73)-(74) pointwise with tower-v1 backward Euler."""
    shape = local_state.field_shape
    alpha = np.zeros(shape, dtype=np.float64)
    alpha_rate = np.zeros(shape, dtype=np.float64)
    beta = np.zeros(shape, dtype=np.float64)
    r_bar = np.zeros(shape, dtype=np.float64)
    r_bar_rate = np.zeros(shape, dtype=np.float64)
    R_bar = np.zeros(shape, dtype=np.float64)

    alpha[0, :] = local_state.alpha[0, :]
    beta[0, :] = kinematic_modulus * alpha[0, :]
    alpha_rate[0, :] = (
        local_state.alpha_rate[0, :]
        - directions.H_beta[0, :]
        * (
            beta[0, :]
            - local_state.beta[0, :]
        )
    )

    r_bar[0, :] = local_state.r_bar[0, :]
    R_bar[0, :] = isotropic_modulus * r_bar[0, :]
    r_bar_rate[0, :] = (
        local_state.r_bar_rate[0, :]
        - directions.H_R_bar[0, :]
        * (
            R_bar[0, :]
            - local_state.R_bar[0, :]
        )
    )

    for step in range(1, local_state.n_time):
        dt = float(
            local_state.time[step]
            - local_state.time[step - 1]
        )

        H_beta = directions.H_beta[step, :]
        alpha_rhs = (
            local_state.alpha_rate[step, :]
            + H_beta * local_state.beta[step, :]
        )
        alpha_denominator = (
            1.0
            + dt
            * H_beta
            * kinematic_modulus
        )
        alpha[step, :] = (
            alpha[step - 1, :]
            + dt * alpha_rhs
        ) / alpha_denominator
        alpha_rate[step, :] = (
            alpha[step, :]
            - alpha[step - 1, :]
        ) / dt
        beta[step, :] = (
            kinematic_modulus * alpha[step, :]
        )

        H_R_bar = directions.H_R_bar[step, :]
        r_bar_rhs = (
            local_state.r_bar_rate[step, :]
            + H_R_bar * local_state.R_bar[step, :]
        )
        r_bar_denominator = (
            1.0
            + dt
            * H_R_bar
            * isotropic_modulus
        )
        r_bar[step, :] = (
            r_bar[step - 1, :]
            + dt * r_bar_rhs
        ) / r_bar_denominator
        r_bar_rate[step, :] = (
            r_bar[step, :]
            - r_bar[step - 1, :]
        ) / dt
        R_bar[step, :] = (
            isotropic_modulus * r_bar[step, :]
        )

    for name, array in (
        ("alpha", alpha),
        ("alpha_rate", alpha_rate),
        ("beta", beta),
        ("r_bar", r_bar),
        ("r_bar_rate", r_bar_rate),
        ("R_bar", R_bar),
    ):
        if np.any(~np.isfinite(array)):
            raise FloatingPointError(
                "Non-finite " + name + " in tower hardening finishing."
            )

    return (
        alpha,
        alpha_rate,
        beta,
        r_bar,
        r_bar_rate,
        R_bar,
    )


def _energy_release_rate(
    stress: FloatArray,
    damage: FloatArray,
    elastic_modulus: FloatArray,
    closure_parameter: FloatArray,
    damage_upper_bound: FloatArray,
) -> FloatArray:
    """Evaluate the scalar-fiber unilateral damage state law after mechanics."""
    if np.any(damage < 0.0):
        raise ValueError(
            "Candidate damage must be non-negative."
        )
    if np.any(
        damage
        > damage_upper_bound[np.newaxis, :]
        + 1.0e-14
    ):
        raise ValueError(
            "Candidate damage exceeds a material-point damage upper bound."
        )

    tensile_denominator = (
        1.0 - damage
    )
    compressive_denominator = (
        1.0
        - closure_parameter[np.newaxis, :]
        * damage
    )
    if np.any(tensile_denominator <= 0.0):
        raise ValueError(
            "Tensile damage denominator must remain positive."
        )
    if np.any(compressive_denominator <= 0.0):
        raise ValueError(
            "Compressive damage denominator must remain positive."
        )

    tensile = stress >= 0.0
    energy = np.where(
        tensile,
        stress**2
        / (
            2.0
            * elastic_modulus[np.newaxis, :]
            * tensile_denominator**2
        ),
        closure_parameter[np.newaxis, :]
        * stress**2
        / (
            2.0
            * elastic_modulus[np.newaxis, :]
            * compressive_denominator**2
        ),
    )
    if np.any(~np.isfinite(energy)):
        raise FloatingPointError(
            "Non-finite energy-release rate in tower global finishing."
        )
    return np.asarray(energy, dtype=np.float64)


def build_unrelaxed_candidate(
    baseline_state: LatinStateTower,
    local_state: LatinStateTower,
    directions: DescentSearchDirections,
    frozen_data: FrozenGlobalData,
    fixed_basis_result: FixedBasisPGDResult,
    materials: MaterialInput,
) -> TowerGlobalCandidate:
    """Build the complete tower global candidate before relaxation."""
    _state_pair_validation(
        baseline_state=baseline_state,
        local_state=local_state,
    )
    _search_direction_validation(
        directions=directions,
        field_shape=baseline_state.field_shape,
    )
    if not isinstance(frozen_data, FrozenGlobalData):
        raise TypeError(
            "frozen_data must be a FrozenGlobalData object."
        )
    _validate_fixed_basis_result(
        fixed_basis_result=fixed_basis_result,
        frozen_data=frozen_data,
        directions=directions,
        field_shape=baseline_state.field_shape,
    )

    (
        elastic_modulus,
        kinematic_modulus,
        isotropic_modulus,
        closure_parameter,
        damage_upper_bound,
    ) = _material_parameter_arrays(
        materials=materials,
        n_material_points=baseline_state.n_material_points,
    )
    modulus_scale = max(
        1.0,
        float(np.max(np.abs(frozen_data.reference_modulus))),
    )
    if not np.allclose(
        elastic_modulus,
        frozen_data.reference_modulus,
        rtol=0.0,
        atol=1.0e-12 * modulus_scale,
    ):
        raise ValueError(
            "materials.E must match the reference modulus used by the "
            "frozen tower equilibrium operator."
        )

    plastic_projection = fixed_basis_result.plastic_projection
    damage_projection = frozen_data.damage_projection

    plastic_strain = (
        baseline_state.plastic_strain
        + fixed_basis_result.plastic_strain_correction
    )
    plastic_strain_rate = (
        baseline_state.plastic_strain_rate
        + fixed_basis_result.plastic_strain_rate_correction
    )
    stress = (
        baseline_state.stress
        + plastic_projection.stress
        + damage_projection.stress
    )
    elastic_strain = (
        baseline_state.elastic_strain
        + plastic_projection.compatible_strain
        - fixed_basis_result.plastic_strain_correction
        + damage_projection.compatible_strain
    )

    (
        alpha,
        alpha_rate,
        beta,
        r_bar,
        r_bar_rate,
        R_bar,
    ) = _finish_hardening(
        local_state=local_state,
        directions=directions,
        kinematic_modulus=kinematic_modulus,
        isotropic_modulus=isotropic_modulus,
    )

    # Paper Eq. (75): D_dot is inherited.  Tower v1 additionally inherits the
    # already-integrated local D history instead of re-integrating nodal rates.
    damage_rate = np.asarray(
        local_state.damage_rate,
        dtype=np.float64,
    ).copy()
    damage = np.asarray(
        local_state.damage,
        dtype=np.float64,
    ).copy()

    energy_release_rate = _energy_release_rate(
        stress=stress,
        damage=damage,
        elastic_modulus=elastic_modulus,
        closure_parameter=closure_parameter,
        damage_upper_bound=damage_upper_bound,
    )

    state = LatinStateTower(
        time=baseline_state.time,
        plastic_strain_rate=plastic_strain_rate,
        elastic_strain=elastic_strain,
        alpha_rate=alpha_rate,
        r_bar_rate=r_bar_rate,
        damage_rate=damage_rate,
        stress=stress,
        beta=beta,
        R_bar=R_bar,
        energy_release_rate=energy_release_rate,
        plastic_strain=plastic_strain,
        alpha=alpha,
        r_bar=r_bar,
        damage=damage,
    )

    plastic_displacement = (
        plastic_projection.displacement_free
    )
    damage_displacement = (
        damage_projection.displacement_free
    )
    total_displacement = (
        plastic_displacement
        + damage_displacement
    )

    return TowerGlobalCandidate(
        state=state,
        total_displacement_correction=total_displacement,
        plastic_displacement_correction=plastic_displacement,
        damage_displacement_correction=damage_displacement,
    )

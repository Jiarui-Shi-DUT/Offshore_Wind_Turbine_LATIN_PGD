# -*- coding: utf-8 -*-
"""
Reduced LATIN global stage based on a reusable one-dimensional PGD basis.

The plastic correction is represented by separated space-time modes.  At each
global stage, the existing spatial basis is reused first and all temporal
functions are updated.  When the resulting mechanical residual remains above
the requested tolerance, residual-driven modes are added greedily.

The damage-dependent correction is not reduced.  It is computed by the same
equilibrium projection as in the validated full-order global stage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from latin.equilibrium_operator import (
    EquilibriumProjection,
    apply_equilibrium_operator,
)
from latin.global_stage import (
    _compute_damage_residual_strain,
    _update_damage_and_energy,
    _update_hardening_variables,
    _validate_global_stage_inputs,
)
from latin.pgd_basis import PGDBasis1D
from latin.pgd_enrichment import (
    PGDEnrichmentResult,
    enrich_pgd_basis_once,
)
from latin.pgd_time_update import update_pgd_time_functions
from latin.search_directions import DescentSearchDirections
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class PGDGlobalStageResult:
    """State, reduced basis and diagnostics produced by one PGD global stage."""

    state: LatinState
    basis: PGDBasis1D
    forcing: FloatArray
    mechanical_residual: FloatArray
    relative_residual: float
    residual_history: FloatArray
    modes_added: int
    time_functions_updated: bool
    reduced_converged: bool
    enrichment_results: Sequence[PGDEnrichmentResult]
    plastic_strain_correction: FloatArray
    plastic_strain_rate_correction: FloatArray
    plastic_projection: EquilibriumProjection
    residual_strain: FloatArray
    damage_projection: EquilibriumProjection
    displacement_correction: FloatArray

    @property
    def n_modes(self) -> int:
        """Dimension of the PGD basis returned by the global stage."""
        return self.basis.n_modes


def _validate_pgd_options(
    basis: PGDBasis1D,
    global_state: LatinState,
    reduced_tolerance: float,
    max_new_modes: int,
    fixed_point_tolerance: float,
    max_fixed_point_iterations: int,
    minimum_spatial_norm: float,
    acceptance_tolerance: float,
    rcond: float,
) -> None:
    """Validate PGD dimensions and numerical control parameters."""
    if basis.field_shape != global_state.field_shape:
        raise ValueError(
            "basis and LATIN states must use the same space-time grid."
        )
    if reduced_tolerance <= 0.0 or not np.isfinite(
        reduced_tolerance
    ):
        raise ValueError(
            "reduced_tolerance must be positive and finite."
        )
    if max_new_modes < 0:
        raise ValueError("max_new_modes must be non-negative.")
    if fixed_point_tolerance <= 0.0 or not np.isfinite(
        fixed_point_tolerance
    ):
        raise ValueError(
            "fixed_point_tolerance must be positive and finite."
        )
    if max_fixed_point_iterations < 1:
        raise ValueError(
            "max_fixed_point_iterations must be at least one."
        )
    if minimum_spatial_norm <= 0.0 or not np.isfinite(
        minimum_spatial_norm
    ):
        raise ValueError(
            "minimum_spatial_norm must be positive and finite."
        )
    if acceptance_tolerance < 0.0 or not np.isfinite(
        acceptance_tolerance
    ):
        raise ValueError(
            "acceptance_tolerance must be non-negative and finite."
        )
    if rcond <= 0.0 or not np.isfinite(rcond):
        raise ValueError("rcond must be positive and finite.")


def _compute_plastic_forcing(
    global_state: LatinState,
    local_state: LatinState,
    directions: DescentSearchDirections,
    damage_stress_correction: FloatArray,
) -> FloatArray:
    """
    Return the known right-hand side of the plastic correction equation.

    The reduced correction satisfies

        Delta eps_p_dot
        - H_sigma Delta sigma_plastic
        = forcing.
    """
    forcing = (
        local_state.plastic_strain_rate
        - global_state.plastic_strain_rate
        - directions.H_sigma
        * (
            local_state.stress
            - global_state.stress
        )
        + directions.H_sigma
        * damage_stress_correction
    )

    if not np.all(np.isfinite(forcing)):
        raise FloatingPointError(
            "The PGD plastic forcing contains non-finite values."
        )

    return forcing


def _weighted_space_time_norm(
    field: FloatArray,
    time: FloatArray,
    directions: DescentSearchDirections,
    element_volumes: FloatArray,
) -> float:
    """Evaluate the H_sigma-inverse weighted residual norm."""
    density = field**2 / directions.H_sigma
    space_integral = density @ element_volumes
    squared_norm = float(
        np.trapz(space_integral, x=time)
    )

    if squared_norm < -1.0e-12:
        raise FloatingPointError(
            "The squared PGD mechanical norm became negative."
        )

    return float(
        np.sqrt(max(squared_norm, 0.0))
    )


def _relative_mechanical_residual(
    residual: FloatArray,
    forcing: FloatArray,
    time: FloatArray,
    directions: DescentSearchDirections,
    element_volumes: FloatArray,
) -> float:
    """Return the reduced mechanical residual relative to the forcing norm."""
    residual_norm = _weighted_space_time_norm(
        field=residual,
        time=time,
        directions=directions,
        element_volumes=element_volumes,
    )
    forcing_norm = _weighted_space_time_norm(
        field=forcing,
        time=time,
        directions=directions,
        element_volumes=element_volumes,
    )

    numerical_zero = float(
        np.finfo(np.float64).eps
    )
    if forcing_norm <= numerical_zero:
        return (
            0.0
            if residual_norm <= numerical_zero
            else np.inf
        )

    return float(residual_norm / forcing_norm)


def solve_pgd_global_stage(
    global_state: LatinState,
    local_state: LatinState,
    directions: DescentSearchDirections,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    basis: PGDBasis1D,
    *,
    reduced_tolerance: float = 1.0e-4,
    max_new_modes: int = 5,
    allow_enrichment: bool = True,
    latin_iteration: int = 0,
    fixed_point_tolerance: float = 1.0e-6,
    max_fixed_point_iterations: int = 30,
    minimum_spatial_norm: float = 1.0e-14,
    acceptance_tolerance: float = 1.0e-8,
    rcond: float = 1.0e-12,
) -> PGDGlobalStageResult:
    """
    Perform one reduced LATIN global stage with adaptive PGD enrichment.

    Existing spatial modes are always reused first by updating their temporal
    functions.  New modes are generated only when ``allow_enrichment`` is true
    and the reduced mechanical residual remains above ``reduced_tolerance``.

    The outer LATIN solver will later decide ``allow_enrichment`` from the
    saturation indicator.  Keeping that decision outside this function avoids
    mixing the global-stage residual with the LATIN convergence indicator.
    """
    _validate_global_stage_inputs(
        global_state=global_state,
        local_state=local_state,
        directions=directions,
        mesh=mesh,
        area=area,
        materials=materials,
    )
    _validate_pgd_options(
        basis=basis,
        global_state=global_state,
        reduced_tolerance=reduced_tolerance,
        max_new_modes=max_new_modes,
        fixed_point_tolerance=fixed_point_tolerance,
        max_fixed_point_iterations=max_fixed_point_iterations,
        minimum_spatial_norm=minimum_spatial_norm,
        acceptance_tolerance=acceptance_tolerance,
        rcond=rcond,
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

    forcing = _compute_plastic_forcing(
        global_state=global_state,
        local_state=local_state,
        directions=directions,
        damage_stress_correction=damage_projection.stress,
    )

    working_basis = basis.copy()
    element_volumes = area * mesh.element_lengths
    time_functions_updated = False

    if working_basis.n_modes > 0:
        time_update = update_pgd_time_functions(
            basis=working_basis,
            time=global_state.time,
            directions=directions,
            forcing=forcing,
            mesh=mesh,
            area=area,
            rcond=rcond,
        )
        working_basis = time_update.basis
        mechanical_residual = time_update.residual
        time_functions_updated = True
    else:
        mechanical_residual = -forcing.copy()

    relative_residual = _relative_mechanical_residual(
        residual=mechanical_residual,
        forcing=forcing,
        time=global_state.time,
        directions=directions,
        element_volumes=element_volumes,
    )
    residual_values: List[float] = [
        relative_residual
    ]
    enrichment_values: List[PGDEnrichmentResult] = []

    while (
        allow_enrichment
        and relative_residual > reduced_tolerance
        and len(enrichment_values) < max_new_modes
    ):
        try:
            enrichment = enrich_pgd_basis_once(
                basis=working_basis,
                time=global_state.time,
                residual=mechanical_residual,
                directions=directions,
                mesh=mesh,
                area=area,
                materials=materials,
                iteration_added=latin_iteration,
                fixed_point_tolerance=(
                    fixed_point_tolerance
                ),
                max_fixed_point_iterations=(
                    max_fixed_point_iterations
                ),
                minimum_spatial_norm=minimum_spatial_norm,
                acceptance_tolerance=acceptance_tolerance,
                rcond=rcond,
            )
        except FloatingPointError:
            break

        enrichment_values.append(enrichment)

        if not enrichment.accepted:
            break

        working_basis.append(enrichment.mode)

        # Re-optimise all temporal coefficients jointly after enrichment.
        time_update = update_pgd_time_functions(
            basis=working_basis,
            time=global_state.time,
            directions=directions,
            forcing=forcing,
            mesh=mesh,
            area=area,
            rcond=rcond,
        )
        working_basis = time_update.basis
        mechanical_residual = time_update.residual
        time_functions_updated = True

        relative_residual = _relative_mechanical_residual(
            residual=mechanical_residual,
            forcing=forcing,
            time=global_state.time,
            directions=directions,
            element_volumes=element_volumes,
        )
        residual_values.append(relative_residual)

    plastic_strain_correction = (
        working_basis.plastic_strain_correction()
    )
    plastic_strain_rate_correction = (
        working_basis.plastic_strain_rate_correction()
    )
    plastic_projection = apply_equilibrium_operator(
        mesh=mesh,
        area=area,
        materials=materials,
        source_strain=plastic_strain_correction,
    )

    # Recompute the final residual from the fields actually used in the state.
    mechanical_residual = (
        plastic_strain_rate_correction
        - directions.H_sigma
        * plastic_projection.stress
        - forcing
    )
    relative_residual = _relative_mechanical_residual(
        residual=mechanical_residual,
        forcing=forcing,
        time=global_state.time,
        directions=directions,
        element_volumes=element_volumes,
    )
    if (
        not residual_values
        or abs(
            residual_values[-1]
            - relative_residual
        ) > 1.0e-15
    ):
        residual_values.append(relative_residual)

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

    return PGDGlobalStageResult(
        state=new_state,
        basis=working_basis,
        forcing=forcing,
        mechanical_residual=mechanical_residual,
        relative_residual=relative_residual,
        residual_history=np.asarray(
            residual_values,
            dtype=np.float64,
        ),
        modes_added=sum(
            1
            for enrichment in enrichment_values
            if enrichment.accepted
        ),
        time_functions_updated=time_functions_updated,
        reduced_converged=bool(
            np.isfinite(relative_residual)
            and relative_residual <= reduced_tolerance
        ),
        enrichment_results=tuple(enrichment_values),
        plastic_strain_correction=plastic_strain_correction,
        plastic_strain_rate_correction=(
            plastic_strain_rate_correction
        ),
        plastic_projection=plastic_projection,
        residual_strain=residual_strain,
        damage_projection=damage_projection,
        displacement_correction=displacement_correction,
    )

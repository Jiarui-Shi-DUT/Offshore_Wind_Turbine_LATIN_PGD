# -*- coding: utf-8 -*-
"""
Nonlinear asymmetric cyclic response of the NREL 5 MW tower.

The periodic tower-top force is

    F(t) = F_mean + F_amplitude * sin(2*pi*t/T),

with R_F = F_min / F_max. The formal benchmark uses R_F = -0.5.

Because the periodic history begins at F_mean rather than zero, this driver
first applies a separate linear preload 0 -> F_mean over one quarter period
using the same time increment as the cyclic history. The preload endpoint is
reused as the first stored periodic state, so cycle diagnostics contain only
complete periodic cycles.
"""

from __future__ import annotations

from typing import List

import numpy as np

from examples.elastic_tapered_tower import (
    TowerConfiguration,
    create_tower_geometry,
)
from examples.nonlinear_tower_reversed_response import (
    NonlinearCyclicResponse,
    locate_critical_fiber,
    snapshot_fiber_fields,
)
from fem.beam_column_2d import create_uniform_vertical_tower_mesh
from fem.tower_loading import AsymmetricCyclicTopForceHistory
from fem.tower_system_2d import top_horizontal_load_vector
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
    ViscoplasticTowerResponse,
    solve_nonlinear_tower_load_step,
)
from material.viscoplastic_damage_1d import MaterialParameters


def _snapshot_response(
    response: ViscoplasticTowerResponse,
    strain_snapshots: List[np.ndarray],
    stress_snapshots: List[np.ndarray],
    state_snapshots: List[np.ndarray],
) -> None:
    """Append one complete fiber-field snapshot."""
    fiber_strains, fiber_stresses, fiber_states = (
        snapshot_fiber_fields(response)
    )
    strain_snapshots.append(fiber_strains)
    stress_snapshots.append(fiber_stresses)
    state_snapshots.append(fiber_states)


def run_nonlinear_asymmetric_analysis(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    loading: AsymmetricCyclicTopForceHistory,
    max_iterations: int = 40,
) -> NonlinearCyclicResponse:
    """
    Run explicit preload followed by complete asymmetric periodic cycles.

    The preload duration is T/4 and contains increments_per_cycle/4 equal
    increments, so the preload and periodic stages use the same time step.
    """
    if not isinstance(
        loading,
        AsymmetricCyclicTopForceHistory,
    ):
        raise TypeError(
            "loading must be an AsymmetricCyclicTopForceHistory."
        )

    mesh = create_uniform_vertical_tower_mesh(
        height=configuration.height,
        n_elements=configuration.n_elements,
    )
    geometry = create_tower_geometry(configuration)
    system = ViscoplasticDamageTowerSystem2D(
        mesh=mesh,
        tower_geometry=geometry,
        material=material,
        n_gauss=configuration.n_gauss,
        n_circumferential=configuration.n_circumferential,
        n_radial=configuration.n_radial,
    )

    preload_increments = loading.increments_per_cycle // 4
    preload_duration = (
        preload_increments * loading.time_increment
    )

    preload_solution = None
    for preload_index in range(1, preload_increments + 1):
        fraction = (
            float(preload_index)
            / float(preload_increments)
        )
        preload_force = fraction * loading.mean_force
        load_vector = top_horizontal_load_vector(
            mesh=mesh,
            horizontal_force=float(preload_force),
        )
        preload_solution = solve_nonlinear_tower_load_step(
            system=system,
            time=float(
                preload_index * loading.time_increment
            ),
            load_vector=load_vector,
            max_iterations=max_iterations,
            relative_residual_tolerance=1.0e-8,
            absolute_residual_tolerance=1.0e-4,
        )

    if preload_solution is None:
        raise RuntimeError(
            "The asymmetric preload produced no solution."
        )

    n_points = loading.n_time_points
    analysis_times = preload_duration + loading.times

    top_displacements = np.empty(n_points, dtype=np.float64)
    top_rotations = np.empty(n_points, dtype=np.float64)
    base_horizontal_reactions = np.empty(
        n_points,
        dtype=np.float64,
    )
    base_moment_reactions = np.empty(
        n_points,
        dtype=np.float64,
    )
    iterations = np.empty(n_points, dtype=np.int64)
    residual_norms = np.empty(n_points, dtype=np.float64)

    strain_snapshots = []
    stress_snapshots = []
    state_snapshots = []

    def record_solution(index, solution):
        top_displacements[index] = solution.displacements[-3]
        top_rotations[index] = solution.displacements[-1]
        base_horizontal_reactions[index] = solution.reactions[0]
        base_moment_reactions[index] = solution.reactions[2]
        iterations[index] = solution.iterations
        residual_norms[index] = solution.residual_norm
        _snapshot_response(
            response=solution.response,
            strain_snapshots=strain_snapshots,
            stress_snapshots=stress_snapshots,
            state_snapshots=state_snapshots,
        )

    record_solution(0, preload_solution)

    for time_index in range(1, n_points):
        load_vector = top_horizontal_load_vector(
            mesh=mesh,
            horizontal_force=float(
                loading.forces[time_index]
            ),
        )
        solution = solve_nonlinear_tower_load_step(
            system=system,
            time=float(analysis_times[time_index]),
            load_vector=load_vector,
            max_iterations=max_iterations,
            relative_residual_tolerance=1.0e-8,
            absolute_residual_tolerance=1.0e-4,
        )
        record_solution(time_index, solution)

    fiber_strain_history = np.stack(
        strain_snapshots,
        axis=0,
    )
    fiber_stress_history = np.stack(
        stress_snapshots,
        axis=0,
    )
    fiber_state_history = np.stack(
        state_snapshots,
        axis=0,
    )

    critical_location = locate_critical_fiber(
        fiber_states=fiber_state_history,
    )
    element_index, gauss_index, fiber_index = critical_location
    critical_element = system.elements[element_index]
    critical_section = critical_element.sections[gauss_index]

    return NonlinearCyclicResponse(
        loading=loading,
        material=material,
        analysis_times=analysis_times,
        top_displacements=top_displacements,
        top_rotations=top_rotations,
        base_horizontal_reactions=base_horizontal_reactions,
        base_moment_reactions=base_moment_reactions,
        iterations=iterations,
        residual_norms=residual_norms,
        fiber_strains=fiber_strain_history,
        fiber_stresses=fiber_stress_history,
        fiber_states=fiber_state_history,
        critical_location=critical_location,
        critical_height=float(
            critical_element.gauss_heights[gauss_index]
        ),
        critical_y_coordinate=float(
            critical_section.section.y_coordinates[
                fiber_index
            ]
        ),
    )

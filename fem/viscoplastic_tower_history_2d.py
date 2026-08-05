# -*- coding: utf-8 -*-
"""
Sequential quasi-static solution of a pulsating offshore-wind-turbine
tower load history with the viscoplastic-damage beam system.

The positive-mean loading history begins at a nonzero force. Therefore,
the first stored loading point is reached through an explicit preload
step. Subsequent history points are advanced with the loading time
increment.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from fem.tower_loading import PulsatingTopForceHistory
from fem.tower_system_2d import (
    cantilever_base_fixed_dofs,
    top_horizontal_load_vector,
)
from fem.viscoplastic_tower_system_2d import (
    ViscoplasticDamageTowerSystem2D,
    solve_nonlinear_tower_load_step,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _finite(value: float, name: str) -> float:
    """Return a finite scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error
    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _positive(value: float, name: str) -> float:
    """Return a finite positive scalar."""
    result = _finite(value, name)
    if result <= 0.0:
        raise ValueError(name + " must be positive.")
    return result


def _nonnegative(value: float, name: str) -> float:
    """Return a finite non-negative scalar."""
    result = _finite(value, name)
    if result < 0.0:
        raise ValueError(name + " must be non-negative.")
    return result


def _integer(value: int, name: str) -> int:
    """Return a validated integer."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


def _fixed_dofs(
    values: Sequence[int],
    n_dof: int,
) -> IntArray:
    """Return validated fixed global DOF indices."""
    array = np.asarray(values, dtype=np.int64)
    if array.ndim != 1:
        raise ValueError(
            "fixed_dofs must be one-dimensional."
        )
    if np.any(array < 0) or np.any(array >= n_dof):
        raise ValueError(
            "fixed_dofs contain an out-of-range index."
        )
    if np.unique(array).size != array.size:
        raise ValueError(
            "fixed_dofs must not contain duplicates."
        )
    return np.sort(array)


@dataclass(frozen=True)
class NonlinearPulsatingTowerHistory:
    """Selected response histories from sequential nonlinear steps."""

    loading: PulsatingTopForceHistory
    preload_duration: float
    analysis_times: FloatArray
    displacements: FloatArray
    reactions: FloatArray
    iterations: IntArray
    residual_norms: FloatArray
    displacement_increment_norms: FloatArray
    maximum_damages: FloatArray
    fixed_dofs: IntArray

    def __post_init__(self) -> None:
        if not isinstance(
            self.loading,
            PulsatingTopForceHistory,
        ):
            raise TypeError(
                "loading must be a PulsatingTopForceHistory."
            )

        preload_duration = _positive(
            self.preload_duration,
            "preload_duration",
        )
        analysis_times = np.asarray(
            self.analysis_times,
            dtype=np.float64,
        )
        displacements = np.asarray(
            self.displacements,
            dtype=np.float64,
        )
        reactions = np.asarray(
            self.reactions,
            dtype=np.float64,
        )
        iterations = np.asarray(
            self.iterations,
            dtype=np.int64,
        )
        residual_norms = np.asarray(
            self.residual_norms,
            dtype=np.float64,
        )
        displacement_increment_norms = np.asarray(
            self.displacement_increment_norms,
            dtype=np.float64,
        )
        maximum_damages = np.asarray(
            self.maximum_damages,
            dtype=np.float64,
        )

        n_points = self.loading.n_time_points
        if analysis_times.shape != (n_points,):
            raise ValueError(
                "analysis_times must match the loading history."
            )
        if displacements.ndim != 2:
            raise ValueError(
                "displacements must have shape "
                "(n_time_points, n_dof)."
            )
        if displacements.shape[0] != n_points:
            raise ValueError(
                "displacements must match the loading history."
            )
        if reactions.shape != displacements.shape:
            raise ValueError(
                "reactions must match displacements."
            )
        if iterations.shape != (n_points,):
            raise ValueError(
                "iterations must match the loading history."
            )
        if residual_norms.shape != (n_points,):
            raise ValueError(
                "residual_norms must match the loading history."
            )
        if displacement_increment_norms.shape != (n_points,):
            raise ValueError(
                "displacement_increment_norms must match "
                "the loading history."
            )
        if maximum_damages.shape != (n_points,):
            raise ValueError(
                "maximum_damages must match the loading history."
            )

        arrays = (
            analysis_times,
            displacements,
            reactions,
            residual_norms,
            displacement_increment_norms,
            maximum_damages,
        )
        if any(np.any(~np.isfinite(array)) for array in arrays):
            raise ValueError(
                "All stored response histories must be finite."
            )
        if np.any(np.diff(analysis_times) <= 0.0):
            raise ValueError(
                "analysis_times must be strictly increasing."
            )
        if np.any(iterations < 1):
            raise ValueError(
                "Every stored iteration count must be positive."
            )
        if np.any(residual_norms < 0.0):
            raise ValueError(
                "residual_norms must be non-negative."
            )
        if np.any(displacement_increment_norms < 0.0):
            raise ValueError(
                "displacement_increment_norms must be "
                "non-negative."
            )
        if np.any(maximum_damages < 0.0):
            raise ValueError(
                "maximum_damages must be non-negative."
            )

        fixed_dofs = _fixed_dofs(
            self.fixed_dofs,
            displacements.shape[1],
        )

        object.__setattr__(
            self,
            "preload_duration",
            preload_duration,
        )
        object.__setattr__(
            self,
            "analysis_times",
            analysis_times.copy(),
        )
        object.__setattr__(
            self,
            "displacements",
            displacements.copy(),
        )
        object.__setattr__(
            self,
            "reactions",
            reactions.copy(),
        )
        object.__setattr__(
            self,
            "iterations",
            iterations.copy(),
        )
        object.__setattr__(
            self,
            "residual_norms",
            residual_norms.copy(),
        )
        object.__setattr__(
            self,
            "displacement_increment_norms",
            displacement_increment_norms.copy(),
        )
        object.__setattr__(
            self,
            "maximum_damages",
            maximum_damages.copy(),
        )
        object.__setattr__(
            self,
            "fixed_dofs",
            fixed_dofs,
        )

    @property
    def physical_times(self) -> FloatArray:
        """Return the original loading-time coordinates."""
        return self.loading.times.copy()

    @property
    def top_horizontal_displacements(self) -> FloatArray:
        """Return horizontal displacement of the tower-top node."""
        return self.displacements[:, -3].copy()

    @property
    def top_rotations(self) -> FloatArray:
        """Return rotation of the tower-top node."""
        return self.displacements[:, -1].copy()

    @property
    def base_horizontal_reactions(self) -> FloatArray:
        """Return horizontal reaction at the fixed tower base."""
        return self.reactions[:, 0].copy()

    @property
    def base_bending_reactions(self) -> FloatArray:
        """Return bending reaction at the fixed tower base."""
        return self.reactions[:, 2].copy()

    @property
    def final_maximum_damage(self) -> float:
        """Return maximum tower damage at the final history point."""
        return float(self.maximum_damages[-1])


def solve_pulsating_tower_history(
    system: ViscoplasticDamageTowerSystem2D,
    loading: PulsatingTopForceHistory,
    fixed_dofs: Optional[Sequence[int]] = None,
    preload_duration: Optional[float] = None,
    max_iterations: int = 30,
    relative_residual_tolerance: float = 1.0e-8,
    absolute_residual_tolerance: float = 1.0e-6,
    strain_perturbation: float = 1.0e-8,
) -> NonlinearPulsatingTowerHistory:
    """
    Advance every point of a pulsating tower-top force history.

    Since ``loading.forces[0]`` is generally nonzero, the first point
    is applied through a preload step. By default, its duration equals
    one regular loading increment. The analysis times are therefore

        committed_time + preload_duration + loading.times.

    The stored physical times remain exactly ``loading.times``.
    """
    if not isinstance(
        system,
        ViscoplasticDamageTowerSystem2D,
    ):
        raise TypeError(
            "system must be a "
            "ViscoplasticDamageTowerSystem2D."
        )
    if not isinstance(
        loading,
        PulsatingTopForceHistory,
    ):
        raise TypeError(
            "loading must be a PulsatingTopForceHistory."
        )

    if fixed_dofs is None:
        fixed = cantilever_base_fixed_dofs(system.mesh)
    else:
        fixed = _fixed_dofs(
            fixed_dofs,
            system.mesh.n_dof,
        )

    if preload_duration is None:
        preload = loading.time_increment
    else:
        preload = _positive(
            preload_duration,
            "preload_duration",
        )

    max_iterations = _integer(
        max_iterations,
        "max_iterations",
    )
    if max_iterations < 1:
        raise ValueError(
            "max_iterations must be at least 1."
        )
    relative_residual_tolerance = _nonnegative(
        relative_residual_tolerance,
        "relative_residual_tolerance",
    )
    absolute_residual_tolerance = _nonnegative(
        absolute_residual_tolerance,
        "absolute_residual_tolerance",
    )
    strain_perturbation = _positive(
        strain_perturbation,
        "strain_perturbation",
    )

    n_points = loading.n_time_points
    n_dof = system.mesh.n_dof
    start_time = system.committed_time
    analysis_times = (
        start_time + preload + loading.times
    )

    displacements = np.empty(
        (n_points, n_dof),
        dtype=np.float64,
    )
    reactions = np.empty(
        (n_points, n_dof),
        dtype=np.float64,
    )
    iterations = np.empty(
        n_points,
        dtype=np.int64,
    )
    residual_norms = np.empty(
        n_points,
        dtype=np.float64,
    )
    displacement_increment_norms = np.empty(
        n_points,
        dtype=np.float64,
    )
    maximum_damages = np.empty(
        n_points,
        dtype=np.float64,
    )

    for time_index in range(n_points):
        load_vector = top_horizontal_load_vector(
            mesh=system.mesh,
            horizontal_force=float(
                loading.forces[time_index]
            ),
        )
        step = solve_nonlinear_tower_load_step(
            system=system,
            time=float(analysis_times[time_index]),
            load_vector=load_vector,
            fixed_dofs=fixed,
            max_iterations=max_iterations,
            relative_residual_tolerance=(
                relative_residual_tolerance
            ),
            absolute_residual_tolerance=(
                absolute_residual_tolerance
            ),
            strain_perturbation=strain_perturbation,
        )

        displacements[time_index] = step.displacements
        reactions[time_index] = step.reactions
        iterations[time_index] = step.iterations
        residual_norms[time_index] = step.residual_norm
        displacement_increment_norms[time_index] = (
            step.displacement_increment_norm
        )
        maximum_damages[time_index] = (
            step.response.maximum_damage
        )

    return NonlinearPulsatingTowerHistory(
        loading=loading,
        preload_duration=preload,
        analysis_times=analysis_times,
        displacements=displacements,
        reactions=reactions,
        iterations=iterations,
        residual_norms=residual_norms,
        displacement_increment_norms=(
            displacement_increment_norms
        ),
        maximum_damages=maximum_damages,
        fixed_dofs=fixed,
    )

# -*- coding: utf-8 -*-
"""
One-dimensional two-node bar finite elements coupled with the
cyclic viscoplastic-damage material model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from material.viscoplastic_damage_1d import (
    MaterialParameters,
    evaluate_state,
    rk4_step,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class BarMesh1D:
    """Uniform or non-uniform one-dimensional bar mesh."""

    coordinates: FloatArray
    connectivity: IntArray

    def __post_init__(self) -> None:
        coordinates = np.asarray(self.coordinates, dtype=np.float64)
        connectivity = np.asarray(self.connectivity, dtype=np.int64)

        if coordinates.ndim != 1:
            raise ValueError("coordinates must be a one-dimensional array.")
        if coordinates.size < 2:
            raise ValueError("The mesh must contain at least two nodes.")
        if np.any(np.diff(coordinates) <= 0.0):
            raise ValueError("Node coordinates must be strictly increasing.")

        if connectivity.ndim != 2 or connectivity.shape[1] != 2:
            raise ValueError("connectivity must have shape (n_elements, 2).")
        if connectivity.shape[0] != coordinates.size - 1:
            raise ValueError(
                "A continuous bar mesh must have n_nodes - 1 elements."
            )

        expected = np.column_stack(
            (
                np.arange(coordinates.size - 1, dtype=np.int64),
                np.arange(1, coordinates.size, dtype=np.int64),
            )
        )
        if not np.array_equal(connectivity, expected):
            raise ValueError(
                "Elements must connect consecutive nodes in ascending order."
            )

        object.__setattr__(self, "coordinates", coordinates)
        object.__setattr__(self, "connectivity", connectivity)

    @property
    def n_nodes(self) -> int:
        return int(self.coordinates.size)

    @property
    def n_elements(self) -> int:
        return int(self.connectivity.shape[0])

    @property
    def element_lengths(self) -> FloatArray:
        return np.diff(self.coordinates)


@dataclass
class BarResponse:
    """Time history returned by the displacement-controlled bar analysis."""

    time: FloatArray
    displacement: FloatArray
    strain: FloatArray
    stress: FloatArray
    state: FloatArray
    reaction_left: FloatArray
    reaction_right: FloatArray
    newton_iterations: NDArray[np.int64]


def create_uniform_bar_mesh(
    length: float,
    n_elements: int,
) -> BarMesh1D:
    """Create a uniform two-node bar mesh on [0, length]."""
    if length <= 0.0:
        raise ValueError("length must be positive.")
    if n_elements < 1:
        raise ValueError("n_elements must be at least 1.")

    coordinates = np.linspace(
        0.0,
        float(length),
        int(n_elements) + 1,
        dtype=np.float64,
    )
    connectivity = np.column_stack(
        (
            np.arange(n_elements, dtype=np.int64),
            np.arange(1, n_elements + 1, dtype=np.int64),
        )
    )
    return BarMesh1D(coordinates, connectivity)


def _integrate_material_trial(
    time_old: float,
    time_step: float,
    strain_old: float,
    strain_new: float,
    state_old: FloatArray,
    material: MaterialParameters,
) -> Tuple[FloatArray, float]:
    """
    Integrate one material point from the committed old state to a trial
    end-of-step strain using a linear strain path inside the time step.
    """

    def strain_function(current_time: float) -> float:
        theta = (current_time - time_old) / time_step
        theta = float(np.clip(theta, 0.0, 1.0))
        return float(strain_old + theta * (strain_new - strain_old))

    state_new = rk4_step(
        time=time_old,
        state=np.asarray(state_old, dtype=np.float64).copy(),
        time_step=time_step,
        strain_function=strain_function,
        material=material,
    )
    stress_new = float(
        evaluate_state(strain_new, state_new, material)[0]
    )
    return state_new, stress_new


def _numerical_material_tangent(
    time_old: float,
    time_step: float,
    strain_old: float,
    strain_new: float,
    state_old: FloatArray,
    material: MaterialParameters,
) -> float:
    """
    Compute the end-of-step algorithmic tangent d(sigma)/d(epsilon)
    by a centered finite difference of the complete material update.
    """
    perturbation = max(
        1.0e-9,
        1.0e-6 * max(abs(strain_new), 1.0e-3),
    )

    _, stress_plus = _integrate_material_trial(
        time_old,
        time_step,
        strain_old,
        strain_new + perturbation,
        state_old,
        material,
    )
    _, stress_minus = _integrate_material_trial(
        time_old,
        time_step,
        strain_old,
        strain_new - perturbation,
        state_old,
        material,
    )

    tangent = (stress_plus - stress_minus) / (2.0 * perturbation)

    if not np.isfinite(tangent):
        raise FloatingPointError(
            "Non-finite material tangent encountered."
        )
    if abs(tangent) <= np.finfo(float).eps:
        raise np.linalg.LinAlgError(
            "The material tangent is numerically zero."
        )

    return float(tangent)


def _solve_tridiagonal(
    lower: FloatArray,
    diagonal: FloatArray,
    upper: FloatArray,
    right_hand_side: FloatArray,
) -> FloatArray:
    """Solve a tridiagonal linear system with the Thomas algorithm."""
    diagonal_work = np.asarray(
        diagonal,
        dtype=np.float64,
    ).copy()
    rhs_work = np.asarray(
        right_hand_side,
        dtype=np.float64,
    ).copy()
    lower = np.asarray(lower, dtype=np.float64)
    upper = np.asarray(upper, dtype=np.float64)

    n = diagonal_work.size
    if n == 0:
        return np.empty(0, dtype=np.float64)
    if lower.size != max(n - 1, 0) or upper.size != max(n - 1, 0):
        raise ValueError("Invalid tridiagonal system dimensions.")

    for row in range(1, n):
        pivot = diagonal_work[row - 1]
        if abs(pivot) <= np.finfo(float).eps:
            raise np.linalg.LinAlgError(
                "Zero pivot in tridiagonal solver."
            )
        multiplier = lower[row - 1] / pivot
        diagonal_work[row] -= multiplier * upper[row - 1]
        rhs_work[row] -= multiplier * rhs_work[row - 1]

    solution = np.empty(n, dtype=np.float64)
    if abs(diagonal_work[-1]) <= np.finfo(float).eps:
        raise np.linalg.LinAlgError(
            "Zero pivot in tridiagonal solver."
        )
    solution[-1] = rhs_work[-1] / diagonal_work[-1]

    for row in range(n - 2, -1, -1):
        if abs(diagonal_work[row]) <= np.finfo(float).eps:
            raise np.linalg.LinAlgError(
                "Zero pivot in tridiagonal solver."
            )
        solution[row] = (
            rhs_work[row] - upper[row] * solution[row + 1]
        ) / diagonal_work[row]

    return solution


def _assemble_trial_system(
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    time_old: float,
    time_step: float,
    displacement_trial: FloatArray,
    strain_old: FloatArray,
    state_old: FloatArray,
) -> Tuple[
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
    FloatArray,
]:
    """Assemble internal force, tangent and trial material quantities."""
    n_nodes = mesh.n_nodes
    n_elements = mesh.n_elements

    internal_force = np.zeros(n_nodes, dtype=np.float64)
    diagonal = np.zeros(n_nodes, dtype=np.float64)
    off_diagonal = np.zeros(n_nodes - 1, dtype=np.float64)

    strain_trial = np.zeros(n_elements, dtype=np.float64)
    stress_trial = np.zeros(n_elements, dtype=np.float64)
    state_trial = np.zeros((n_elements, 4), dtype=np.float64)

    for element in range(n_elements):
        node_left, node_right = mesh.connectivity[element]
        element_length = (
            mesh.coordinates[node_right] - mesh.coordinates[node_left]
        )

        strain = (
            displacement_trial[node_right]
            - displacement_trial[node_left]
        ) / element_length

        state_new, stress = _integrate_material_trial(
            time_old=time_old,
            time_step=time_step,
            strain_old=float(strain_old[element]),
            strain_new=float(strain),
            state_old=state_old[element],
            material=materials[element],
        )
        tangent = _numerical_material_tangent(
            time_old=time_old,
            time_step=time_step,
            strain_old=float(strain_old[element]),
            strain_new=float(strain),
            state_old=state_old[element],
            material=materials[element],
        )

        axial_force = area * stress
        element_stiffness = area * tangent / element_length

        internal_force[node_left] -= axial_force
        internal_force[node_right] += axial_force

        diagonal[node_left] += element_stiffness
        diagonal[node_right] += element_stiffness
        off_diagonal[node_left] -= element_stiffness

        strain_trial[element] = strain
        stress_trial[element] = stress
        state_trial[element] = state_new

    return (
        internal_force,
        diagonal,
        off_diagonal,
        strain_trial,
        stress_trial,
        state_trial,
    )


def solve_displacement_controlled_bar(
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    time: FloatArray,
    right_displacement: FloatArray,
    *,
    newton_tolerance: float = 1.0e-8,
    newton_absolute_tolerance: float = 1.0e-6,
    maximum_newton_iterations: int = 30,
) -> BarResponse:
    """
    Solve a quasi-static one-dimensional bar.

    Boundary conditions:
        u(0, t) = 0
        u(L, t) = right_displacement(t)

    There are no distributed loads or nodal forces. Interior equilibrium is
    enforced by Newton iterations at every time step.
    """
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if newton_tolerance <= 0.0:
        raise ValueError("newton_tolerance must be positive.")
    if newton_absolute_tolerance <= 0.0:
        raise ValueError(
            "newton_absolute_tolerance must be positive."
        )
    if maximum_newton_iterations < 1:
        raise ValueError(
            "maximum_newton_iterations must be at least 1."
        )

    time = np.asarray(time, dtype=np.float64)
    right_displacement = np.asarray(
        right_displacement,
        dtype=np.float64,
    )

    if time.ndim != 1 or right_displacement.ndim != 1:
        raise ValueError(
            "time and right_displacement must be one-dimensional arrays."
        )
    if time.size < 2:
        raise ValueError("At least two time points are required.")
    if right_displacement.shape != time.shape:
        raise ValueError(
            "right_displacement must have the same shape as time."
        )
    if np.any(np.diff(time) <= 0.0):
        raise ValueError("time must be strictly increasing.")
    if abs(right_displacement[0]) > 1.0e-14:
        raise ValueError(
            "The initial right-end displacement must be zero."
        )

    n_steps = time.size
    n_nodes = mesh.n_nodes
    n_elements = mesh.n_elements

    displacement_history = np.zeros(
        (n_steps, n_nodes),
        dtype=np.float64,
    )
    strain_history = np.zeros(
        (n_steps, n_elements),
        dtype=np.float64,
    )
    stress_history = np.zeros(
        (n_steps, n_elements),
        dtype=np.float64,
    )
    state_history = np.zeros(
        (n_steps, n_elements, 4),
        dtype=np.float64,
    )
    reaction_left = np.zeros(n_steps, dtype=np.float64)
    reaction_right = np.zeros(n_steps, dtype=np.float64)
    newton_iterations = np.zeros(n_steps, dtype=np.int64)

    displacement = np.zeros(n_nodes, dtype=np.float64)
    strain_committed = np.zeros(n_elements, dtype=np.float64)
    state_committed = np.zeros((n_elements, 4), dtype=np.float64)

    normalized_coordinate = (
        mesh.coordinates - mesh.coordinates[0]
    ) / (
        mesh.coordinates[-1] - mesh.coordinates[0]
    )

    for step in range(1, n_steps):
        time_old = float(time[step - 1])
        time_step = float(time[step] - time[step - 1])

        displacement_trial = displacement.copy()
        displacement_increment = (
            right_displacement[step]
            - right_displacement[step - 1]
        )
        displacement_trial += (
            normalized_coordinate * displacement_increment
        )
        displacement_trial[0] = 0.0
        displacement_trial[-1] = right_displacement[step]

        # Initialize trial quantities before the Newton loop so that their
        # existence is explicit to both Python readers and static analyzers.
        internal_force = np.zeros(n_nodes, dtype=np.float64)
        strain_trial = strain_committed.copy()
        stress_trial = np.zeros(n_elements, dtype=np.float64)
        state_trial = state_committed.copy()

        converged = False

        for iteration in range(1, maximum_newton_iterations + 1):
            (
                internal_force,
                diagonal,
                off_diagonal,
                strain_trial,
                stress_trial,
                state_trial,
            ) = _assemble_trial_system(
                mesh=mesh,
                area=area,
                materials=materials,
                time_old=time_old,
                time_step=time_step,
                displacement_trial=displacement_trial,
                strain_old=strain_committed,
                state_old=state_committed,
            )

            residual = internal_force[1:-1]
            residual_norm = (
                float(np.linalg.norm(residual, ord=np.inf))
                if residual.size > 0
                else 0.0
            )
            force_scale = max(
                float(np.linalg.norm(internal_force, ord=np.inf)),
                1.0,
            )
            convergence_limit = (
                newton_absolute_tolerance
                + newton_tolerance * force_scale
            )

            if residual_norm <= convergence_limit:
                converged = True
                newton_iterations[step] = iteration
                break

            diagonal_free = diagonal[1:-1]
            off_diagonal_free = off_diagonal[1:-1]

            displacement_correction = _solve_tridiagonal(
                lower=off_diagonal_free.copy(),
                diagonal=diagonal_free,
                upper=off_diagonal_free.copy(),
                right_hand_side=-residual,
            )
            displacement_trial[1:-1] += displacement_correction

        if not converged:
            raise RuntimeError(
                "Newton iteration did not converge at "
                f"time step {step}, t={time[step]:.6g}."
            )

        displacement = displacement_trial
        strain_committed = strain_trial
        state_committed = state_trial

        displacement_history[step] = displacement
        strain_history[step] = strain_committed
        stress_history[step] = stress_trial
        state_history[step] = state_committed
        reaction_left[step] = internal_force[0]
        reaction_right[step] = internal_force[-1]

    return BarResponse(
        time=time.copy(),
        displacement=displacement_history,
        strain=strain_history,
        stress=stress_history,
        state=state_history,
        reaction_left=reaction_left,
        reaction_right=reaction_right,
        newton_iterations=newton_iterations,
    )

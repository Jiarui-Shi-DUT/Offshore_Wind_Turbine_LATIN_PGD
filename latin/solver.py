# -*- coding: utf-8 -*-
"""
Full-order outer iteration driver for the one-dimensional LATIN solver.

This module assembles the previously validated building blocks into the
alternating LATIN sequence

    global state s_i
        -> nonlinear local stage s_hat_(i+1/2)
        -> descent search directions
        -> full-order global candidate s_breve_(i+1)
        -> relaxed global state s_(i+1)
        -> relative LATIN indicator xi_(i+1).

The driver is deliberately full-order. It provides the reference iterative
solution against which the later PGD-reduced global stage will be compared.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np
from numpy.typing import NDArray

from fem.bar_1d import BarMesh1D
from latin.global_stage import solve_full_order_global_stage
from latin.iteration_control import (
    relative_latin_indicator,
    relax_global_state,
)
from latin.local_stage import solve_local_stage
from latin.search_directions import (
    DescentSearchDirections,
    compute_descent_search_directions,
)
from latin.state import LatinState
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass(frozen=True)
class FullOrderLatinResult:
    """Result and convergence history of the full-order LATIN iteration."""

    state: LatinState
    local_state: LatinState
    directions: DescentSearchDirections
    indicator_history: FloatArray
    converged: bool
    iterations: int
    tolerance: float
    relaxation: float

    @property
    def final_indicator(self) -> float:
        """Return the final relative LATIN indicator."""
        return float(self.indicator_history[-1])


def _validate_solver_inputs(
    initial_state: LatinState,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    tolerance: float,
    max_iterations: int,
    relaxation: float,
) -> None:
    """Validate the common discretisation and iteration parameters."""
    if area <= 0.0:
        raise ValueError("area must be positive.")
    if initial_state.n_elements != mesh.n_elements:
        raise ValueError(
            "initial_state and mesh must contain the same "
            "number of elements."
        )
    if len(materials) != mesh.n_elements:
        raise ValueError(
            "One MaterialParameters object is required per element."
        )
    if tolerance <= 0.0 or not np.isfinite(tolerance):
        raise ValueError("tolerance must be positive and finite.")
    if max_iterations < 1:
        raise ValueError("max_iterations must be at least one.")
    if not 0.0 < relaxation <= 1.0:
        raise ValueError(
            "relaxation must satisfy 0 < relaxation <= 1."
        )


def solve_full_order_latin(
    initial_state: LatinState,
    mesh: BarMesh1D,
    area: float,
    materials: Sequence[MaterialParameters],
    *,
    tolerance: float = 1.0e-3,
    max_iterations: int = 50,
    relaxation: float = 0.8,
) -> FullOrderLatinResult:
    """
    Run the full-order LATIN iterations over the complete time-space domain.

    Parameters
    ----------
    initial_state:
        Globally admissible initial state, normally the elastic solution.
    mesh:
        One-dimensional finite-element mesh.
    area:
        Constant bar cross-sectional area.
    materials:
        One material-parameter object for each element.
    tolerance:
        Stopping tolerance for the relative LATIN indicator.
    max_iterations:
        Maximum number of complete local-global iterations.
    relaxation:
        Relaxation parameter applied after every global stage. The reference
        paper uses 0.8.

    Returns
    -------
    FullOrderLatinResult
        Final relaxed global state, last local state, last search directions,
        convergence history and termination information.
    """
    _validate_solver_inputs(
        initial_state=initial_state,
        mesh=mesh,
        area=area,
        materials=materials,
        tolerance=tolerance,
        max_iterations=max_iterations,
        relaxation=relaxation,
    )

    global_state = initial_state.copy()
    indicator_values = []

    last_local_state = initial_state.copy()
    last_directions: Optional[DescentSearchDirections] = None

    for iteration in range(1, max_iterations + 1):
        local_state = solve_local_stage(
            global_state=global_state,
            materials=materials,
        )
        directions = compute_descent_search_directions(
            local_state=local_state,
            materials=materials,
        )
        candidate_state = solve_full_order_global_stage(
            global_state=global_state,
            local_state=local_state,
            directions=directions,
            mesh=mesh,
            area=area,
            materials=materials,
        ).state
        next_global_state = relax_global_state(
            previous_state=global_state,
            candidate_state=candidate_state,
            relaxation=relaxation,
        )
        indicator = relative_latin_indicator(
            local_state=local_state,
            global_state=next_global_state,
            directions=directions,
            mesh=mesh,
            area=area,
            materials=materials,
        )

        if not np.isfinite(indicator):
            raise FloatingPointError(
                "The LATIN indicator became non-finite at iteration "
                f"{iteration}."
            )

        indicator_values.append(indicator)
        global_state = next_global_state
        last_local_state = local_state
        last_directions = directions

        if indicator <= tolerance:
            return FullOrderLatinResult(
                state=global_state,
                local_state=last_local_state,
                directions=last_directions,
                indicator_history=np.asarray(
                    indicator_values,
                    dtype=np.float64,
                ),
                converged=True,
                iterations=iteration,
                tolerance=float(tolerance),
                relaxation=float(relaxation),
            )

    if last_directions is None:
        raise RuntimeError(
            "The LATIN iteration terminated before its first iteration."
        )

    return FullOrderLatinResult(
        state=global_state,
        local_state=last_local_state,
        directions=last_directions,
        indicator_history=np.asarray(
            indicator_values,
            dtype=np.float64,
        ),
        converged=False,
        iterations=max_iterations,
        tolerance=float(tolerance),
        relaxation=float(relaxation),
    )

# -*- coding: utf-8 -*-
"""
Stateful viscoplastic-damage response of an annular fiber section.

Each fiber owns an independent material state

    [plastic_strain, alpha, r_bar, damage].

Section conventions
-------------------
    epsilon_f = epsilon_0 - kappa * y_f
    N = sum(sigma_f * A_f)
    M = -sum(sigma_f * y_f * A_f)

The existing material model uses MPa by default, whereas the tower geometry
uses metres. Consequently, ``stress_to_force_factor`` defaults to 1.0e6 so
that MPa times m^2 gives N and MPa times m^3 gives N m.

Trial states are always integrated from the last committed state. Repeated
trial calls therefore do not accumulate plasticity or damage. State variables
become permanent only after ``commit_state()``.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.fiber_section import AnnularFiberSection, compute_fiber_strains
from material.viscoplastic_damage_1d import (
    MaterialParameters,
    MaterialState,
    rk4_step,
    stress_from_state,
)


FloatArray = NDArray[np.float64]


def _finite_scalar(value: float, name: str) -> float:
    """Return a validated finite scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error
    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _positive_scalar(value: float, name: str) -> float:
    """Return a validated positive scalar."""
    result = _finite_scalar(value, name)
    if result <= 0.0:
        raise ValueError(name + " must be positive.")
    return result


def _state_array(
    values: FloatArray,
    n_fibers: int,
    name: str,
) -> FloatArray:
    """Validate and copy an array of fiber states."""
    result = np.asarray(values, dtype=np.float64)
    if result.shape != (n_fibers, 4):
        raise ValueError(
            name + " must have shape (section.n_fibers, 4)."
        )
    if np.any(~np.isfinite(result)):
        raise ValueError(name + " must contain finite values.")
    if np.any(result[:, 3] < 0.0) or np.any(result[:, 3] >= 1.0):
        raise ValueError("Fiber damage must satisfy 0 <= D < 1.")
    return result.copy()


@dataclass(frozen=True)
class ViscoplasticSectionResponse:
    """Trial or committed response of one annular fiber section."""

    time: float
    axial_strain: float
    curvature: float
    fiber_strains: FloatArray
    fiber_stresses: FloatArray
    fiber_states: FloatArray
    axial_force: float
    bending_moment: float
    tangent: Optional[FloatArray] = None

    def __post_init__(self) -> None:
        time = _finite_scalar(self.time, "time")
        axial_strain = _finite_scalar(
            self.axial_strain,
            "axial_strain",
        )
        curvature = _finite_scalar(self.curvature, "curvature")
        axial_force = _finite_scalar(self.axial_force, "axial_force")
        bending_moment = _finite_scalar(
            self.bending_moment,
            "bending_moment",
        )

        strains = np.asarray(self.fiber_strains, dtype=np.float64)
        stresses = np.asarray(self.fiber_stresses, dtype=np.float64)
        states = np.asarray(self.fiber_states, dtype=np.float64)

        if strains.ndim != 1:
            raise ValueError("fiber_strains must be one-dimensional.")
        if stresses.shape != strains.shape:
            raise ValueError("fiber_stresses must match fiber_strains.")
        if states.shape != (strains.size, 4):
            raise ValueError(
                "fiber_states must have shape (n_fibers, 4)."
            )
        if np.any(~np.isfinite(strains)):
            raise ValueError("fiber_strains must be finite.")
        if np.any(~np.isfinite(stresses)):
            raise ValueError("fiber_stresses must be finite.")
        if np.any(~np.isfinite(states)):
            raise ValueError("fiber_states must be finite.")

        tangent = self.tangent
        if tangent is not None:
            tangent = np.asarray(tangent, dtype=np.float64)
            if tangent.shape != (2, 2):
                raise ValueError("tangent must have shape (2, 2).")
            if np.any(~np.isfinite(tangent)):
                raise ValueError("tangent must be finite.")
            tangent = tangent.copy()

        object.__setattr__(self, "time", time)
        object.__setattr__(self, "axial_strain", axial_strain)
        object.__setattr__(self, "curvature", curvature)
        object.__setattr__(self, "fiber_strains", strains.copy())
        object.__setattr__(self, "fiber_stresses", stresses.copy())
        object.__setattr__(self, "fiber_states", states.copy())
        object.__setattr__(self, "axial_force", axial_force)
        object.__setattr__(self, "bending_moment", bending_moment)
        object.__setattr__(self, "tangent", tangent)

    @property
    def resultants(self) -> FloatArray:
        """Return [N, M]."""
        return np.array(
            [self.axial_force, self.bending_moment],
            dtype=np.float64,
        )

    @property
    def plastic_strains(self) -> FloatArray:
        """Return the plastic strain of every fiber."""
        return self.fiber_states[:, 0].copy()

    @property
    def alphas(self) -> FloatArray:
        """Return the kinematic internal variable of every fiber."""
        return self.fiber_states[:, 1].copy()

    @property
    def r_bars(self) -> FloatArray:
        """Return the isotropic internal variable of every fiber."""
        return self.fiber_states[:, 2].copy()

    @property
    def damages(self) -> FloatArray:
        """Return the damage of every fiber."""
        return self.fiber_states[:, 3].copy()

    @property
    def maximum_damage(self) -> float:
        """Return the largest fiber damage."""
        return float(np.max(self.fiber_states[:, 3]))


class ViscoplasticDamageFiberSection:
    """Annular fiber section with independent 1D material states."""

    def __init__(
        self,
        section: AnnularFiberSection,
        material: MaterialParameters,
        stress_to_force_factor: float = 1.0e6,
        initial_state: Optional[MaterialState] = None,
        initial_states: Optional[FloatArray] = None,
        initial_time: float = 0.0,
        initial_axial_strain: float = 0.0,
        initial_curvature: float = 0.0,
    ) -> None:
        if not isinstance(section, AnnularFiberSection):
            raise TypeError(
                "section must be an AnnularFiberSection."
            )
        if not isinstance(material, MaterialParameters):
            raise TypeError(
                "material must be a MaterialParameters object."
            )
        if initial_state is not None and initial_states is not None:
            raise ValueError(
                "Specify initial_state or initial_states, not both."
            )
        if initial_state is not None and not isinstance(
            initial_state,
            MaterialState,
        ):
            raise TypeError("initial_state must be a MaterialState.")

        self.section = section
        self.material = material
        self.stress_to_force_factor = _positive_scalar(
            stress_to_force_factor,
            "stress_to_force_factor",
        )
        self._initial_time = _finite_scalar(
            initial_time,
            "initial_time",
        )
        self._initial_axial_strain = _finite_scalar(
            initial_axial_strain,
            "initial_axial_strain",
        )
        self._initial_curvature = _finite_scalar(
            initial_curvature,
            "initial_curvature",
        )

        if initial_states is not None:
            states = _state_array(
                initial_states,
                section.n_fibers,
                "initial_states",
            )
        else:
            common_state = (
                MaterialState()
                if initial_state is None
                else initial_state
            )
            states = np.tile(
                common_state.to_array(),
                (section.n_fibers, 1),
            ).astype(np.float64, copy=False)

        self._initial_states = states.copy()
        self._committed_states = states.copy()
        self._trial_states = states.copy()
        self._committed_time = self._initial_time
        self._trial_time = self._initial_time
        self._committed_axial_strain = self._initial_axial_strain
        self._trial_axial_strain = self._initial_axial_strain
        self._committed_curvature = self._initial_curvature
        self._trial_curvature = self._initial_curvature
        self._has_uncommitted_trial = False

        self._committed_response = self._response_at_known_state(
            time=self._committed_time,
            axial_strain=self._committed_axial_strain,
            curvature=self._committed_curvature,
            states=self._committed_states,
            tangent=None,
        )
        self._trial_response = self._committed_response

    @property
    def committed_time(self) -> float:
        """Return the time of the last committed state."""
        return float(self._committed_time)

    @property
    def trial_time(self) -> float:
        """Return the current trial time."""
        return float(self._trial_time)

    @property
    def has_uncommitted_trial(self) -> bool:
        """Return whether a trial state awaits commit or rollback."""
        return bool(self._has_uncommitted_trial)

    @property
    def committed_states(self) -> FloatArray:
        """Return a defensive copy of committed fiber states."""
        return self._committed_states.copy()

    @property
    def trial_states(self) -> FloatArray:
        """Return a defensive copy of trial fiber states."""
        return self._trial_states.copy()

    @property
    def committed_response(self) -> ViscoplasticSectionResponse:
        """Return the last committed section response."""
        return self._committed_response

    @property
    def trial_response(self) -> ViscoplasticSectionResponse:
        """Return the current trial section response."""
        return self._trial_response

    def _integrate_fibers(
        self,
        time: float,
        axial_strain: float,
        curvature: float,
    ) -> Tuple[FloatArray, FloatArray, FloatArray]:
        """Integrate every fiber from the last committed state."""
        trial_time = _finite_scalar(time, "time")
        trial_axial_strain = _finite_scalar(
            axial_strain,
            "axial_strain",
        )
        trial_curvature = _finite_scalar(curvature, "curvature")
        time_step = trial_time - self._committed_time
        if time_step <= 0.0:
            raise ValueError(
                "Trial time must be greater than committed time."
            )

        strain_start = compute_fiber_strains(
            section=self.section,
            axial_strain=self._committed_axial_strain,
            curvature=self._committed_curvature,
        )
        strain_end = compute_fiber_strains(
            section=self.section,
            axial_strain=trial_axial_strain,
            curvature=trial_curvature,
        )

        states = np.zeros_like(self._committed_states)
        stresses = np.zeros(
            self.section.n_fibers,
            dtype=np.float64,
        )
        committed_time = self._committed_time

        for fiber_index in range(self.section.n_fibers):
            start_value = float(strain_start[fiber_index])
            end_value = float(strain_end[fiber_index])

            def strain_function(
                current_time: float,
                start: float = start_value,
                end: float = end_value,
            ) -> float:
                fraction = (
                    (float(current_time) - committed_time)
                    / time_step
                )
                return float(start + fraction * (end - start))

            state = rk4_step(
                time=committed_time,
                state=self._committed_states[fiber_index].copy(),
                time_step=time_step,
                strain_function=strain_function,
                material=self.material,
            )
            states[fiber_index] = state
            stresses[fiber_index] = stress_from_state(
                total_strain=end_value,
                plastic_strain=float(state[0]),
                damage=float(state[3]),
                material=self.material,
            )

        return strain_end, stresses, states

    def _resultants(
        self,
        fiber_stresses: FloatArray,
    ) -> FloatArray:
        """Integrate material stresses into [N, M]."""
        stresses = (
            self.stress_to_force_factor
            * np.asarray(fiber_stresses, dtype=np.float64)
        )
        axial_force = float(np.dot(self.section.areas, stresses))
        bending_moment = float(
            -np.dot(
                self.section.areas * self.section.y_coordinates,
                stresses,
            )
        )
        return np.array(
            [axial_force, bending_moment],
            dtype=np.float64,
        )

    def _integrate_response_data(
        self,
        time: float,
        axial_strain: float,
        curvature: float,
    ) -> Tuple[FloatArray, FloatArray, FloatArray, FloatArray]:
        """Return strains, stresses, states, and resultants."""
        strains, stresses, states = self._integrate_fibers(
            time=time,
            axial_strain=axial_strain,
            curvature=curvature,
        )
        return strains, stresses, states, self._resultants(stresses)

    def _numerical_tangent(
        self,
        time: float,
        axial_strain: float,
        curvature: float,
        strain_perturbation: float,
    ) -> FloatArray:
        """Return a central-difference algorithmic section tangent."""
        delta_strain = _positive_scalar(
            strain_perturbation,
            "strain_perturbation",
        )
        radius = float(
            np.max(np.abs(self.section.y_coordinates))
        )
        if radius <= 0.0:
            raise ValueError("Section has no usable bending radius.")
        delta_curvature = delta_strain / radius

        axial_plus = self._integrate_response_data(
            time,
            axial_strain + delta_strain,
            curvature,
        )[3]
        axial_minus = self._integrate_response_data(
            time,
            axial_strain - delta_strain,
            curvature,
        )[3]
        curvature_plus = self._integrate_response_data(
            time,
            axial_strain,
            curvature + delta_curvature,
        )[3]
        curvature_minus = self._integrate_response_data(
            time,
            axial_strain,
            curvature - delta_curvature,
        )[3]

        return np.column_stack(
            (
                (axial_plus - axial_minus) / (2.0 * delta_strain),
                (curvature_plus - curvature_minus)
                / (2.0 * delta_curvature),
            )
        ).astype(np.float64, copy=False)

    def _response_at_known_state(
        self,
        time: float,
        axial_strain: float,
        curvature: float,
        states: FloatArray,
        tangent: Optional[FloatArray],
    ) -> ViscoplasticSectionResponse:
        """Evaluate a response without advancing the material state."""
        strains = compute_fiber_strains(
            section=self.section,
            axial_strain=axial_strain,
            curvature=curvature,
        )
        stresses = np.zeros(
            self.section.n_fibers,
            dtype=np.float64,
        )
        for fiber_index in range(self.section.n_fibers):
            stresses[fiber_index] = stress_from_state(
                total_strain=float(strains[fiber_index]),
                plastic_strain=float(states[fiber_index, 0]),
                damage=float(states[fiber_index, 3]),
                material=self.material,
            )
        resultants = self._resultants(stresses)
        return ViscoplasticSectionResponse(
            time=time,
            axial_strain=axial_strain,
            curvature=curvature,
            fiber_strains=strains,
            fiber_stresses=stresses,
            fiber_states=states,
            axial_force=float(resultants[0]),
            bending_moment=float(resultants[1]),
            tangent=tangent,
        )

    def set_trial_deformation(
        self,
        time: float,
        axial_strain: float,
        curvature: float,
        compute_tangent: bool = False,
        strain_perturbation: float = 1.0e-8,
    ) -> ViscoplasticSectionResponse:
        """Integrate one trial step from the last committed state."""
        if not isinstance(compute_tangent, (bool, np.bool_)):
            raise TypeError("compute_tangent must be Boolean.")

        trial_time = _finite_scalar(time, "time")
        trial_axial_strain = _finite_scalar(
            axial_strain,
            "axial_strain",
        )
        trial_curvature = _finite_scalar(curvature, "curvature")
        strains, stresses, states, resultants = (
            self._integrate_response_data(
                trial_time,
                trial_axial_strain,
                trial_curvature,
            )
        )

        tangent = None
        if compute_tangent:
            tangent = self._numerical_tangent(
                time=trial_time,
                axial_strain=trial_axial_strain,
                curvature=trial_curvature,
                strain_perturbation=strain_perturbation,
            )

        response = ViscoplasticSectionResponse(
            time=trial_time,
            axial_strain=trial_axial_strain,
            curvature=trial_curvature,
            fiber_strains=strains,
            fiber_stresses=stresses,
            fiber_states=states,
            axial_force=float(resultants[0]),
            bending_moment=float(resultants[1]),
            tangent=tangent,
        )
        self._trial_time = trial_time
        self._trial_axial_strain = trial_axial_strain
        self._trial_curvature = trial_curvature
        self._trial_states = states.copy()
        self._trial_response = response
        self._has_uncommitted_trial = True
        return response

    def commit_state(self) -> ViscoplasticSectionResponse:
        """Make the current trial state permanent."""
        if not self._has_uncommitted_trial:
            raise RuntimeError("No uncommitted trial state is available.")
        self._committed_time = self._trial_time
        self._committed_axial_strain = self._trial_axial_strain
        self._committed_curvature = self._trial_curvature
        self._committed_states = self._trial_states.copy()
        self._committed_response = self._response_at_known_state(
            time=self._committed_time,
            axial_strain=self._committed_axial_strain,
            curvature=self._committed_curvature,
            states=self._committed_states,
            tangent=self._trial_response.tangent,
        )
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False
        return self._committed_response

    def revert_to_last_commit(self) -> ViscoplasticSectionResponse:
        """Discard the current trial state."""
        self._trial_time = self._committed_time
        self._trial_axial_strain = self._committed_axial_strain
        self._trial_curvature = self._committed_curvature
        self._trial_states = self._committed_states.copy()
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False
        return self._trial_response

    def revert_to_start(self) -> ViscoplasticSectionResponse:
        """Restore the initial state and deformation."""
        self._committed_time = self._initial_time
        self._committed_axial_strain = self._initial_axial_strain
        self._committed_curvature = self._initial_curvature
        self._committed_states = self._initial_states.copy()
        self._committed_response = self._response_at_known_state(
            time=self._committed_time,
            axial_strain=self._committed_axial_strain,
            curvature=self._committed_curvature,
            states=self._committed_states,
            tangent=None,
        )
        self._trial_time = self._committed_time
        self._trial_axial_strain = self._committed_axial_strain
        self._trial_curvature = self._committed_curvature
        self._trial_states = self._committed_states.copy()
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False
        return self._committed_response

# -*- coding: utf-8 -*-
"""
Global assembly and Newton solution for a two-dimensional tapered tower
built from viscoplastic-damage fiber beam elements.

Every Newton trial is evaluated from the last committed material state.
Only a converged load step is committed.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.beam_column_2d import (
    BeamMesh2D,
    LinearTaperedTowerGeometry,
)
from fem.tower_system_2d import (
    cantilever_base_fixed_dofs,
    element_dof_indices,
    free_dofs_from_fixed,
)
from fem.viscoplastic_beam_column_2d import (
    ViscoplasticBeamElementResponse,
    ViscoplasticDamageBeamElement2D,
)
from material.viscoplastic_damage_1d import MaterialParameters


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


def _vector(
    values: FloatArray,
    size: int,
    name: str,
) -> FloatArray:
    """Return a finite vector with the requested size."""
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (size,):
        raise ValueError(
            name + " must have shape (" + str(size) + ",)."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain finite values.")
    return array.copy()


def _dof_array(
    values: Sequence[int],
    n_dof: int,
    name: str,
) -> IntArray:
    """Validate a one-dimensional set of DOF indices."""
    array = np.asarray(values, dtype=np.int64)
    if array.ndim != 1:
        raise ValueError(name + " must be one-dimensional.")
    if np.any(array < 0) or np.any(array >= n_dof):
        raise ValueError(name + " contains an out-of-range index.")
    if np.unique(array).size != array.size:
        raise ValueError(name + " must not contain duplicates.")
    return np.sort(array)


class NonlinearTowerConvergenceError(RuntimeError):
    """Raised when one nonlinear tower load step does not converge."""


@dataclass(frozen=True)
class ViscoplasticTowerResponse:
    """Trial or committed global response of the nonlinear tower."""

    time: float
    displacements: FloatArray
    internal_force: FloatArray
    tangent: Optional[FloatArray]
    element_dofs: IntArray
    element_responses: Tuple[
        ViscoplasticBeamElementResponse,
        ...,
    ]

    def __post_init__(self) -> None:
        time = _finite(self.time, "time")
        displacements = np.asarray(
            self.displacements,
            dtype=np.float64,
        )
        internal_force = np.asarray(
            self.internal_force,
            dtype=np.float64,
        )
        element_dofs = np.asarray(
            self.element_dofs,
            dtype=np.int64,
        )
        element_responses = tuple(self.element_responses)

        if displacements.ndim != 1:
            raise ValueError(
                "displacements must be one-dimensional."
            )
        if internal_force.shape != displacements.shape:
            raise ValueError(
                "internal_force must match displacements."
            )
        if np.any(~np.isfinite(displacements)):
            raise ValueError("displacements must be finite.")
        if np.any(~np.isfinite(internal_force)):
            raise ValueError("internal_force must be finite.")

        if element_dofs.shape != (len(element_responses), 6):
            raise ValueError(
                "element_dofs must have shape (n_elements, 6)."
            )
        if np.any(element_dofs < 0):
            raise ValueError("element_dofs must be non-negative.")
        if np.any(element_dofs >= displacements.size):
            raise ValueError(
                "element_dofs contain an out-of-range index."
            )
        for response in element_responses:
            if not isinstance(
                response,
                ViscoplasticBeamElementResponse,
            ):
                raise TypeError(
                    "element_responses must contain "
                    "ViscoplasticBeamElementResponse objects."
                )

        tangent = self.tangent
        if tangent is not None:
            tangent_array = np.asarray(
                tangent,
                dtype=np.float64,
            )
            expected_shape = (
                displacements.size,
                displacements.size,
            )
            if tangent_array.shape != expected_shape:
                raise ValueError(
                    "tangent must have shape "
                    + str(expected_shape)
                    + "."
                )
            if np.any(~np.isfinite(tangent_array)):
                raise ValueError("tangent must be finite.")
            tangent = tangent_array.copy()

        object.__setattr__(self, "time", time)
        object.__setattr__(
            self,
            "displacements",
            displacements.copy(),
        )
        object.__setattr__(
            self,
            "internal_force",
            internal_force.copy(),
        )
        object.__setattr__(self, "tangent", tangent)
        object.__setattr__(
            self,
            "element_dofs",
            element_dofs.copy(),
        )
        object.__setattr__(
            self,
            "element_responses",
            element_responses,
        )

    @property
    def maximum_damage(self) -> float:
        """Return the largest fiber damage in the tower."""
        if not self.element_responses:
            return 0.0
        return float(
            max(
                response.maximum_damage
                for response in self.element_responses
            )
        )


@dataclass(frozen=True)
class NonlinearTowerStepSolution:
    """Converged result of one force-controlled tower load step."""

    time: float
    load_vector: FloatArray
    displacements: FloatArray
    reactions: FloatArray
    residual: FloatArray
    fixed_dofs: IntArray
    free_dofs: IntArray
    iterations: int
    residual_norm: float
    displacement_increment_norm: float
    response: ViscoplasticTowerResponse

    def __post_init__(self) -> None:
        time = _finite(self.time, "time")
        displacements = np.asarray(
            self.displacements,
            dtype=np.float64,
        )
        if displacements.ndim != 1:
            raise ValueError(
                "displacements must be one-dimensional."
            )
        n_dof = displacements.size

        load_vector = _vector(
            self.load_vector,
            n_dof,
            "load_vector",
        )
        reactions = _vector(
            self.reactions,
            n_dof,
            "reactions",
        )
        residual = _vector(
            self.residual,
            n_dof,
            "residual",
        )
        fixed_dofs = _dof_array(
            self.fixed_dofs,
            n_dof,
            "fixed_dofs",
        )
        free_dofs = _dof_array(
            self.free_dofs,
            n_dof,
            "free_dofs",
        )

        if np.intersect1d(fixed_dofs, free_dofs).size != 0:
            raise ValueError(
                "fixed_dofs and free_dofs must not overlap."
            )
        if fixed_dofs.size + free_dofs.size != n_dof:
            raise ValueError(
                "fixed_dofs and free_dofs must partition all DOFs."
            )

        iterations = _integer(self.iterations, "iterations")
        if iterations < 1:
            raise ValueError("iterations must be at least 1.")

        residual_norm = _nonnegative(
            self.residual_norm,
            "residual_norm",
        )
        displacement_increment_norm = _nonnegative(
            self.displacement_increment_norm,
            "displacement_increment_norm",
        )
        if not isinstance(
            self.response,
            ViscoplasticTowerResponse,
        ):
            raise TypeError(
                "response must be a ViscoplasticTowerResponse."
            )

        object.__setattr__(self, "time", time)
        object.__setattr__(self, "load_vector", load_vector)
        object.__setattr__(
            self,
            "displacements",
            displacements.copy(),
        )
        object.__setattr__(self, "reactions", reactions)
        object.__setattr__(self, "residual", residual)
        object.__setattr__(self, "fixed_dofs", fixed_dofs)
        object.__setattr__(self, "free_dofs", free_dofs)
        object.__setattr__(self, "iterations", iterations)
        object.__setattr__(
            self,
            "residual_norm",
            residual_norm,
        )
        object.__setattr__(
            self,
            "displacement_increment_norm",
            displacement_increment_norm,
        )


class ViscoplasticDamageTowerSystem2D:
    """Stateful assembly of nonlinear fiber beam elements."""

    def __init__(
        self,
        mesh: BeamMesh2D,
        tower_geometry: LinearTaperedTowerGeometry,
        material: MaterialParameters,
        n_gauss: int = 4,
        n_circumferential: int = 32,
        n_radial: int = 2,
        stress_to_force_factor: float = 1.0e6,
    ) -> None:
        if not isinstance(mesh, BeamMesh2D):
            raise TypeError("mesh must be a BeamMesh2D.")
        if not isinstance(
            tower_geometry,
            LinearTaperedTowerGeometry,
        ):
            raise TypeError(
                "tower_geometry must be a "
                "LinearTaperedTowerGeometry."
            )
        if not isinstance(material, MaterialParameters):
            raise TypeError(
                "material must be a MaterialParameters object."
            )

        n_gauss = _integer(n_gauss, "n_gauss")
        n_circumferential = _integer(
            n_circumferential,
            "n_circumferential",
        )
        n_radial = _integer(n_radial, "n_radial")
        stress_to_force_factor = _positive(
            stress_to_force_factor,
            "stress_to_force_factor",
        )

        if n_gauss < 1:
            raise ValueError("n_gauss must be at least 1.")
        if n_circumferential < 4:
            raise ValueError(
                "n_circumferential must be at least 4."
            )
        if n_circumferential % 4 != 0:
            raise ValueError(
                "n_circumferential must be divisible by 4."
            )
        if n_radial < 1:
            raise ValueError("n_radial must be at least 1.")

        self.mesh = mesh
        self.tower_geometry = tower_geometry
        self.material = material
        self.n_gauss = n_gauss
        self.n_circumferential = n_circumferential
        self.n_radial = n_radial
        self.stress_to_force_factor = stress_to_force_factor

        elements = []
        element_dofs = np.empty(
            (mesh.n_elements, 6),
            dtype=np.int64,
        )

        for element_index in range(mesh.n_elements):
            nodes = mesh.connectivity[element_index]
            dofs = element_dof_indices(nodes)
            element_dofs[element_index] = dofs

            elements.append(
                ViscoplasticDamageBeamElement2D(
                    node_coordinates=mesh.coordinates[nodes],
                    tower_axis_start=float(
                        mesh.tower_axis_coordinates[nodes[0]]
                    ),
                    tower_axis_end=float(
                        mesh.tower_axis_coordinates[nodes[1]]
                    ),
                    tower_geometry=tower_geometry,
                    material=material,
                    n_gauss=n_gauss,
                    n_circumferential=(
                        n_circumferential
                    ),
                    n_radial=n_radial,
                    stress_to_force_factor=(
                        stress_to_force_factor
                    ),
                )
            )

        self.elements = tuple(elements)
        self.element_dofs = element_dofs
        self._committed_displacements = np.zeros(
            mesh.n_dof,
            dtype=np.float64,
        )
        self._trial_displacements = (
            self._committed_displacements.copy()
        )
        self._has_uncommitted_trial = False

        self._committed_response = self._assemble_response(
            time=0.0,
            displacements=self._committed_displacements,
            element_responses=tuple(
                element.committed_response
                for element in self.elements
            ),
            compute_tangent=False,
        )
        self._trial_response = self._committed_response

    @property
    def committed_time(self) -> float:
        """Return the common committed time of all elements."""
        if not self.elements:
            return 0.0
        return float(self.elements[0].committed_time)

    @property
    def has_uncommitted_trial(self) -> bool:
        """Return whether a global trial state awaits resolution."""
        return bool(self._has_uncommitted_trial)

    @property
    def committed_displacements(self) -> FloatArray:
        """Return the last committed displacement vector."""
        return self._committed_displacements.copy()

    @property
    def trial_displacements(self) -> FloatArray:
        """Return the current trial displacement vector."""
        return self._trial_displacements.copy()

    @property
    def committed_response(self) -> ViscoplasticTowerResponse:
        """Return the last committed tower response."""
        return self._committed_response

    @property
    def trial_response(self) -> ViscoplasticTowerResponse:
        """Return the current tower trial response."""
        return self._trial_response


    def _assemble_response(
        self,
        time: float,
        displacements: FloatArray,
        element_responses: Tuple[
            ViscoplasticBeamElementResponse,
            ...,
        ],
        compute_tangent: bool,
    ) -> ViscoplasticTowerResponse:
        """Assemble global internal force and tangent."""
        displacements = _vector(
            displacements,
            self.mesh.n_dof,
            "displacements",
        )
        element_responses = tuple(element_responses)

        if len(element_responses) != self.mesh.n_elements:
            raise ValueError(
                "element_responses must contain one response "
                "per element."
            )

        internal_force = np.zeros(
            self.mesh.n_dof,
            dtype=np.float64,
        )
        tangent = (
            np.zeros(
                (self.mesh.n_dof, self.mesh.n_dof),
                dtype=np.float64,
            )
            if compute_tangent
            else None
        )

        for element_index, response in enumerate(
            element_responses
        ):
            dofs = self.element_dofs[element_index]
            internal_force[dofs] += (
                response.global_internal_force
            )

            if compute_tangent:
                if response.global_tangent is None:
                    raise RuntimeError(
                        "An element tangent is required for "
                        "global tangent assembly."
                    )
                tangent[np.ix_(dofs, dofs)] += (
                    response.global_tangent
                )

        if tangent is not None:
            tangent = 0.5 * (tangent + tangent.T)

        return ViscoplasticTowerResponse(
            time=time,
            displacements=displacements,
            internal_force=internal_force,
            tangent=tangent,
            element_dofs=self.element_dofs,
            element_responses=element_responses,
        )

    def set_trial_displacements(
        self,
        time: float,
        displacements: FloatArray,
        compute_tangent: bool = True,
        strain_perturbation: float = 1.0e-8,
    ) -> ViscoplasticTowerResponse:
        """Evaluate a global trial state from the last commit."""
        trial_time = _finite(time, "time")
        trial_displacements = _vector(
            displacements,
            self.mesh.n_dof,
            "displacements",
        )
        if not isinstance(
            compute_tangent,
            (bool, np.bool_),
        ):
            raise TypeError(
                "compute_tangent must be Boolean."
            )
        strain_perturbation = _positive(
            strain_perturbation,
            "strain_perturbation",
        )

        element_responses = []

        try:
            for element_index, element in enumerate(
                self.elements
            ):
                dofs = self.element_dofs[element_index]
                element_responses.append(
                    element.set_trial_displacements(
                        time=trial_time,
                        global_displacements=(
                            trial_displacements[dofs]
                        ),
                        compute_tangent=compute_tangent,
                        strain_perturbation=(
                            strain_perturbation
                        ),
                    )
                )
        except Exception:
            for element in self.elements:
                element.revert_to_last_commit()
            raise

        response = self._assemble_response(
            time=trial_time,
            displacements=trial_displacements,
            element_responses=tuple(
                element_responses
            ),
            compute_tangent=compute_tangent,
        )

        self._trial_displacements = (
            trial_displacements.copy()
        )
        self._trial_response = response
        self._has_uncommitted_trial = True

        return response

    def commit_state(self) -> ViscoplasticTowerResponse:
        """Commit every element after global convergence."""
        if not self._has_uncommitted_trial:
            raise RuntimeError(
                "No uncommitted tower trial state is available."
            )

        element_responses = tuple(
            element.commit_state()
            for element in self.elements
        )
        self._committed_displacements = (
            self._trial_displacements.copy()
        )
        self._committed_response = self._assemble_response(
            time=self.committed_time,
            displacements=self._committed_displacements,
            element_responses=element_responses,
            compute_tangent=(
                self._trial_response.tangent is not None
            ),
        )
        self._trial_displacements = (
            self._committed_displacements.copy()
        )
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False

        return self._committed_response

    def revert_to_last_commit(
        self,
    ) -> ViscoplasticTowerResponse:
        """Discard all current element trial states."""
        for element in self.elements:
            element.revert_to_last_commit()

        self._trial_displacements = (
            self._committed_displacements.copy()
        )
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False

        return self._trial_response

    def revert_to_start(self) -> ViscoplasticTowerResponse:
        """Restore zero displacement and all initial states."""
        element_responses = tuple(
            element.revert_to_start()
            for element in self.elements
        )
        self._committed_displacements = np.zeros(
            self.mesh.n_dof,
            dtype=np.float64,
        )
        self._trial_displacements = (
            self._committed_displacements.copy()
        )
        self._has_uncommitted_trial = False
        self._committed_response = self._assemble_response(
            time=0.0,
            displacements=self._committed_displacements,
            element_responses=element_responses,
            compute_tangent=False,
        )
        self._trial_response = self._committed_response

        return self._committed_response


def solve_nonlinear_tower_load_step(
    system: ViscoplasticDamageTowerSystem2D,
    time: float,
    load_vector: FloatArray,
    fixed_dofs: Optional[Sequence[int]] = None,
    initial_displacements: Optional[FloatArray] = None,
    max_iterations: int = 30,
    relative_residual_tolerance: float = 1.0e-8,
    absolute_residual_tolerance: float = 1.0e-6,
    strain_perturbation: float = 1.0e-8,
) -> NonlinearTowerStepSolution:
    """
    Solve and commit one force-controlled quasi-static load step.

    The Newton equation on the free DOFs is

        K_t * delta_u = f_ext - f_int.

    On failure, all trial element states are rolled back to the last
    committed state.
    """
    if not isinstance(
        system,
        ViscoplasticDamageTowerSystem2D,
    ):
        raise TypeError(
            "system must be a "
            "ViscoplasticDamageTowerSystem2D."
        )

    trial_time = _finite(time, "time")
    if trial_time <= system.committed_time:
        raise ValueError(
            "time must be greater than the committed time."
        )

    loads = _vector(
        load_vector,
        system.mesh.n_dof,
        "load_vector",
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

    if fixed_dofs is None:
        fixed = cantilever_base_fixed_dofs(system.mesh)
    else:
        fixed = _dof_array(
            fixed_dofs,
            system.mesh.n_dof,
            "fixed_dofs",
        )
    free = free_dofs_from_fixed(
        system.mesh.n_dof,
        fixed,
    )

    if initial_displacements is None:
        displacements = system.committed_displacements
    else:
        displacements = _vector(
            initial_displacements,
            system.mesh.n_dof,
            "initial_displacements",
        )
    displacements[fixed] = 0.0

    reference_force_norm = max(
        1.0,
        float(np.linalg.norm(loads[free])),
    )
    residual_limit = (
        absolute_residual_tolerance
        + relative_residual_tolerance
        * reference_force_norm
    )
    displacement_increment_norm = 0.0

    try:
        for iteration in range(
            1,
            max_iterations + 1,
        ):
            response = system.set_trial_displacements(
                time=trial_time,
                displacements=displacements,
                compute_tangent=True,
                strain_perturbation=(
                    strain_perturbation
                ),
            )
            residual = (
                loads - response.internal_force
            )
            residual_norm = float(
                np.linalg.norm(residual[free])
            )

            if residual_norm <= residual_limit:
                committed_response = (
                    system.commit_state()
                )
                reactions = (
                    committed_response.internal_force
                    - loads
                )
                final_residual = (
                    loads
                    - committed_response.internal_force
                )

                return NonlinearTowerStepSolution(
                    time=trial_time,
                    load_vector=loads,
                    displacements=(
                        committed_response.displacements
                    ),
                    reactions=reactions,
                    residual=final_residual,
                    fixed_dofs=fixed,
                    free_dofs=free,
                    iterations=iteration,
                    residual_norm=float(
                        np.linalg.norm(
                            final_residual[free]
                        )
                    ),
                    displacement_increment_norm=(
                        displacement_increment_norm
                    ),
                    response=committed_response,
                )

            if response.tangent is None:
                raise RuntimeError(
                    "A global tangent is required by "
                    "Newton's method."
                )

            reduced_tangent = (
                response.tangent[
                    np.ix_(free, free)
                ]
            )
            try:
                increment = np.linalg.solve(
                    reduced_tangent,
                    residual[free],
                )
            except np.linalg.LinAlgError as error:
                raise NonlinearTowerConvergenceError(
                    "The reduced nonlinear tangent is "
                    "singular or numerically "
                    "ill-conditioned."
                ) from error

            displacement_increment_norm = float(
                np.linalg.norm(increment)
            )
            displacements[free] += increment
            displacements[fixed] = 0.0

    except Exception:
        system.revert_to_last_commit()
        raise

    system.revert_to_last_commit()
    raise NonlinearTowerConvergenceError(
        "The nonlinear tower load step did not "
        "converge within "
        + str(max_iterations)
        + " Newton iterations."
    )

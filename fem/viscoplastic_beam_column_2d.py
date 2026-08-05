# -*- coding: utf-8 -*-
"""
Stateful 2D Euler-Bernoulli beam element with viscoplastic-damage
fiber sections at Gauss integration points.

Each Gauss point owns one independent
``ViscoplasticDamageFiberSection``. Repeated trial evaluations start
from the last committed state, so Newton or LATIN iterations do not
repeatedly accumulate plasticity or damage.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.beam_column_2d import (
    LinearTaperedTowerGeometry,
    euler_bernoulli_strain_displacement,
    gauss_legendre_rule,
    planar_frame_transformation,
)
from fem.fiber_section import create_annular_fiber_section
from fem.viscoplastic_fiber_section import (
    ViscoplasticDamageFiberSection,
    ViscoplasticSectionResponse,
)
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


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


def _integer(value: int, name: str) -> int:
    """Return a validated integer."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


def _vector6(values: FloatArray, name: str) -> FloatArray:
    """Return a finite six-component vector."""
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (6,):
        raise ValueError(name + " must have shape (6,).")
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain finite values.")
    return array.copy()


@dataclass(frozen=True)
class ViscoplasticBeamElementResponse:
    """Trial or committed response of one nonlinear beam element."""

    time: float
    global_displacements: FloatArray
    local_displacements: FloatArray
    generalized_strains: FloatArray
    section_resultants: FloatArray
    section_responses: Tuple[ViscoplasticSectionResponse, ...]
    local_internal_force: FloatArray
    global_internal_force: FloatArray
    local_tangent: Optional[FloatArray]
    global_tangent: Optional[FloatArray]

    def __post_init__(self) -> None:
        object.__setattr__(self, "time", _finite(self.time, "time"))

        global_displacements = _vector6(
            self.global_displacements,
            "global_displacements",
        )
        local_displacements = _vector6(
            self.local_displacements,
            "local_displacements",
        )
        local_internal_force = _vector6(
            self.local_internal_force,
            "local_internal_force",
        )
        global_internal_force = _vector6(
            self.global_internal_force,
            "global_internal_force",
        )

        generalized_strains = np.asarray(
            self.generalized_strains,
            dtype=np.float64,
        )
        section_resultants = np.asarray(
            self.section_resultants,
            dtype=np.float64,
        )

        if (
            generalized_strains.ndim != 2
            or generalized_strains.shape[1] != 2
        ):
            raise ValueError(
                "generalized_strains must have shape (n_gauss, 2)."
            )
        if section_resultants.shape != generalized_strains.shape:
            raise ValueError(
                "section_resultants must match generalized_strains."
            )
        if np.any(~np.isfinite(generalized_strains)):
            raise ValueError("generalized_strains must be finite.")
        if np.any(~np.isfinite(section_resultants)):
            raise ValueError("section_resultants must be finite.")

        section_responses = tuple(self.section_responses)
        if len(section_responses) != generalized_strains.shape[0]:
            raise ValueError(
                "section_responses must contain one response "
                "per Gauss point."
            )
        for response in section_responses:
            if not isinstance(response, ViscoplasticSectionResponse):
                raise TypeError(
                    "section_responses must contain "
                    "ViscoplasticSectionResponse objects."
                )

        local_tangent = self.local_tangent
        global_tangent = self.global_tangent
        if (local_tangent is None) != (global_tangent is None):
            raise ValueError(
                "local_tangent and global_tangent must both be "
                "present or both be None."
            )
        if local_tangent is not None:
            local_tangent = np.asarray(
                local_tangent,
                dtype=np.float64,
            )
            global_tangent = np.asarray(
                global_tangent,
                dtype=np.float64,
            )
            if local_tangent.shape != (6, 6):
                raise ValueError(
                    "local_tangent must have shape (6, 6)."
                )
            if global_tangent.shape != (6, 6):
                raise ValueError(
                    "global_tangent must have shape (6, 6)."
                )
            if np.any(~np.isfinite(local_tangent)):
                raise ValueError("local_tangent must be finite.")
            if np.any(~np.isfinite(global_tangent)):
                raise ValueError("global_tangent must be finite.")
            local_tangent = local_tangent.copy()
            global_tangent = global_tangent.copy()

        object.__setattr__(
            self,
            "global_displacements",
            global_displacements,
        )
        object.__setattr__(
            self,
            "local_displacements",
            local_displacements,
        )
        object.__setattr__(
            self,
            "generalized_strains",
            generalized_strains.copy(),
        )
        object.__setattr__(
            self,
            "section_resultants",
            section_resultants.copy(),
        )
        object.__setattr__(
            self,
            "section_responses",
            section_responses,
        )
        object.__setattr__(
            self,
            "local_internal_force",
            local_internal_force,
        )
        object.__setattr__(
            self,
            "global_internal_force",
            global_internal_force,
        )
        object.__setattr__(self, "local_tangent", local_tangent)
        object.__setattr__(self, "global_tangent", global_tangent)

    @property
    def axial_strains(self) -> FloatArray:
        """Return centroidal axial strains at all Gauss points."""
        return self.generalized_strains[:, 0].copy()

    @property
    def curvatures(self) -> FloatArray:
        """Return curvatures at all Gauss points."""
        return self.generalized_strains[:, 1].copy()

    @property
    def axial_forces(self) -> FloatArray:
        """Return axial forces at all Gauss points."""
        return self.section_resultants[:, 0].copy()

    @property
    def bending_moments(self) -> FloatArray:
        """Return bending moments at all Gauss points."""
        return self.section_resultants[:, 1].copy()

    @property
    def maximum_damage(self) -> float:
        """Return the largest fiber damage over all Gauss points."""
        return float(
            max(
                response.maximum_damage
                for response in self.section_responses
            )
        )


class ViscoplasticDamageBeamElement2D:
    """
    Stateful two-node planar Euler-Bernoulli beam element.

    Local element DOF ordering is

        [u1, v1, theta1, u2, v2, theta2].

    The transformation convention is

        u_local = T @ u_global.
    """

    def __init__(
        self,
        node_coordinates: FloatArray,
        tower_axis_start: float,
        tower_axis_end: float,
        tower_geometry: LinearTaperedTowerGeometry,
        material: MaterialParameters,
        n_gauss: int = 4,
        n_circumferential: int = 32,
        n_radial: int = 2,
        stress_to_force_factor: float = 1.0e6,
    ) -> None:
        coordinates = np.asarray(
            node_coordinates,
            dtype=np.float64,
        )
        if coordinates.shape != (2, 2):
            raise ValueError(
                "node_coordinates must have shape (2, 2)."
            )
        if np.any(~np.isfinite(coordinates)):
            raise ValueError(
                "node_coordinates must contain finite values."
            )
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

        tower_axis_start = _finite(
            tower_axis_start,
            "tower_axis_start",
        )
        tower_axis_end = _finite(
            tower_axis_end,
            "tower_axis_end",
        )
        if tower_axis_end <= tower_axis_start:
            raise ValueError(
                "tower_axis_end must be greater than "
                "tower_axis_start."
            )

        tolerance = (
            1.0e-12 * max(1.0, tower_geometry.height)
        )
        if tower_axis_start < -tolerance:
            raise ValueError(
                "tower_axis_start must lie within the tower."
            )
        if (
            tower_axis_end
            > tower_geometry.height + tolerance
        ):
            raise ValueError(
                "tower_axis_end must lie within the tower."
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

        length = float(
            np.linalg.norm(
                coordinates[1] - coordinates[0]
            )
        )
        if length <= 0.0:
            raise ValueError(
                "Beam-element node coordinates must be distinct."
            )

        self.node_coordinates = coordinates.copy()
        self.tower_axis_start = tower_axis_start
        self.tower_axis_end = tower_axis_end
        self.tower_geometry = tower_geometry
        self.material = material
        self.n_gauss = n_gauss
        self.n_circumferential = n_circumferential
        self.n_radial = n_radial
        self.stress_to_force_factor = (
            stress_to_force_factor
        )
        self.length = length
        self.jacobian = 0.5 * length

        self.transformation = (
            planar_frame_transformation(coordinates)
        )
        (
            self.gauss_coordinates,
            self.gauss_weights,
        ) = gauss_legendre_rule(n_gauss)

        self.gauss_heights = (
            0.5
            * (
                tower_axis_start
                + tower_axis_end
            )
            + 0.5
            * (
                tower_axis_end
                - tower_axis_start
            )
            * self.gauss_coordinates
        )

        sections = []
        b_matrices = []

        for gauss_index in range(n_gauss):
            height_coordinate = float(
                self.gauss_heights[gauss_index]
            )
            geometry_section = (
                create_annular_fiber_section(
                    outer_diameter=(
                        tower_geometry.outer_diameter_at(
                            height_coordinate
                        )
                    ),
                    thickness=(
                        tower_geometry.thickness_at(
                            height_coordinate
                        )
                    ),
                    n_circumferential=(
                        n_circumferential
                    ),
                    n_radial=n_radial,
                )
            )
            sections.append(
                ViscoplasticDamageFiberSection(
                    section=geometry_section,
                    material=material,
                    stress_to_force_factor=(
                        stress_to_force_factor
                    ),
                )
            )
            b_matrices.append(
                euler_bernoulli_strain_displacement(
                    length=length,
                    natural_coordinate=float(
                        self.gauss_coordinates[
                            gauss_index
                        ]
                    ),
                )
            )

        self.sections = tuple(sections)
        self.b_matrices = np.stack(
            b_matrices,
            axis=0,
        )

        zero_displacements = np.zeros(
            6,
            dtype=np.float64,
        )
        self._committed_global_displacements = (
            zero_displacements.copy()
        )
        self._trial_global_displacements = (
            zero_displacements.copy()
        )
        self._has_uncommitted_trial = False

        self._committed_response = (
            self._response_from_sections(
                time=0.0,
                global_displacements=zero_displacements,
                section_responses=tuple(
                    section.committed_response
                    for section in self.sections
                ),
                compute_tangent=False,
            )
        )
        self._trial_response = self._committed_response

    @property
    def committed_time(self) -> float:
        """Return the common committed time of all sections."""
        return float(self.sections[0].committed_time)

    @property
    def trial_time(self) -> float:
        """Return the common trial time of all sections."""
        return float(self.sections[0].trial_time)

    @property
    def has_uncommitted_trial(self) -> bool:
        """Return whether a trial state awaits commit or rollback."""
        return bool(self._has_uncommitted_trial)

    @property
    def committed_global_displacements(self) -> FloatArray:
        """Return the last committed displacement vector."""
        return self._committed_global_displacements.copy()

    @property
    def trial_global_displacements(self) -> FloatArray:
        """Return the current trial displacement vector."""
        return self._trial_global_displacements.copy()

    @property
    def committed_response(
        self,
    ) -> ViscoplasticBeamElementResponse:
        """Return the last committed element response."""
        return self._committed_response

    @property
    def trial_response(
        self,
    ) -> ViscoplasticBeamElementResponse:
        """Return the current trial element response."""
        return self._trial_response


    def _response_from_sections(
        self,
        time: float,
        global_displacements: FloatArray,
        section_responses: Tuple[
            ViscoplasticSectionResponse,
            ...,
        ],
        compute_tangent: bool,
    ) -> ViscoplasticBeamElementResponse:
        """Integrate known section responses over the element."""
        global_displacements = _vector6(
            global_displacements,
            "global_displacements",
        )
        section_responses = tuple(section_responses)

        if len(section_responses) != self.n_gauss:
            raise ValueError(
                "section_responses must contain one response "
                "per Gauss point."
            )

        local_displacements = (
            self.transformation @ global_displacements
        )

        generalized_strains = np.zeros(
            (self.n_gauss, 2),
            dtype=np.float64,
        )
        section_resultants = np.zeros(
            (self.n_gauss, 2),
            dtype=np.float64,
        )
        local_internal_force = np.zeros(
            6,
            dtype=np.float64,
        )
        local_tangent = (
            np.zeros((6, 6), dtype=np.float64)
            if compute_tangent
            else None
        )

        for gauss_index in range(self.n_gauss):
            b_matrix = self.b_matrices[gauss_index]
            response = section_responses[gauss_index]
            integration_factor = float(
                self.gauss_weights[gauss_index]
                * self.jacobian
            )

            generalized_strains[gauss_index] = (
                b_matrix @ local_displacements
            )
            section_resultants[gauss_index] = (
                response.resultants
            )
            local_internal_force += (
                b_matrix.T
                @ response.resultants
                * integration_factor
            )

            if compute_tangent:
                if response.tangent is None:
                    raise RuntimeError(
                        "A section tangent is required for "
                        "element tangent integration."
                    )
                local_tangent += (
                    b_matrix.T
                    @ response.tangent
                    @ b_matrix
                    * integration_factor
                )

        if local_tangent is not None:
            local_tangent = 0.5 * (
                local_tangent
                + local_tangent.T
            )

        global_internal_force = (
            self.transformation.T
            @ local_internal_force
        )

        global_tangent = None
        if local_tangent is not None:
            global_tangent = (
                self.transformation.T
                @ local_tangent
                @ self.transformation
            )
            global_tangent = 0.5 * (
                global_tangent
                + global_tangent.T
            )

        return ViscoplasticBeamElementResponse(
            time=time,
            global_displacements=global_displacements,
            local_displacements=local_displacements,
            generalized_strains=generalized_strains,
            section_resultants=section_resultants,
            section_responses=section_responses,
            local_internal_force=local_internal_force,
            global_internal_force=global_internal_force,
            local_tangent=local_tangent,
            global_tangent=global_tangent,
        )

    def set_trial_displacements(
        self,
        time: float,
        global_displacements: FloatArray,
        compute_tangent: bool = True,
        strain_perturbation: float = 1.0e-8,
    ) -> ViscoplasticBeamElementResponse:
        """
        Evaluate one trial step from the last committed state.

        Repeated calls before ``commit_state()`` restart every
        Gauss-point fiber section from the same committed state.
        """
        trial_time = _finite(time, "time")
        trial_global_displacements = _vector6(
            global_displacements,
            "global_displacements",
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

        local_displacements = (
            self.transformation
            @ trial_global_displacements
        )

        section_responses = []

        try:
            for gauss_index in range(self.n_gauss):
                generalized_strain = (
                    self.b_matrices[gauss_index]
                    @ local_displacements
                )
                section_responses.append(
                    self.sections[
                        gauss_index
                    ].set_trial_deformation(
                        time=trial_time,
                        axial_strain=float(
                            generalized_strain[0]
                        ),
                        curvature=float(
                            generalized_strain[1]
                        ),
                        compute_tangent=compute_tangent,
                        strain_perturbation=(
                            strain_perturbation
                        ),
                    )
                )
        except Exception:
            for section in self.sections:
                section.revert_to_last_commit()
            raise

        response = self._response_from_sections(
            time=trial_time,
            global_displacements=(
                trial_global_displacements
            ),
            section_responses=tuple(
                section_responses
            ),
            compute_tangent=compute_tangent,
        )

        self._trial_global_displacements = (
            trial_global_displacements.copy()
        )
        self._trial_response = response
        self._has_uncommitted_trial = True

        return response


    def commit_state(
        self,
    ) -> ViscoplasticBeamElementResponse:
        """Commit the current trial state at every Gauss point."""
        if not self._has_uncommitted_trial:
            raise RuntimeError(
                "No uncommitted trial element state is available."
            )

        committed_section_responses = tuple(
            section.commit_state()
            for section in self.sections
        )

        self._committed_global_displacements = (
            self._trial_global_displacements.copy()
        )
        self._committed_response = (
            self._response_from_sections(
                time=self.committed_time,
                global_displacements=(
                    self._committed_global_displacements
                ),
                section_responses=(
                    committed_section_responses
                ),
                compute_tangent=(
                    self._trial_response.local_tangent
                    is not None
                ),
            )
        )

        self._trial_global_displacements = (
            self._committed_global_displacements.copy()
        )
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False

        return self._committed_response

    def revert_to_last_commit(
        self,
    ) -> ViscoplasticBeamElementResponse:
        """Discard the current trial state at every Gauss point."""
        for section in self.sections:
            section.revert_to_last_commit()

        self._trial_global_displacements = (
            self._committed_global_displacements.copy()
        )
        self._trial_response = self._committed_response
        self._has_uncommitted_trial = False

        return self._trial_response

    def revert_to_start(
        self,
    ) -> ViscoplasticBeamElementResponse:
        """Restore zero displacement and all initial fiber states."""
        for section in self.sections:
            section.revert_to_start()

        zero_displacements = np.zeros(
            6,
            dtype=np.float64,
        )
        self._committed_global_displacements = (
            zero_displacements.copy()
        )
        self._trial_global_displacements = (
            zero_displacements.copy()
        )
        self._has_uncommitted_trial = False

        self._committed_response = (
            self._response_from_sections(
                time=0.0,
                global_displacements=zero_displacements,
                section_responses=tuple(
                    section.committed_response
                    for section in self.sections
                ),
                compute_tangent=False,
            )
        )
        self._trial_response = self._committed_response

        return self._committed_response

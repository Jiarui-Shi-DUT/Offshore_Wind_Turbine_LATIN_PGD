# -*- coding: utf-8 -*-
"""
Global assembly and linear-static solution utilities for the
two-dimensional tapered offshore-wind-turbine tower model.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from fem.beam_column_2d import (
    BeamMesh2D,
    ElasticBeamElementResult,
    LinearTaperedTowerGeometry,
    compute_elastic_beam_element_stiffness,
)


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _validated_finite_scalar(value: float, name: str) -> float:
    """Return a finite floating-point scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error

    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _validated_positive_scalar(value: float, name: str) -> float:
    """Return a finite positive floating-point scalar."""
    result = _validated_finite_scalar(value, name)
    if result <= 0.0:
        raise ValueError(name + " must be positive.")
    return result


def _validated_integer(value: int, name: str) -> int:
    """Return a validated integer parameter."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


def element_dof_indices(
    node_indices: Sequence[int],
) -> IntArray:
    """
    Return the six global DOF indices of a two-node planar beam element.

    The global nodal DOF ordering is

        [u_x, u_y, theta_z].
    """
    nodes = np.asarray(node_indices, dtype=np.int64)

    if nodes.shape != (2,):
        raise ValueError("node_indices must contain exactly two nodes.")
    if np.any(nodes < 0):
        raise ValueError("node_indices must be non-negative.")
    if nodes[0] == nodes[1]:
        raise ValueError("The two element nodes must be distinct.")

    dofs = np.empty(6, dtype=np.int64)
    dofs[:3] = 3 * nodes[0] + np.arange(3, dtype=np.int64)
    dofs[3:] = 3 * nodes[1] + np.arange(3, dtype=np.int64)
    return dofs


@dataclass(frozen=True)
class ElasticTowerAssembly:
    """Global stiffness matrix and element-level data."""

    mesh: BeamMesh2D
    global_stiffness: FloatArray
    element_dofs: IntArray
    element_results: Tuple[ElasticBeamElementResult, ...]

    def __post_init__(self) -> None:
        if not isinstance(self.mesh, BeamMesh2D):
            raise TypeError("mesh must be a BeamMesh2D.")

        global_stiffness = np.asarray(
            self.global_stiffness,
            dtype=np.float64,
        )
        element_dofs = np.asarray(
            self.element_dofs,
            dtype=np.int64,
        )
        element_results = tuple(self.element_results)

        expected_matrix_shape = (self.mesh.n_dof, self.mesh.n_dof)
        if global_stiffness.shape != expected_matrix_shape:
            raise ValueError(
                "global_stiffness must have shape "
                + str(expected_matrix_shape)
                + "."
            )
        if np.any(~np.isfinite(global_stiffness)):
            raise ValueError("global_stiffness must be finite.")

        if element_dofs.shape != (self.mesh.n_elements, 6):
            raise ValueError(
                "element_dofs must have shape (n_elements, 6)."
            )
        if np.any(element_dofs < 0):
            raise ValueError("element_dofs must be non-negative.")
        if np.any(element_dofs >= self.mesh.n_dof):
            raise ValueError(
                "element_dofs contain an out-of-range DOF index."
            )

        if len(element_results) != self.mesh.n_elements:
            raise ValueError(
                "element_results must contain one result per element."
            )
        for result in element_results:
            if not isinstance(result, ElasticBeamElementResult):
                raise TypeError(
                    "Every element result must be an "
                    "ElasticBeamElementResult."
                )

        object.__setattr__(
            self,
            "global_stiffness",
            global_stiffness,
        )
        object.__setattr__(self, "element_dofs", element_dofs)
        object.__setattr__(
            self,
            "element_results",
            element_results,
        )


@dataclass(frozen=True)
class LinearStaticSolution:
    """Displacements, reactions, and DOF partitions."""

    displacements: FloatArray
    reactions: FloatArray
    fixed_dofs: IntArray
    free_dofs: IntArray

    def __post_init__(self) -> None:
        displacements = np.asarray(
            self.displacements,
            dtype=np.float64,
        )
        reactions = np.asarray(
            self.reactions,
            dtype=np.float64,
        )
        fixed_dofs = np.asarray(
            self.fixed_dofs,
            dtype=np.int64,
        )
        free_dofs = np.asarray(
            self.free_dofs,
            dtype=np.int64,
        )

        if displacements.ndim != 1:
            raise ValueError("displacements must be one-dimensional.")
        if reactions.shape != displacements.shape:
            raise ValueError(
                "reactions must have the same shape as displacements."
            )
        if np.any(~np.isfinite(displacements)):
            raise ValueError("displacements must be finite.")
        if np.any(~np.isfinite(reactions)):
            raise ValueError("reactions must be finite.")

        if fixed_dofs.ndim != 1 or free_dofs.ndim != 1:
            raise ValueError(
                "fixed_dofs and free_dofs must be one-dimensional."
            )
        if np.intersect1d(fixed_dofs, free_dofs).size != 0:
            raise ValueError(
                "fixed_dofs and free_dofs must not overlap."
            )

        all_dofs = np.concatenate((fixed_dofs, free_dofs))
        if all_dofs.size != displacements.size:
            raise ValueError(
                "fixed_dofs and free_dofs must partition all DOFs."
            )
        if np.unique(all_dofs).size != displacements.size:
            raise ValueError(
                "fixed_dofs and free_dofs must contain unique DOFs."
            )
        if np.any(all_dofs < 0):
            raise ValueError("DOF indices must be non-negative.")
        if np.any(all_dofs >= displacements.size):
            raise ValueError("A DOF index is out of range.")

        object.__setattr__(self, "displacements", displacements)
        object.__setattr__(self, "reactions", reactions)
        object.__setattr__(self, "fixed_dofs", fixed_dofs)
        object.__setattr__(self, "free_dofs", free_dofs)


def assemble_elastic_tower_stiffness(
    mesh: BeamMesh2D,
    tower_geometry: LinearTaperedTowerGeometry,
    elastic_modulus: float,
    n_gauss: int = 4,
    n_circumferential: int = 32,
    n_radial: int = 2,
) -> ElasticTowerAssembly:
    """
    Assemble the global elastic stiffness matrix of the tower.

    Every beam element is integrated independently. Local annular
    fiber sections are rebuilt at each axial Gauss point.
    """
    if not isinstance(mesh, BeamMesh2D):
        raise TypeError("mesh must be a BeamMesh2D.")
    if not isinstance(
        tower_geometry,
        LinearTaperedTowerGeometry,
    ):
        raise TypeError(
            "tower_geometry must be a LinearTaperedTowerGeometry."
        )

    elastic_modulus = _validated_positive_scalar(
        elastic_modulus,
        "elastic_modulus",
    )
    n_gauss = _validated_integer(n_gauss, "n_gauss")
    n_circumferential = _validated_integer(
        n_circumferential,
        "n_circumferential",
    )
    n_radial = _validated_integer(n_radial, "n_radial")

    if n_gauss < 1:
        raise ValueError("n_gauss must be at least 1.")
    if n_circumferential < 1:
        raise ValueError(
            "n_circumferential must be at least 1."
        )
    if n_radial < 1:
        raise ValueError("n_radial must be at least 1.")

    global_stiffness = np.zeros(
        (mesh.n_dof, mesh.n_dof),
        dtype=np.float64,
    )
    element_dofs = np.empty(
        (mesh.n_elements, 6),
        dtype=np.int64,
    )
    element_results = []

    for element_index in range(mesh.n_elements):
        nodes = mesh.connectivity[element_index]
        dofs = element_dof_indices(nodes)
        node_coordinates = mesh.coordinates[nodes]

        tower_axis_start = float(
            mesh.tower_axis_coordinates[nodes[0]]
        )
        tower_axis_end = float(
            mesh.tower_axis_coordinates[nodes[1]]
        )

        result = compute_elastic_beam_element_stiffness(
            node_coordinates=node_coordinates,
            tower_axis_start=tower_axis_start,
            tower_axis_end=tower_axis_end,
            tower_geometry=tower_geometry,
            elastic_modulus=elastic_modulus,
            n_gauss=n_gauss,
            n_circumferential=n_circumferential,
            n_radial=n_radial,
        )

        global_stiffness[np.ix_(dofs, dofs)] += (
            result.global_stiffness
        )
        element_dofs[element_index] = dofs
        element_results.append(result)

    global_stiffness = 0.5 * (
        global_stiffness + global_stiffness.T
    )

    return ElasticTowerAssembly(
        mesh=mesh,
        global_stiffness=global_stiffness,
        element_dofs=element_dofs,
        element_results=tuple(element_results),
    )


def cantilever_base_fixed_dofs(
    mesh: BeamMesh2D,
) -> IntArray:
    """
    Return the fully fixed DOFs of the first tower node.

    The first node is the tower base, so its horizontal displacement,
    vertical displacement, and rotation are all constrained.
    """
    if not isinstance(mesh, BeamMesh2D):
        raise TypeError("mesh must be a BeamMesh2D.")

    return np.array([0, 1, 2], dtype=np.int64)


def top_horizontal_load_vector(
    mesh: BeamMesh2D,
    horizontal_force: float,
) -> FloatArray:
    """
    Create a global nodal load vector with a tower-top horizontal force.

    Positive force acts in the positive global x direction.
    """
    if not isinstance(mesh, BeamMesh2D):
        raise TypeError("mesh must be a BeamMesh2D.")

    horizontal_force = _validated_finite_scalar(
        horizontal_force,
        "horizontal_force",
    )

    load_vector = np.zeros(mesh.n_dof, dtype=np.float64)
    top_node = mesh.n_nodes - 1
    top_horizontal_dof = 3 * top_node
    load_vector[top_horizontal_dof] = horizontal_force
    return load_vector


def free_dofs_from_fixed(
    n_dof: int,
    fixed_dofs: Sequence[int],
) -> IntArray:
    """Return the sorted complement of the prescribed DOFs."""
    n_dof = _validated_integer(n_dof, "n_dof")
    if n_dof < 1:
        raise ValueError("n_dof must be at least 1.")

    fixed = np.asarray(fixed_dofs, dtype=np.int64)
    if fixed.ndim != 1:
        raise ValueError("fixed_dofs must be one-dimensional.")
    if np.any(fixed < 0) or np.any(fixed >= n_dof):
        raise ValueError("fixed_dofs contain an out-of-range index.")
    if np.unique(fixed).size != fixed.size:
        raise ValueError("fixed_dofs must not contain duplicates.")

    mask = np.ones(n_dof, dtype=bool)
    mask[fixed] = False
    return np.flatnonzero(mask).astype(np.int64)


def solve_linear_static_system(
    global_stiffness: FloatArray,
    load_vector: FloatArray,
    fixed_dofs: Sequence[int],
) -> LinearStaticSolution:
    """
    Solve a linear static finite-element system by direct elimination.

    The prescribed displacement value is zero at every fixed DOF.
    Reactions are recovered from K u - f after solving the reduced
    free-DOF system.
    """
    stiffness = np.asarray(
        global_stiffness,
        dtype=np.float64,
    )
    loads = np.asarray(load_vector, dtype=np.float64)

    if stiffness.ndim != 2:
        raise ValueError(
            "global_stiffness must be a two-dimensional matrix."
        )
    if stiffness.shape[0] != stiffness.shape[1]:
        raise ValueError(
            "global_stiffness must be square."
        )
    if np.any(~np.isfinite(stiffness)):
        raise ValueError("global_stiffness must be finite.")

    n_dof = stiffness.shape[0]
    if loads.shape != (n_dof,):
        raise ValueError(
            "load_vector must have shape (n_dof,)."
        )
    if np.any(~np.isfinite(loads)):
        raise ValueError("load_vector must be finite.")

    fixed = np.asarray(fixed_dofs, dtype=np.int64)
    free = free_dofs_from_fixed(n_dof, fixed)

    displacements = np.zeros(n_dof, dtype=np.float64)

    if free.size > 0:
        reduced_stiffness = stiffness[np.ix_(free, free)]
        reduced_loads = loads[free]

        try:
            displacements[free] = np.linalg.solve(
                reduced_stiffness,
                reduced_loads,
            )
        except np.linalg.LinAlgError as error:
            raise np.linalg.LinAlgError(
                "The reduced stiffness matrix is singular or "
                "numerically ill-conditioned."
            ) from error

    reactions = stiffness @ displacements - loads

    return LinearStaticSolution(
        displacements=displacements,
        reactions=reactions,
        fixed_dofs=np.sort(fixed),
        free_dofs=free,
    )


def solve_elastic_tower_top_force(
    mesh: BeamMesh2D,
    tower_geometry: LinearTaperedTowerGeometry,
    elastic_modulus: float,
    horizontal_force: float,
    n_gauss: int = 4,
    n_circumferential: int = 32,
    n_radial: int = 2,
) -> Tuple[ElasticTowerAssembly, LinearStaticSolution]:
    """
    Assemble and solve the fixed-base tower under one top horizontal force.
    """
    assembly = assemble_elastic_tower_stiffness(
        mesh=mesh,
        tower_geometry=tower_geometry,
        elastic_modulus=elastic_modulus,
        n_gauss=n_gauss,
        n_circumferential=n_circumferential,
        n_radial=n_radial,
    )
    loads = top_horizontal_load_vector(
        mesh=mesh,
        horizontal_force=horizontal_force,
    )
    fixed_dofs = cantilever_base_fixed_dofs(mesh)
    solution = solve_linear_static_system(
        global_stiffness=assembly.global_stiffness,
        load_vector=loads,
        fixed_dofs=fixed_dofs,
    )
    return assembly, solution

# -*- coding: utf-8 -*-
"""
Mode-wise SVD diagnostics for LATIN-PGD tower cycle-phase snapshots.

This module performs empirical low-rank diagnostics only.  It does not solve
the LATIN equations and it does not construct a PGD approximation.

For one field q(n, tau, x), the trailing spatial dimensions are first flattened
into one spatial coordinate s.  Three matrix unfoldings are then inspected:

    slow-cycle mode : n   x (tau, s)
    fast-phase mode : tau x (n, s)
    spatial mode    : s   x (n, tau)

The singular-value decay of these unfoldings provides an HOSVD-style
diagnostic of multilinear compressibility.  These matrix ranks are not the
same as the CP/PGD separation rank, so they should be interpreted as empirical
evidence about separability rather than as a direct PGD rank estimate.

Two versions of each tower field are diagnosed:

    raw
        q(n, tau, s)

    cycle_increment
        q(n, tau, s) - q(n, 0, s)

The cycle-increment form removes the slowly drifting cycle-start baseline while
retaining the complete within-cycle path, including the non-zero cycle-end
increment caused by ratcheting or damage accumulation.  This is especially
useful for distinguishing slow evolution from changes in the fast cyclic
waveform.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Tuple

import numpy as np
from numpy.typing import NDArray

from examples.nonlinear_tower_snapshot_tensor import (
    TowerCyclePhaseSnapshots,
)


FloatArray = NDArray[np.float64]
Shape3D = Tuple[int, int, int]


def flatten_cycle_phase_field(
    values: FloatArray,
    name: str = "field",
) -> FloatArray:
    """
    Return one field as (n_cycles, n_phase_points, n_spatial).

    The first two axes are preserved exactly.  All trailing axes are flattened
    into one spatial coordinate.
    """
    array = np.asarray(values, dtype=np.float64)

    if array.ndim < 3:
        raise ValueError(
            name
            + " must contain cycle, phase, and at least one spatial axis."
        )
    if array.shape[0] < 1:
        raise ValueError(name + " must contain at least one cycle.")
    if array.shape[1] < 2:
        raise ValueError(
            name + " must contain at least two phase points."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain only finite values.")

    return np.asarray(
        array.reshape(
            array.shape[0],
            array.shape[1],
            -1,
        ),
        dtype=np.float64,
    )


def cycle_increment_field(
    values: FloatArray,
    name: str = "field",
) -> FloatArray:
    """
    Return q(n, tau, s) - q(n, 0, s).

    The cycle-end value is not forced to zero.  Consequently, cycle-end drift
    and irreversible accumulation remain present in the result.
    """
    tensor = flatten_cycle_phase_field(values=values, name=name)
    return tensor - tensor[:, 0:1, :]


def mode_unfolding(
    tensor: FloatArray,
    mode: str,
) -> FloatArray:
    """
    Unfold a three-dimensional (cycle, phase, space) tensor along one mode.

    Supported modes are ``cycle``, ``phase``, and ``space``.
    """
    array = np.asarray(tensor, dtype=np.float64)

    if array.ndim != 3:
        raise ValueError(
            "tensor must have shape (n_cycles, n_phase_points, n_spatial)."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError("tensor must contain only finite values.")

    n_cycles, n_phase_points, n_spatial = array.shape

    if mode == "cycle":
        return np.asarray(
            array.reshape(n_cycles, n_phase_points * n_spatial),
            dtype=np.float64,
        )

    if mode == "phase":
        return np.asarray(
            np.transpose(array, (1, 0, 2)).reshape(
                n_phase_points,
                n_cycles * n_spatial,
            ),
            dtype=np.float64,
        )

    if mode == "space":
        return np.asarray(
            np.transpose(array, (2, 0, 1)).reshape(
                n_spatial,
                n_cycles * n_phase_points,
            ),
            dtype=np.float64,
        )

    raise ValueError(
        "mode must be one of: 'cycle', 'phase', or 'space'."
    )


@dataclass(frozen=True)
class SvdSpectrum:
    """Singular-value energy spectrum of one matrix unfolding."""

    singular_values: FloatArray
    energy_fractions: FloatArray
    cumulative_energy: FloatArray
    relative_truncation_errors: FloatArray
    frobenius_norm: float

    def __post_init__(self) -> None:
        singular_values = np.asarray(
            self.singular_values,
            dtype=np.float64,
        )
        energy_fractions = np.asarray(
            self.energy_fractions,
            dtype=np.float64,
        )
        cumulative_energy = np.asarray(
            self.cumulative_energy,
            dtype=np.float64,
        )
        relative_truncation_errors = np.asarray(
            self.relative_truncation_errors,
            dtype=np.float64,
        )
        frobenius_norm = float(self.frobenius_norm)

        if singular_values.ndim != 1:
            raise ValueError(
                "singular_values must be one-dimensional."
            )
        if singular_values.size < 1:
            raise ValueError(
                "At least one singular value is required."
            )
        expected_shape = singular_values.shape
        for array, name in (
            (energy_fractions, "energy_fractions"),
            (cumulative_energy, "cumulative_energy"),
            (
                relative_truncation_errors,
                "relative_truncation_errors",
            ),
        ):
            if array.shape != expected_shape:
                raise ValueError(
                    name + " must match singular_values."
                )

        if np.any(~np.isfinite(singular_values)):
            raise ValueError(
                "singular_values must contain finite values."
            )
        if np.any(singular_values < 0.0):
            raise ValueError(
                "singular_values must be non-negative."
            )
        if np.any(~np.isfinite(energy_fractions)):
            raise ValueError(
                "energy_fractions must contain finite values."
            )
        if np.any(~np.isfinite(cumulative_energy)):
            raise ValueError(
                "cumulative_energy must contain finite values."
            )
        if np.any(~np.isfinite(relative_truncation_errors)):
            raise ValueError(
                "relative_truncation_errors must contain finite values."
            )
        if not np.isfinite(frobenius_norm):
            raise ValueError("frobenius_norm must be finite.")
        if frobenius_norm < 0.0:
            raise ValueError("frobenius_norm must be non-negative.")

        object.__setattr__(
            self,
            "singular_values",
            singular_values.copy(),
        )
        object.__setattr__(
            self,
            "energy_fractions",
            energy_fractions.copy(),
        )
        object.__setattr__(
            self,
            "cumulative_energy",
            cumulative_energy.copy(),
        )
        object.__setattr__(
            self,
            "relative_truncation_errors",
            relative_truncation_errors.copy(),
        )
        object.__setattr__(
            self,
            "frobenius_norm",
            frobenius_norm,
        )

    @property
    def is_zero(self) -> bool:
        """Return whether the analyzed matrix is exactly zero in norm."""
        return self.frobenius_norm == 0.0

    def rank_for_energy(self, target_energy: float) -> int:
        """
        Return the smallest rank reaching the requested squared-norm energy.

        For an exactly zero matrix, rank zero reconstructs the field exactly.
        """
        try:
            target = float(target_energy)
        except (TypeError, ValueError) as error:
            raise TypeError(
                "target_energy must be a real scalar."
            ) from error

        if not np.isfinite(target):
            raise ValueError("target_energy must be finite.")
        if target <= 0.0 or target > 1.0:
            raise ValueError(
                "target_energy must satisfy 0 < target_energy <= 1."
            )

        if self.is_zero:
            return 0

        index = int(
            np.searchsorted(
                self.cumulative_energy,
                target,
                side="left",
            )
        )
        return min(index + 1, self.singular_values.size)

    def error_after_rank(self, rank: int) -> float:
        """
        Return the optimal relative Frobenius truncation error after ``rank``.

        Rank zero is allowed and corresponds to the zero approximation.
        """
        if isinstance(rank, (bool, np.bool_)):
            raise TypeError("rank must be an integer.")
        if not isinstance(rank, (int, np.integer)):
            raise TypeError("rank must be an integer.")

        rank = int(rank)
        if rank < 0 or rank > self.singular_values.size:
            raise ValueError(
                "rank must lie between zero and the number "
                "of singular values."
            )

        if self.is_zero:
            return 0.0
        if rank == 0:
            return 1.0
        return float(self.relative_truncation_errors[rank - 1])


def compute_svd_spectrum(matrix: FloatArray) -> SvdSpectrum:
    """Compute singular values and squared-Frobenius energy diagnostics."""
    array = np.asarray(matrix, dtype=np.float64)

    if array.ndim != 2:
        raise ValueError("matrix must be two-dimensional.")
    if array.shape[0] < 1 or array.shape[1] < 1:
        raise ValueError("matrix dimensions must be positive.")
    if np.any(~np.isfinite(array)):
        raise ValueError("matrix must contain only finite values.")

    singular_values = np.linalg.svd(
        array,
        full_matrices=False,
        compute_uv=False,
    )
    singular_values = np.asarray(
        singular_values,
        dtype=np.float64,
    )

    squared = singular_values * singular_values
    total_energy = float(np.sum(squared))
    frobenius_norm = float(np.sqrt(total_energy))

    if total_energy == 0.0:
        energy_fractions = np.zeros_like(singular_values)
        cumulative_energy = np.zeros_like(singular_values)
        relative_errors = np.zeros_like(singular_values)
    else:
        energy_fractions = squared / total_energy
        cumulative_energy = np.cumsum(energy_fractions)
        cumulative_energy[-1] = 1.0
        relative_errors = np.sqrt(
            np.maximum(0.0, 1.0 - cumulative_energy)
        )

    return SvdSpectrum(
        singular_values=singular_values,
        energy_fractions=energy_fractions,
        cumulative_energy=cumulative_energy,
        relative_truncation_errors=relative_errors,
        frobenius_norm=frobenius_norm,
    )


@dataclass(frozen=True)
class TensorSvdDiagnostics:
    """Mode-wise SVD diagnostics of one (cycle, phase, space) tensor."""

    tensor_shape: Shape3D
    cycle_mode: SvdSpectrum
    phase_mode: SvdSpectrum
    space_mode: SvdSpectrum

    def __post_init__(self) -> None:
        shape = tuple(int(value) for value in self.tensor_shape)
        if len(shape) != 3:
            raise ValueError(
                "tensor_shape must contain exactly three dimensions."
            )
        if any(value < 1 for value in shape):
            raise ValueError(
                "tensor_shape dimensions must all be positive."
            )

        object.__setattr__(self, "tensor_shape", shape)

    def mode(self, name: str) -> SvdSpectrum:
        """Return one named mode spectrum."""
        if name == "cycle":
            return self.cycle_mode
        if name == "phase":
            return self.phase_mode
        if name == "space":
            return self.space_mode
        raise ValueError(
            "name must be one of: 'cycle', 'phase', or 'space'."
        )


def analyze_tensor_low_rank(
    values: FloatArray,
    name: str = "field",
) -> TensorSvdDiagnostics:
    """Compute mode-wise SVD diagnostics for one snapshot field."""
    tensor = flatten_cycle_phase_field(
        values=values,
        name=name,
    )

    return TensorSvdDiagnostics(
        tensor_shape=(
            int(tensor.shape[0]),
            int(tensor.shape[1]),
            int(tensor.shape[2]),
        ),
        cycle_mode=compute_svd_spectrum(
            mode_unfolding(tensor=tensor, mode="cycle")
        ),
        phase_mode=compute_svd_spectrum(
            mode_unfolding(tensor=tensor, mode="phase")
        ),
        space_mode=compute_svd_spectrum(
            mode_unfolding(tensor=tensor, mode="space")
        ),
    )


@dataclass(frozen=True)
class FieldLowRankDiagnostics:
    """Raw and cycle-increment diagnostics for one physical field."""

    raw: TensorSvdDiagnostics
    cycle_increment: TensorSvdDiagnostics


def analyze_field_low_rank(
    values: FloatArray,
    name: str = "field",
) -> FieldLowRankDiagnostics:
    """Analyze both raw and cycle-start-referenced versions of one field."""
    raw_tensor = flatten_cycle_phase_field(
        values=values,
        name=name,
    )
    increment_tensor = cycle_increment_field(
        values=values,
        name=name,
    )

    return FieldLowRankDiagnostics(
        raw=analyze_tensor_low_rank(
            values=raw_tensor,
            name=name + " raw",
        ),
        cycle_increment=analyze_tensor_low_rank(
            values=increment_tensor,
            name=name + " cycle increment",
        ),
    )


@dataclass(frozen=True)
class TowerLowRankDiagnostics:
    """
    Low-rank diagnostics for the four LATIN-PGD target tower fields.

    The analyzed fields are:
        u       nodal displacement vector
        sigma   fiber stress
        eps_p   fiber plastic strain
        D       fiber damage
    """

    nodal_displacements: FieldLowRankDiagnostics
    fiber_stresses: FieldLowRankDiagnostics
    fiber_plastic_strains: FieldLowRankDiagnostics
    fiber_damages: FieldLowRankDiagnostics

    def as_dict(self) -> Dict[str, FieldLowRankDiagnostics]:
        """Return the target fields under compact physical names."""
        return {
            "u": self.nodal_displacements,
            "sigma": self.fiber_stresses,
            "eps_p": self.fiber_plastic_strains,
            "D": self.fiber_damages,
        }


def analyze_tower_low_rank(
    snapshots: TowerCyclePhaseSnapshots,
) -> TowerLowRankDiagnostics:
    """Analyze u, sigma, eps_p, and D from tower cycle-phase snapshots."""
    if not isinstance(snapshots, TowerCyclePhaseSnapshots):
        raise TypeError(
            "snapshots must be a TowerCyclePhaseSnapshots object."
        )

    return TowerLowRankDiagnostics(
        nodal_displacements=analyze_field_low_rank(
            values=snapshots.nodal_displacements,
            name="nodal_displacements",
        ),
        fiber_stresses=analyze_field_low_rank(
            values=snapshots.fiber_stresses,
            name="fiber_stresses",
        ),
        fiber_plastic_strains=analyze_field_low_rank(
            values=snapshots.fiber_plastic_strains,
            name="fiber_plastic_strains",
        ),
        fiber_damages=analyze_field_low_rank(
            values=snapshots.fiber_damages,
            name="fiber_damages",
        ),
    )

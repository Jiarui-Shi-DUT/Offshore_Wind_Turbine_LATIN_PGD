# -*- coding: utf-8 -*-
"""
Separated space-time representation used by the one-dimensional LATIN-PGD
global stage.

For each PGD mode j, the plastic part of the global correction is stored as

    Delta eps_p_dot_j(x, t)
        = lambda_dot_j(t) * eps_p_bar_j(x)

    Delta sigma_prime_j(x, t)
        = lambda_j(t) * sigma_bar_j(x)

where sigma_bar_j is obtained from the spatial equilibrium operator applied
to eps_p_bar_j.

This module only defines and reconstructs the reduced representation.  It does
not yet calculate new modes or update their temporal functions.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


def _as_finite_vector(
    name: str,
    values: FloatArray,
    expected_size: Optional[int] = None,
) -> FloatArray:
    """Convert an input to a finite one-dimensional float64 array."""
    array = np.asarray(values, dtype=np.float64)

    if array.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional.")
    if array.size < 1:
        raise ValueError(f"{name} must not be empty.")
    if expected_size is not None and array.size != expected_size:
        raise ValueError(
            f"{name} has size {array.size}; expected {expected_size}."
        )
    if not np.all(np.isfinite(array)):
        raise ValueError(f"{name} contains non-finite values.")

    return array.copy()


@dataclass
class PGDMode1D:
    """
    One separated space-time mode of the plastic LATIN correction.

    The spatial plastic-strain and stress fields are elementwise quantities.
    The temporal amplitude and its rate are stored at the LATIN time points.
    """

    spatial_plastic_strain: FloatArray
    spatial_stress: FloatArray
    temporal_amplitude: FloatArray
    temporal_rate: FloatArray
    iteration_added: int = 0

    def __post_init__(self) -> None:
        self.spatial_plastic_strain = _as_finite_vector(
            "spatial_plastic_strain",
            self.spatial_plastic_strain,
        )
        self.spatial_stress = _as_finite_vector(
            "spatial_stress",
            self.spatial_stress,
            expected_size=self.spatial_plastic_strain.size,
        )
        self.temporal_amplitude = _as_finite_vector(
            "temporal_amplitude",
            self.temporal_amplitude,
        )
        self.temporal_rate = _as_finite_vector(
            "temporal_rate",
            self.temporal_rate,
            expected_size=self.temporal_amplitude.size,
        )

        if self.iteration_added < 0:
            raise ValueError("iteration_added must be non-negative.")

    @property
    def n_elements(self) -> int:
        """Number of spatial finite elements represented by the mode."""
        return int(self.spatial_plastic_strain.size)

    @property
    def n_time(self) -> int:
        """Number of time points represented by the mode."""
        return int(self.temporal_amplitude.size)

    def plastic_strain_correction(self) -> FloatArray:
        """Reconstruct lambda(t) * eps_p_bar(x)."""
        return np.outer(
            self.temporal_amplitude,
            self.spatial_plastic_strain,
        )

    def plastic_strain_rate_correction(self) -> FloatArray:
        """Reconstruct lambda_dot(t) * eps_p_bar(x)."""
        return np.outer(
            self.temporal_rate,
            self.spatial_plastic_strain,
        )

    def stress_correction(self) -> FloatArray:
        """Reconstruct lambda(t) * sigma_bar(x)."""
        return np.outer(
            self.temporal_amplitude,
            self.spatial_stress,
        )

    def copy(self) -> "PGDMode1D":
        """Return a deep copy of the separated mode."""
        return PGDMode1D(
            spatial_plastic_strain=self.spatial_plastic_strain.copy(),
            spatial_stress=self.spatial_stress.copy(),
            temporal_amplitude=self.temporal_amplitude.copy(),
            temporal_rate=self.temporal_rate.copy(),
            iteration_added=self.iteration_added,
        )


@dataclass
class PGDBasis1D:
    """
    Collection of separated PGD modes sharing one space-time discretisation.
    """

    n_elements: int
    n_time: int
    modes: List[PGDMode1D] = field(default_factory=list)

    def __post_init__(self) -> None:
        if self.n_elements < 1:
            raise ValueError("n_elements must be at least one.")
        if self.n_time < 2:
            raise ValueError("n_time must be at least two.")

        initial_modes = list(self.modes)
        self.modes = []

        for mode in initial_modes:
            self.append(mode)

    @property
    def n_modes(self) -> int:
        """Current reduced-basis dimension."""
        return len(self.modes)

    @property
    def field_shape(self) -> Tuple[int, int]:
        """Shape of reconstructed time-by-element fields."""
        return self.n_time, self.n_elements

    def append(self, mode: PGDMode1D) -> None:
        """Append one compatible mode to the reduced basis."""
        if mode.n_elements != self.n_elements:
            raise ValueError(
                "The mode and basis must use the same number of elements."
            )
        if mode.n_time != self.n_time:
            raise ValueError(
                "The mode and basis must use the same number of time points."
            )

        self.modes.append(mode.copy())

    def clear(self) -> None:
        """Remove every PGD mode while preserving the discretisation."""
        self.modes.clear()

    def copy(self) -> "PGDBasis1D":
        """Return a deep copy of the complete reduced basis."""
        return PGDBasis1D(
            n_elements=self.n_elements,
            n_time=self.n_time,
            modes=[mode.copy() for mode in self.modes],
        )

    def spatial_plastic_strain_matrix(self) -> FloatArray:
        """
        Return the spatial basis matrix with one mode per column.

        Shape: (n_elements, n_modes).
        """
        if not self.modes:
            return np.zeros(
                (self.n_elements, 0),
                dtype=np.float64,
            )

        return np.column_stack(
            [
                mode.spatial_plastic_strain
                for mode in self.modes
            ]
        )

    def spatial_stress_matrix(self) -> FloatArray:
        """
        Return the spatial stress matrix with one mode per column.

        Shape: (n_elements, n_modes).
        """
        if not self.modes:
            return np.zeros(
                (self.n_elements, 0),
                dtype=np.float64,
            )

        return np.column_stack(
            [
                mode.spatial_stress
                for mode in self.modes
            ]
        )

    def temporal_amplitude_matrix(self) -> FloatArray:
        """
        Return the temporal amplitudes with one mode per column.

        Shape: (n_time, n_modes).
        """
        if not self.modes:
            return np.zeros(
                (self.n_time, 0),
                dtype=np.float64,
            )

        return np.column_stack(
            [
                mode.temporal_amplitude
                for mode in self.modes
            ]
        )

    def temporal_rate_matrix(self) -> FloatArray:
        """
        Return the temporal rates with one mode per column.

        Shape: (n_time, n_modes).
        """
        if not self.modes:
            return np.zeros(
                (self.n_time, 0),
                dtype=np.float64,
            )

        return np.column_stack(
            [
                mode.temporal_rate
                for mode in self.modes
            ]
        )

    def plastic_strain_correction(self) -> FloatArray:
        """
        Reconstruct the total PGD plastic-strain correction.

        Sum_j lambda_j(t) * eps_p_bar_j(x).
        """
        if not self.modes:
            return np.zeros(self.field_shape, dtype=np.float64)

        return (
            self.temporal_amplitude_matrix()
            @ self.spatial_plastic_strain_matrix().T
        )

    def plastic_strain_rate_correction(self) -> FloatArray:
        """
        Reconstruct the total PGD plastic-strain-rate correction.

        Sum_j lambda_dot_j(t) * eps_p_bar_j(x).
        """
        if not self.modes:
            return np.zeros(self.field_shape, dtype=np.float64)

        return (
            self.temporal_rate_matrix()
            @ self.spatial_plastic_strain_matrix().T
        )

    def stress_correction(self) -> FloatArray:
        """
        Reconstruct the total plastic stress correction.

        Sum_j lambda_j(t) * sigma_bar_j(x).
        """
        if not self.modes:
            return np.zeros(self.field_shape, dtype=np.float64)

        return (
            self.temporal_amplitude_matrix()
            @ self.spatial_stress_matrix().T
        )

# ---------------------------------------------------------------------------
# Tower LATIN-PGD value-style basis
# ---------------------------------------------------------------------------

def _tower_readonly_finite_vector(
    name: str,
    values: FloatArray,
    expected_size: Optional[int] = None,
) -> FloatArray:
    """Return a detached finite read-only float64 vector."""
    array = np.array(values, dtype=np.float64, copy=True)

    if array.ndim != 1:
        raise ValueError(name + " must be one-dimensional.")
    if array.size < 1:
        raise ValueError(name + " must not be empty.")
    if expected_size is not None and array.size != expected_size:
        raise ValueError(
            name
            + " has size "
            + str(array.size)
            + "; expected "
            + str(expected_size)
            + "."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " contains non-finite values.")

    array.setflags(write=False)
    return array


def _tower_nonnegative_integer(value: int, name: str) -> int:
    """Return a validated non-negative integer."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")

    result = int(value)
    if result < 0:
        raise ValueError(name + " must be non-negative.")
    return result


def _tower_positive_integer(value: int, name: str) -> int:
    """Return a validated strictly positive integer."""
    result = _tower_nonnegative_integer(value, name)
    if result < 1:
        raise ValueError(name + " must be at least one.")
    return result


@dataclass(frozen=True)
class PGDModeTower:
    """
    One immutable separated tower PGD mode on the material-point space.

    The mode stores

        p_j(q), s_j(q), lambda_j(t), lambda_dot_j(t),

    where p_j is the plastic-strain spatial mode and s_j is its associated
    equilibrated reference stress mode.
    """

    spatial_plastic_strain: FloatArray
    spatial_stress: FloatArray
    temporal_amplitude: FloatArray
    temporal_rate: FloatArray
    iteration_added: int = 0

    def __post_init__(self) -> None:
        spatial_plastic = _tower_readonly_finite_vector(
            "spatial_plastic_strain",
            self.spatial_plastic_strain,
        )
        spatial_stress = _tower_readonly_finite_vector(
            "spatial_stress",
            self.spatial_stress,
            expected_size=spatial_plastic.size,
        )
        temporal_amplitude = _tower_readonly_finite_vector(
            "temporal_amplitude",
            self.temporal_amplitude,
        )
        temporal_rate = _tower_readonly_finite_vector(
            "temporal_rate",
            self.temporal_rate,
            expected_size=temporal_amplitude.size,
        )
        iteration_added = _tower_nonnegative_integer(
            self.iteration_added,
            "iteration_added",
        )

        object.__setattr__(
            self,
            "spatial_plastic_strain",
            spatial_plastic,
        )
        object.__setattr__(
            self,
            "spatial_stress",
            spatial_stress,
        )
        object.__setattr__(
            self,
            "temporal_amplitude",
            temporal_amplitude,
        )
        object.__setattr__(
            self,
            "temporal_rate",
            temporal_rate,
        )
        object.__setattr__(
            self,
            "iteration_added",
            iteration_added,
        )

    @property
    def n_material_points(self) -> int:
        """Return the number of canonical tower material points."""
        return int(self.spatial_plastic_strain.size)

    @property
    def n_time(self) -> int:
        """Return the number of LATIN time points."""
        return int(self.temporal_amplitude.size)

    def plastic_strain_correction(self) -> FloatArray:
        """Return lambda(t) p(q)."""
        return np.outer(
            self.temporal_amplitude,
            self.spatial_plastic_strain,
        )

    def plastic_strain_rate_correction(self) -> FloatArray:
        """Return lambda_dot(t) p(q)."""
        return np.outer(
            self.temporal_rate,
            self.spatial_plastic_strain,
        )

    def stress_correction(self) -> FloatArray:
        """Return lambda(t) s(q)."""
        return np.outer(
            self.temporal_amplitude,
            self.spatial_stress,
        )

    def copy(self) -> "PGDModeTower":
        """Return a fully detached mode value."""
        return PGDModeTower(
            spatial_plastic_strain=self.spatial_plastic_strain,
            spatial_stress=self.spatial_stress,
            temporal_amplitude=self.temporal_amplitude,
            temporal_rate=self.temporal_rate,
            iteration_added=self.iteration_added,
        )

    def with_temporal_coordinates(
        self,
        temporal_amplitude: FloatArray,
        temporal_rate: FloatArray,
    ) -> "PGDModeTower":
        """Return the same spatial pair with new temporal coordinates."""
        return PGDModeTower(
            spatial_plastic_strain=self.spatial_plastic_strain,
            spatial_stress=self.spatial_stress,
            temporal_amplitude=temporal_amplitude,
            temporal_rate=temporal_rate,
            iteration_added=self.iteration_added,
        )


@dataclass(frozen=True)
class PGDBasisTower:
    """
    Immutable value-style tower PGD basis.

    No public in-place append or temporal-coordinate mutation is provided.
    Basis-changing operations return a new PGDBasisTower value.
    """

    n_material_points: int
    n_time: int
    modes: Tuple[PGDModeTower, ...] = ()

    def __post_init__(self) -> None:
        n_material_points = _tower_positive_integer(
            self.n_material_points,
            "n_material_points",
        )
        n_time = _tower_positive_integer(
            self.n_time,
            "n_time",
        )
        if n_time < 2:
            raise ValueError("n_time must be at least two.")

        copied_modes = []
        for mode in tuple(self.modes):
            if not isinstance(mode, PGDModeTower):
                raise TypeError(
                    "modes must contain only PGDModeTower objects."
                )
            if mode.n_material_points != n_material_points:
                raise ValueError(
                    "Every mode must use the basis material-point count."
                )
            if mode.n_time != n_time:
                raise ValueError(
                    "Every mode must use the basis time-point count."
                )
            copied_modes.append(mode.copy())

        object.__setattr__(
            self,
            "n_material_points",
            n_material_points,
        )
        object.__setattr__(
            self,
            "n_time",
            n_time,
        )
        object.__setattr__(
            self,
            "modes",
            tuple(copied_modes),
        )

    @property
    def n_modes(self) -> int:
        """Return the current reduced-basis dimension."""
        return int(len(self.modes))

    @property
    def field_shape(self) -> Tuple[int, int]:
        """Return (n_time, n_material_points)."""
        return self.n_time, self.n_material_points

    def copy(self) -> "PGDBasisTower":
        """Return a fully detached basis value."""
        return PGDBasisTower(
            n_material_points=self.n_material_points,
            n_time=self.n_time,
            modes=self.modes,
        )

    def with_appended(
        self,
        mode: PGDModeTower,
    ) -> "PGDBasisTower":
        """Return a new basis containing one additional compatible mode."""
        if not isinstance(mode, PGDModeTower):
            raise TypeError("mode must be a PGDModeTower.")
        if mode.n_material_points != self.n_material_points:
            raise ValueError(
                "mode and basis must use the same material-point count."
            )
        if mode.n_time != self.n_time:
            raise ValueError(
                "mode and basis must use the same time-point count."
            )

        return PGDBasisTower(
            n_material_points=self.n_material_points,
            n_time=self.n_time,
            modes=self.modes + (mode,),
        )

    def with_temporal_coordinates(
        self,
        temporal_amplitudes: FloatArray,
        temporal_rates: FloatArray,
    ) -> "PGDBasisTower":
        """
        Return the same spatial basis with jointly replaced temporal functions.

        Both matrices have shape (n_time, n_modes).
        """
        amplitudes = np.asarray(
            temporal_amplitudes,
            dtype=np.float64,
        )
        rates = np.asarray(
            temporal_rates,
            dtype=np.float64,
        )
        expected_shape = (
            self.n_time,
            self.n_modes,
        )

        if amplitudes.shape != expected_shape:
            raise ValueError(
                "temporal_amplitudes must have shape "
                + str(expected_shape)
                + "."
            )
        if rates.shape != expected_shape:
            raise ValueError(
                "temporal_rates must have shape "
                + str(expected_shape)
                + "."
            )
        if np.any(~np.isfinite(amplitudes)):
            raise ValueError(
                "temporal_amplitudes contain non-finite values."
            )
        if np.any(~np.isfinite(rates)):
            raise ValueError(
                "temporal_rates contain non-finite values."
            )

        new_modes = tuple(
            mode.with_temporal_coordinates(
                temporal_amplitude=amplitudes[:, index],
                temporal_rate=rates[:, index],
            )
            for index, mode in enumerate(self.modes)
        )

        return PGDBasisTower(
            n_material_points=self.n_material_points,
            n_time=self.n_time,
            modes=new_modes,
        )

    def spatial_plastic_strain_matrix(self) -> FloatArray:
        """Return P with shape (n_material_points, n_modes)."""
        if not self.modes:
            return np.zeros(
                (self.n_material_points, 0),
                dtype=np.float64,
            )
        return np.column_stack(
            [
                mode.spatial_plastic_strain
                for mode in self.modes
            ]
        )

    def spatial_stress_matrix(self) -> FloatArray:
        """Return S with shape (n_material_points, n_modes)."""
        if not self.modes:
            return np.zeros(
                (self.n_material_points, 0),
                dtype=np.float64,
            )
        return np.column_stack(
            [
                mode.spatial_stress
                for mode in self.modes
            ]
        )

    def temporal_amplitude_matrix(self) -> FloatArray:
        """Return Lambda with shape (n_time, n_modes)."""
        if not self.modes:
            return np.zeros(
                (self.n_time, 0),
                dtype=np.float64,
            )
        return np.column_stack(
            [
                mode.temporal_amplitude
                for mode in self.modes
            ]
        )

    def temporal_rate_matrix(self) -> FloatArray:
        """Return Lambda_dot with shape (n_time, n_modes)."""
        if not self.modes:
            return np.zeros(
                (self.n_time, 0),
                dtype=np.float64,
            )
        return np.column_stack(
            [
                mode.temporal_rate
                for mode in self.modes
            ]
        )

    def plastic_strain_correction(self) -> FloatArray:
        """Return Lambda P^T on the complete time-material grid."""
        if not self.modes:
            return np.zeros(
                self.field_shape,
                dtype=np.float64,
            )
        return (
            self.temporal_amplitude_matrix()
            @ self.spatial_plastic_strain_matrix().T
        )

    def plastic_strain_rate_correction(self) -> FloatArray:
        """Return Lambda_dot P^T on the complete time-material grid."""
        if not self.modes:
            return np.zeros(
                self.field_shape,
                dtype=np.float64,
            )
        return (
            self.temporal_rate_matrix()
            @ self.spatial_plastic_strain_matrix().T
        )

    def stress_correction(self) -> FloatArray:
        """Return Lambda S^T on the complete time-material grid."""
        if not self.modes:
            return np.zeros(
                self.field_shape,
                dtype=np.float64,
            )
        return (
            self.temporal_amplitude_matrix()
            @ self.spatial_stress_matrix().T
        )

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

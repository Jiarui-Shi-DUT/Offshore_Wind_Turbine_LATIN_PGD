# -*- coding: utf-8 -*-
"""
Per-cycle diagnostics for fully reversed nonlinear tower responses.

This module converts one multi-cycle ``NonlinearReversedResponse`` into
cycle-by-cycle scalar diagnostics. The extraction is deliberately separated
from the nonlinear solver so that the same full-order response can later be
used for:

    - cyclic stabilization assessment;
    - cyclic hardening or softening assessment;
    - residual-displacement and ratcheting assessment;
    - damage-increment tracking;
    - force-displacement work tracking;
    - LATIN-PGD full-order reference comparisons.

One cycle is indexed by the five exact loading checkpoints

    0 -> +F_a -> 0 -> -F_a -> 0.

The force-displacement line integral is stored as external work. Before a
closed stabilized loop is reached, it must not automatically be interpreted
as pure material dissipation because stored internal energy can also change
between the beginning and end of a cycle.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from examples.nonlinear_tower_reversed_response import (
    NonlinearReversedResponse,
)
from fem.tower_loading import ReversedTopForceHistory


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _finite_scalar(value: float, name: str) -> float:
    """Return a validated finite scalar."""
    result = float(value)
    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _nonnegative_scalar(value: float, name: str) -> float:
    """Return a validated non-negative scalar."""
    result = _finite_scalar(value, name)
    if result < 0.0:
        raise ValueError(name + " must be non-negative.")
    return result


def _validated_cycle_number(
    cycle_number: int,
    n_cycles: int,
) -> int:
    """Return a validated one-based cycle number."""
    if isinstance(cycle_number, bool):
        raise TypeError("cycle_number must be an integer.")
    try:
        result = int(cycle_number)
    except (TypeError, ValueError) as error:
        raise TypeError(
            "cycle_number must be an integer."
        ) from error

    if result != cycle_number:
        raise TypeError("cycle_number must be an integer.")
    if result < 1 or result > n_cycles:
        raise ValueError(
            "cycle_number must lie between 1 and loading.n_cycles."
        )
    return result


@dataclass(frozen=True)
class CycleIndices:
    """Exact history indices for one complete reversed cycle."""

    cycle_number: int
    start: int
    positive_peak: int
    first_zero: int
    negative_peak: int
    end: int

    def __post_init__(self) -> None:
        values = (
            self.cycle_number,
            self.start,
            self.positive_peak,
            self.first_zero,
            self.negative_peak,
            self.end,
        )
        if any(isinstance(value, bool) for value in values):
            raise TypeError("Cycle indices must be integers.")

        converted = tuple(int(value) for value in values)
        if converted != values:
            raise TypeError("Cycle indices must be integers.")

        (
            cycle_number,
            start,
            positive_peak,
            first_zero,
            negative_peak,
            end,
        ) = converted

        if cycle_number < 1:
            raise ValueError("cycle_number must be at least one.")
        if start < 0:
            raise ValueError("start must be non-negative.")
        if not (
            start
            < positive_peak
            < first_zero
            < negative_peak
            < end
        ):
            raise ValueError(
                "Cycle indices must be strictly ordered."
            )

        object.__setattr__(
            self,
            "cycle_number",
            cycle_number,
        )
        object.__setattr__(self, "start", start)
        object.__setattr__(
            self,
            "positive_peak",
            positive_peak,
        )
        object.__setattr__(
            self,
            "first_zero",
            first_zero,
        )
        object.__setattr__(
            self,
            "negative_peak",
            negative_peak,
        )
        object.__setattr__(self, "end", end)

    @property
    def slice(self) -> slice:
        """Return an inclusive Python slice for the cycle."""
        return slice(self.start, self.end + 1)


@dataclass(frozen=True)
class CycleDiagnostics:
    """Scalar diagnostics extracted from one complete cycle."""

    cycle_number: int
    indices: CycleIndices

    displacement_at_positive_peak: float
    displacement_at_negative_peak: float
    maximum_displacement: float
    minimum_displacement: float
    residual_displacement: float
    displacement_range: float

    critical_stress_at_positive_peak: float
    critical_stress_at_negative_peak: float
    maximum_critical_stress: float
    minimum_critical_stress: float
    critical_stress_range: float

    critical_plastic_strain_at_positive_peak: float
    critical_plastic_strain_at_negative_peak: float
    critical_plastic_strain_at_end: float
    critical_plastic_strain_increment: float
    critical_plastic_strain_range: float

    critical_backstress_at_positive_peak: float
    critical_backstress_at_negative_peak: float
    maximum_critical_backstress: float
    minimum_critical_backstress: float
    critical_backstress_at_end: float
    critical_backstress_range: float

    critical_r_bar_at_end: float
    critical_r_bar_increment: float

    maximum_damage_at_start: float
    maximum_damage_at_end: float
    maximum_damage_increment: float
    critical_damage_at_start: float
    critical_damage_at_end: float
    critical_damage_increment: float

    signed_external_work: float
    external_work_magnitude: float

    maximum_newton_iterations: int
    maximum_residual_norm: float

    def __post_init__(self) -> None:
        if isinstance(self.cycle_number, bool):
            raise TypeError("cycle_number must be an integer.")
        cycle_number = int(self.cycle_number)
        if cycle_number != self.cycle_number:
            raise TypeError("cycle_number must be an integer.")
        if cycle_number < 1:
            raise ValueError("cycle_number must be at least one.")
        if not isinstance(self.indices, CycleIndices):
            raise TypeError("indices must be a CycleIndices object.")
        if self.indices.cycle_number != cycle_number:
            raise ValueError(
                "indices.cycle_number must match cycle_number."
            )

        scalar_names = (
            "displacement_at_positive_peak",
            "displacement_at_negative_peak",
            "maximum_displacement",
            "minimum_displacement",
            "residual_displacement",
            "displacement_range",
            "critical_stress_at_positive_peak",
            "critical_stress_at_negative_peak",
            "maximum_critical_stress",
            "minimum_critical_stress",
            "critical_stress_range",
            "critical_plastic_strain_at_positive_peak",
            "critical_plastic_strain_at_negative_peak",
            "critical_plastic_strain_at_end",
            "critical_plastic_strain_increment",
            "critical_plastic_strain_range",
            "critical_backstress_at_positive_peak",
            "critical_backstress_at_negative_peak",
            "maximum_critical_backstress",
            "minimum_critical_backstress",
            "critical_backstress_at_end",
            "critical_backstress_range",
            "critical_r_bar_at_end",
            "critical_r_bar_increment",
            "maximum_damage_at_start",
            "maximum_damage_at_end",
            "maximum_damage_increment",
            "critical_damage_at_start",
            "critical_damage_at_end",
            "critical_damage_increment",
            "signed_external_work",
            "external_work_magnitude",
            "maximum_residual_norm",
        )
        nonnegative_names = (
            "displacement_range",
            "critical_stress_range",
            "critical_plastic_strain_range",
            "critical_backstress_range",
            "maximum_damage_at_start",
            "maximum_damage_at_end",
            "maximum_damage_increment",
            "critical_damage_at_start",
            "critical_damage_at_end",
            "critical_damage_increment",
            "external_work_magnitude",
            "maximum_residual_norm",
        )

        for name in scalar_names:
            value = _finite_scalar(
                getattr(self, name),
                name,
            )
            object.__setattr__(self, name, value)

        for name in nonnegative_names:
            _nonnegative_scalar(
                getattr(self, name),
                name,
            )

        if self.maximum_displacement < self.minimum_displacement:
            raise ValueError(
                "maximum_displacement must not be smaller "
                "than minimum_displacement."
            )
        if (
            self.maximum_critical_stress
            < self.minimum_critical_stress
        ):
            raise ValueError(
                "maximum_critical_stress must not be smaller "
                "than minimum_critical_stress."
            )
        if (
            self.maximum_critical_backstress
            < self.minimum_critical_backstress
        ):
            raise ValueError(
                "maximum_critical_backstress must not be smaller "
                "than minimum_critical_backstress."
            )
        if (
            self.maximum_damage_at_end
            + 1.0e-14
            < self.maximum_damage_at_start
        ):
            raise ValueError(
                "Maximum damage must be non-decreasing."
            )
        if (
            self.critical_damage_at_end
            + 1.0e-14
            < self.critical_damage_at_start
        ):
            raise ValueError(
                "Critical-fiber damage must be non-decreasing."
            )

        expected_displacement_range = (
            self.maximum_displacement
            - self.minimum_displacement
        )
        if not np.isclose(
            self.displacement_range,
            expected_displacement_range,
            rtol=1.0e-12,
            atol=1.0e-14,
        ):
            raise ValueError(
                "displacement_range is inconsistent with extrema."
            )

        expected_stress_range = (
            self.maximum_critical_stress
            - self.minimum_critical_stress
        )
        if not np.isclose(
            self.critical_stress_range,
            expected_stress_range,
            rtol=1.0e-12,
            atol=1.0e-12,
        ):
            raise ValueError(
                "critical_stress_range is inconsistent with extrema."
            )

        expected_backstress_range = (
            self.maximum_critical_backstress
            - self.minimum_critical_backstress
        )
        if not np.isclose(
            self.critical_backstress_range,
            expected_backstress_range,
            rtol=1.0e-12,
            atol=1.0e-12,
        ):
            raise ValueError(
                "critical_backstress_range is inconsistent with extrema."
            )

        if isinstance(self.maximum_newton_iterations, bool):
            raise TypeError(
                "maximum_newton_iterations must be an integer."
            )
        maximum_newton_iterations = int(
            self.maximum_newton_iterations
        )
        if (
            maximum_newton_iterations
            != self.maximum_newton_iterations
        ):
            raise TypeError(
                "maximum_newton_iterations must be an integer."
            )
        if maximum_newton_iterations < 1:
            raise ValueError(
                "maximum_newton_iterations must be positive."
            )
        object.__setattr__(
            self,
            "maximum_newton_iterations",
            maximum_newton_iterations,
        )
        object.__setattr__(
            self,
            "cycle_number",
            cycle_number,
        )


@dataclass(frozen=True)
class MulticycleDiagnostics:
    """Ordered diagnostics for every complete loading cycle."""

    cycles: Tuple[CycleDiagnostics, ...]

    def __post_init__(self) -> None:
        cycles = tuple(self.cycles)
        if len(cycles) < 1:
            raise ValueError(
                "cycles must contain at least one cycle."
            )
        for expected_number, cycle in enumerate(
            cycles,
            start=1,
        ):
            if not isinstance(cycle, CycleDiagnostics):
                raise TypeError(
                    "Every item must be a CycleDiagnostics object."
                )
            if cycle.cycle_number != expected_number:
                raise ValueError(
                    "cycles must be ordered and consecutively numbered."
                )
        object.__setattr__(self, "cycles", cycles)

    @property
    def n_cycles(self) -> int:
        """Return the number of stored cycles."""
        return len(self.cycles)

    @property
    def cycle_numbers(self) -> IntArray:
        """Return one-based cycle numbers."""
        return np.asarray(
            [cycle.cycle_number for cycle in self.cycles],
            dtype=np.int64,
        )

    def _float_history(self, attribute: str) -> FloatArray:
        """Return one scalar diagnostic across all cycles."""
        return np.asarray(
            [
                float(getattr(cycle, attribute))
                for cycle in self.cycles
            ],
            dtype=np.float64,
        )

    @property
    def residual_displacements(self) -> FloatArray:
        """Return the zero-force end displacement of every cycle."""
        return self._float_history("residual_displacement")

    @property
    def displacement_ranges(self) -> FloatArray:
        """Return the displacement range of every cycle."""
        return self._float_history("displacement_range")

    @property
    def critical_stress_ranges(self) -> FloatArray:
        """Return the fixed-fiber stress range of every cycle."""
        return self._float_history("critical_stress_range")

    @property
    def critical_plastic_strain_ends(self) -> FloatArray:
        """Return fixed-fiber plastic strain at every cycle end."""
        return self._float_history(
            "critical_plastic_strain_at_end"
        )

    @property
    def critical_plastic_strain_increments(self) -> FloatArray:
        """Return fixed-fiber plastic-strain increment per cycle."""
        return self._float_history(
            "critical_plastic_strain_increment"
        )

    @property
    def critical_backstress_ends(self) -> FloatArray:
        """Return fixed-fiber backstress at every cycle end."""
        return self._float_history(
            "critical_backstress_at_end"
        )

    @property
    def critical_r_bar_ends(self) -> FloatArray:
        """Return fixed-fiber isotropic variable at every cycle end."""
        return self._float_history("critical_r_bar_at_end")

    @property
    def maximum_damage_ends(self) -> FloatArray:
        """Return global maximum damage at every cycle end."""
        return self._float_history("maximum_damage_at_end")

    @property
    def maximum_damage_increments(self) -> FloatArray:
        """Return global maximum-damage increment per cycle."""
        return self._float_history(
            "maximum_damage_increment"
        )

    @property
    def critical_damage_ends(self) -> FloatArray:
        """Return fixed-fiber damage at every cycle end."""
        return self._float_history("critical_damage_at_end")

    @property
    def external_work_magnitudes(self) -> FloatArray:
        """Return force-displacement work magnitude per cycle."""
        return self._float_history("external_work_magnitude")


def cycle_indices(
    loading: ReversedTopForceHistory,
    cycle_number: int,
) -> CycleIndices:
    """Return exact checkpoint indices for one one-based cycle."""
    if not isinstance(loading, ReversedTopForceHistory):
        raise TypeError(
            "loading must be a ReversedTopForceHistory."
        )
    cycle_number = _validated_cycle_number(
        cycle_number,
        loading.n_cycles,
    )

    increments = loading.increments_per_cycle
    quarter = increments // 4
    start = (cycle_number - 1) * increments

    return CycleIndices(
        cycle_number=cycle_number,
        start=start,
        positive_peak=start + quarter,
        first_zero=start + 2 * quarter,
        negative_peak=start + 3 * quarter,
        end=start + increments,
    )


def signed_force_displacement_work(
    forces: FloatArray,
    displacements: FloatArray,
) -> float:
    """
    Return the discrete line integral of force with respect to displacement.

    The result has units of joules when force is in newtons and displacement
    is in metres. It is external work over the supplied path. Only after the
    response approaches a closed stabilized loop can it be interpreted as
    approximately equal to hysteretic dissipation.
    """
    force_values = np.asarray(
        forces,
        dtype=np.float64,
    )
    displacement_values = np.asarray(
        displacements,
        dtype=np.float64,
    )

    if force_values.ndim != 1:
        raise ValueError("forces must be one-dimensional.")
    if displacement_values.shape != force_values.shape:
        raise ValueError(
            "displacements must match forces."
        )
    if force_values.size < 2:
        raise ValueError(
            "At least two path points are required."
        )
    if np.any(~np.isfinite(force_values)):
        raise ValueError("forces must be finite.")
    if np.any(~np.isfinite(displacement_values)):
        raise ValueError("displacements must be finite.")

    displacement_increments = np.diff(
        displacement_values
    )
    average_forces = 0.5 * (
        force_values[:-1]
        + force_values[1:]
    )
    return float(
        np.sum(
            average_forces * displacement_increments
        )
    )


def extract_cycle_diagnostics(
    response: NonlinearReversedResponse,
    cycle_number: int,
) -> CycleDiagnostics:
    """Extract all scalar diagnostics for one complete cycle."""
    if not isinstance(
        response,
        NonlinearReversedResponse,
    ):
        raise TypeError(
            "response must be a NonlinearReversedResponse."
        )

    indices = cycle_indices(
        loading=response.loading,
        cycle_number=cycle_number,
    )
    cycle_slice = indices.slice

    displacements = response.top_displacements[
        cycle_slice
    ]
    forces = response.loading.forces[cycle_slice]
    critical_stresses = (
        response.critical_fiber_stresses[cycle_slice]
    )
    critical_plastic_strains = (
        response.critical_fiber_plastic_strains[
            cycle_slice
        ]
    )
    critical_backstresses = (
        response.critical_fiber_backstresses[
            cycle_slice
        ]
    )

    maximum_damages = response.maximum_damages
    critical_damages = (
        response.critical_fiber_damages
    )
    critical_r_bars = (
        response.critical_fiber_r_bars
    )

    signed_work = signed_force_displacement_work(
        forces=forces,
        displacements=displacements,
    )

    maximum_displacement = float(
        np.max(displacements)
    )
    minimum_displacement = float(
        np.min(displacements)
    )
    maximum_critical_stress = float(
        np.max(critical_stresses)
    )
    minimum_critical_stress = float(
        np.min(critical_stresses)
    )
    maximum_critical_backstress = float(
        np.max(critical_backstresses)
    )
    minimum_critical_backstress = float(
        np.min(critical_backstresses)
    )
    critical_plastic_range = float(
        np.max(critical_plastic_strains)
        - np.min(critical_plastic_strains)
    )

    maximum_damage_start = float(
        maximum_damages[indices.start]
    )
    maximum_damage_end = float(
        maximum_damages[indices.end]
    )
    critical_damage_start = float(
        critical_damages[indices.start]
    )
    critical_damage_end = float(
        critical_damages[indices.end]
    )

    return CycleDiagnostics(
        cycle_number=indices.cycle_number,
        indices=indices,
        displacement_at_positive_peak=float(
            response.top_displacements[
                indices.positive_peak
            ]
        ),
        displacement_at_negative_peak=float(
            response.top_displacements[
                indices.negative_peak
            ]
        ),
        maximum_displacement=maximum_displacement,
        minimum_displacement=minimum_displacement,
        residual_displacement=float(
            response.top_displacements[indices.end]
        ),
        displacement_range=float(
            maximum_displacement
            - minimum_displacement
        ),
        critical_stress_at_positive_peak=float(
            response.critical_fiber_stresses[
                indices.positive_peak
            ]
        ),
        critical_stress_at_negative_peak=float(
            response.critical_fiber_stresses[
                indices.negative_peak
            ]
        ),
        maximum_critical_stress=(
            maximum_critical_stress
        ),
        minimum_critical_stress=(
            minimum_critical_stress
        ),
        critical_stress_range=float(
            maximum_critical_stress
            - minimum_critical_stress
        ),
        critical_plastic_strain_at_positive_peak=float(
            response.critical_fiber_plastic_strains[
                indices.positive_peak
            ]
        ),
        critical_plastic_strain_at_negative_peak=float(
            response.critical_fiber_plastic_strains[
                indices.negative_peak
            ]
        ),
        critical_plastic_strain_at_end=float(
            response.critical_fiber_plastic_strains[
                indices.end
            ]
        ),
        critical_plastic_strain_increment=float(
            response.critical_fiber_plastic_strains[
                indices.end
            ]
            - response.critical_fiber_plastic_strains[
                indices.start
            ]
        ),
        critical_plastic_strain_range=(
            critical_plastic_range
        ),
        critical_backstress_at_positive_peak=float(
            response.critical_fiber_backstresses[
                indices.positive_peak
            ]
        ),
        critical_backstress_at_negative_peak=float(
            response.critical_fiber_backstresses[
                indices.negative_peak
            ]
        ),
        maximum_critical_backstress=(
            maximum_critical_backstress
        ),
        minimum_critical_backstress=(
            minimum_critical_backstress
        ),
        critical_backstress_at_end=float(
            response.critical_fiber_backstresses[
                indices.end
            ]
        ),
        critical_backstress_range=float(
            maximum_critical_backstress
            - minimum_critical_backstress
        ),
        critical_r_bar_at_end=float(
            critical_r_bars[indices.end]
        ),
        critical_r_bar_increment=float(
            critical_r_bars[indices.end]
            - critical_r_bars[indices.start]
        ),
        maximum_damage_at_start=(
            maximum_damage_start
        ),
        maximum_damage_at_end=maximum_damage_end,
        maximum_damage_increment=float(
            maximum_damage_end
            - maximum_damage_start
        ),
        critical_damage_at_start=(
            critical_damage_start
        ),
        critical_damage_at_end=critical_damage_end,
        critical_damage_increment=float(
            critical_damage_end
            - critical_damage_start
        ),
        signed_external_work=signed_work,
        external_work_magnitude=abs(signed_work),
        maximum_newton_iterations=int(
            np.max(response.iterations[cycle_slice])
        ),
        maximum_residual_norm=float(
            np.max(response.residual_norms[cycle_slice])
        ),
    )


def extract_multicycle_diagnostics(
    response: NonlinearReversedResponse,
) -> MulticycleDiagnostics:
    """Extract ordered diagnostics for all complete cycles."""
    if not isinstance(
        response,
        NonlinearReversedResponse,
    ):
        raise TypeError(
            "response must be a NonlinearReversedResponse."
        )

    cycles = tuple(
        extract_cycle_diagnostics(
            response=response,
            cycle_number=cycle_number,
        )
        for cycle_number in range(
            1,
            response.loading.n_cycles + 1,
        )
    )
    return MulticycleDiagnostics(cycles=cycles)

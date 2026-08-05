# -*- coding: utf-8 -*-
"""
Time histories for horizontal loading of the offshore-wind-turbine tower.

The initial fatigue-loading definition is a positive-mean pulsating force,

    F(t) = F_mean + F_amplitude * sin(2*pi*t/T),

with

    F_min / F_max = R_F.

For R_F = 0.1,

    F_mean = 0.55 * F_max,
    F_amplitude = 0.45 * F_max.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]


def _validated_finite_scalar(
    value: float,
    name: str,
) -> float:
    """Return a finite floating-point scalar."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise TypeError(name + " must be a real scalar.") from error

    if not np.isfinite(result):
        raise ValueError(name + " must be finite.")
    return result


def _validated_positive_scalar(
    value: float,
    name: str,
) -> float:
    """Return a finite positive floating-point scalar."""
    result = _validated_finite_scalar(value, name)
    if result <= 0.0:
        raise ValueError(name + " must be positive.")
    return result


def _validated_integer(
    value: int,
    name: str,
) -> int:
    """Return a validated integer."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


@dataclass(frozen=True)
class PulsatingTopForceHistory:
    """
    Discrete positive-mean sinusoidal tower-top force history.

    Parameters
    ----------
    maximum_force
        Maximum horizontal force F_max.
    force_ratio
        Ratio R_F = F_min / F_max. The pulsating case requires
        0 <= R_F <= 1.
    period
        Duration of one force cycle.
    n_cycles
        Number of complete cycles.
    increments_per_cycle
        Number of equal time increments in one cycle. It must be
        divisible by four so that the discrete history contains the
        exact maximum and minimum forces.
    times
        Time coordinates, including both the initial and final points.
    forces
        Horizontal-force values evaluated at ``times``.
    """

    maximum_force: float
    force_ratio: float
    period: float
    n_cycles: int
    increments_per_cycle: int
    times: FloatArray
    forces: FloatArray

    def __post_init__(self) -> None:
        maximum_force = _validated_positive_scalar(
            self.maximum_force,
            "maximum_force",
        )
        force_ratio = _validated_finite_scalar(
            self.force_ratio,
            "force_ratio",
        )
        period = _validated_positive_scalar(
            self.period,
            "period",
        )
        n_cycles = _validated_integer(
            self.n_cycles,
            "n_cycles",
        )
        increments_per_cycle = _validated_integer(
            self.increments_per_cycle,
            "increments_per_cycle",
        )

        if force_ratio < 0.0 or force_ratio > 1.0:
            raise ValueError(
                "force_ratio must lie in [0, 1] for a "
                "positive-mean pulsating force."
            )
        if n_cycles < 1:
            raise ValueError("n_cycles must be at least 1.")
        if increments_per_cycle < 4:
            raise ValueError(
                "increments_per_cycle must be at least 4."
            )
        if increments_per_cycle % 4 != 0:
            raise ValueError(
                "increments_per_cycle must be divisible by 4."
            )

        times = np.asarray(self.times, dtype=np.float64)
        forces = np.asarray(self.forces, dtype=np.float64)

        expected_size = (
            n_cycles * increments_per_cycle + 1
        )
        if times.shape != (expected_size,):
            raise ValueError(
                "times must contain "
                "n_cycles * increments_per_cycle + 1 values."
            )
        if forces.shape != times.shape:
            raise ValueError("forces must match times.")
        if np.any(~np.isfinite(times)):
            raise ValueError("times must be finite.")
        if np.any(~np.isfinite(forces)):
            raise ValueError("forces must be finite.")
        if np.any(np.diff(times) <= 0.0):
            raise ValueError(
                "times must be strictly increasing."
            )

        object.__setattr__(
            self,
            "maximum_force",
            maximum_force,
        )
        object.__setattr__(
            self,
            "force_ratio",
            force_ratio,
        )
        object.__setattr__(self, "period", period)
        object.__setattr__(self, "n_cycles", n_cycles)
        object.__setattr__(
            self,
            "increments_per_cycle",
            increments_per_cycle,
        )
        object.__setattr__(self, "times", times)
        object.__setattr__(self, "forces", forces)

    @property
    def minimum_force(self) -> float:
        """Return F_min = R_F * F_max."""
        return float(self.force_ratio * self.maximum_force)

    @property
    def mean_force(self) -> float:
        """Return the sinusoidal mean force."""
        return float(
            0.5 * (self.maximum_force + self.minimum_force)
        )

    @property
    def force_amplitude(self) -> float:
        """Return the sinusoidal force amplitude."""
        return float(
            0.5 * (self.maximum_force - self.minimum_force)
        )

    @property
    def angular_frequency(self) -> float:
        """Return omega = 2*pi/T."""
        return float(2.0 * np.pi / self.period)

    @property
    def time_increment(self) -> float:
        """Return the uniform time increment."""
        return float(self.period / self.increments_per_cycle)

    @property
    def n_time_points(self) -> int:
        """Return the number of stored time points."""
        return int(self.times.size)


def evaluate_pulsating_top_force(
    time: FloatArray,
    maximum_force: float,
    force_ratio: float = 0.1,
    period: float = 1.0,
) -> FloatArray:
    """
    Evaluate the positive-mean sinusoidal force at arbitrary times.

    The force definition is

        F(t) = F_mean + F_amplitude * sin(2*pi*t/T).
    """
    maximum_force = _validated_positive_scalar(
        maximum_force,
        "maximum_force",
    )
    force_ratio = _validated_finite_scalar(
        force_ratio,
        "force_ratio",
    )
    period = _validated_positive_scalar(period, "period")

    if force_ratio < 0.0 or force_ratio > 1.0:
        raise ValueError(
            "force_ratio must lie in [0, 1] for a "
            "positive-mean pulsating force."
        )

    times = np.asarray(time, dtype=np.float64)
    if np.any(~np.isfinite(times)):
        raise ValueError("time must be finite.")

    minimum_force = force_ratio * maximum_force
    mean_force = 0.5 * (maximum_force + minimum_force)
    force_amplitude = 0.5 * (
        maximum_force - minimum_force
    )
    angular_frequency = 2.0 * np.pi / period

    return np.asarray(
        mean_force
        + force_amplitude
        * np.sin(angular_frequency * times),
        dtype=np.float64,
    )


def create_pulsating_top_force_history(
    maximum_force: float,
    force_ratio: float = 0.1,
    period: float = 1.0,
    n_cycles: int = 1,
    increments_per_cycle: int = 100,
) -> PulsatingTopForceHistory:
    """
    Create a uniformly sampled pulsating tower-top force history.

    The time vector contains the initial point and the final point of the
    last complete cycle. Therefore,

        n_time_points = n_cycles * increments_per_cycle + 1.
    """
    maximum_force = _validated_positive_scalar(
        maximum_force,
        "maximum_force",
    )
    force_ratio = _validated_finite_scalar(
        force_ratio,
        "force_ratio",
    )
    period = _validated_positive_scalar(period, "period")
    n_cycles = _validated_integer(n_cycles, "n_cycles")
    increments_per_cycle = _validated_integer(
        increments_per_cycle,
        "increments_per_cycle",
    )

    if force_ratio < 0.0 or force_ratio > 1.0:
        raise ValueError(
            "force_ratio must lie in [0, 1] for a "
            "positive-mean pulsating force."
        )
    if n_cycles < 1:
        raise ValueError("n_cycles must be at least 1.")
    if increments_per_cycle < 4:
        raise ValueError(
            "increments_per_cycle must be at least 4."
        )
    if increments_per_cycle % 4 != 0:
        raise ValueError(
            "increments_per_cycle must be divisible by 4."
        )

    n_increments = n_cycles * increments_per_cycle
    times = np.linspace(
        0.0,
        n_cycles * period,
        n_increments + 1,
        dtype=np.float64,
    )
    forces = evaluate_pulsating_top_force(
        time=times,
        maximum_force=maximum_force,
        force_ratio=force_ratio,
        period=period,
    )

    return PulsatingTopForceHistory(
        maximum_force=maximum_force,
        force_ratio=force_ratio,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
        times=times,
        forces=forces,
    )

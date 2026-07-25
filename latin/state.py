# -*- coding: utf-8 -*-
"""
Space-time state representation for the one-dimensional LATIN solver.

The paper defines the LATIN state by the fields
    {eps_p_dot, eps_e, X_dot, D_dot, sigma, Z, Y}.
For the present one-dimensional mixed-hardening model:
    X = {alpha, r_bar}
    Z = {beta, R_bar}

The integrated internal variables are also stored explicitly because they are
needed by the nonlinear elastic-damage law and by numerical verification.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence, Tuple

import numpy as np
from numpy.typing import NDArray

from latin.initialization import ElasticInitialization
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


@dataclass
class LatinState:
    """Complete one-dimensional LATIN state on the time-space grid."""

    time: FloatArray

    # Primary rate and state fields used by the LATIN formulation
    plastic_strain_rate: FloatArray
    elastic_strain: FloatArray
    alpha_rate: FloatArray
    r_bar_rate: FloatArray
    damage_rate: FloatArray
    stress: FloatArray
    beta: FloatArray
    R_bar: FloatArray
    energy_release_rate: FloatArray

    # Time-integrated internal variables retained explicitly
    plastic_strain: FloatArray
    alpha: FloatArray
    r_bar: FloatArray
    damage: FloatArray

    def __post_init__(self) -> None:
        """Normalize arrays and enforce one common time-space shape."""
        self.time = np.asarray(self.time, dtype=np.float64)

        if self.time.ndim != 1:
            raise ValueError("time must be a one-dimensional array.")
        if self.time.size < 2:
            raise ValueError("At least two time points are required.")
        if np.any(np.diff(self.time) <= 0.0):
            raise ValueError("time must be strictly increasing.")
        if not np.all(np.isfinite(self.time)):
            raise ValueError("time contains non-finite values.")

        arrays = (
            ("plastic_strain_rate", self.plastic_strain_rate),
            ("elastic_strain", self.elastic_strain),
            ("alpha_rate", self.alpha_rate),
            ("r_bar_rate", self.r_bar_rate),
            ("damage_rate", self.damage_rate),
            ("stress", self.stress),
            ("beta", self.beta),
            ("R_bar", self.R_bar),
            ("energy_release_rate", self.energy_release_rate),
            ("plastic_strain", self.plastic_strain),
            ("alpha", self.alpha),
            ("r_bar", self.r_bar),
            ("damage", self.damage),
        )

        expected_shape = None

        for field_name, field_value in arrays:
            array = np.asarray(field_value, dtype=np.float64)

            if array.ndim != 2:
                raise ValueError(
                    f"{field_name} must be a two-dimensional "
                    "time-by-element array."
                )

            if expected_shape is None:
                expected_shape = array.shape
            elif array.shape != expected_shape:
                raise ValueError(
                    f"{field_name} has shape {array.shape}; "
                    f"expected {expected_shape}."
                )

            if array.shape[0] != self.time.size:
                raise ValueError(
                    f"{field_name} must contain one row per time point."
                )

            if not np.all(np.isfinite(array)):
                raise ValueError(
                    f"{field_name} contains non-finite values."
                )

            setattr(self, field_name, array)

        if expected_shape is None or expected_shape[1] < 1:
            raise ValueError("The LATIN state must contain at least one element.")

        if np.any(self.damage < 0.0) or np.any(self.damage >= 1.0):
            raise ValueError("damage must satisfy 0 <= damage < 1.")

    @property
    def n_time(self) -> int:
        """Number of time points."""
        return int(self.time.size)

    @property
    def n_elements(self) -> int:
        """Number of finite elements or material points."""
        return int(self.stress.shape[1])

    @property
    def field_shape(self) -> Tuple[int, int]:
        """Common shape (n_time, n_elements) of all space-time fields."""
        return self.stress.shape

    def copy(self) -> "LatinState":
        """Return a deep copy suitable for a new LATIN half-iteration."""
        return LatinState(
            time=self.time.copy(),
            plastic_strain_rate=self.plastic_strain_rate.copy(),
            elastic_strain=self.elastic_strain.copy(),
            alpha_rate=self.alpha_rate.copy(),
            r_bar_rate=self.r_bar_rate.copy(),
            damage_rate=self.damage_rate.copy(),
            stress=self.stress.copy(),
            beta=self.beta.copy(),
            R_bar=self.R_bar.copy(),
            energy_release_rate=self.energy_release_rate.copy(),
            plastic_strain=self.plastic_strain.copy(),
            alpha=self.alpha.copy(),
            r_bar=self.r_bar.copy(),
            damage=self.damage.copy(),
        )

    @classmethod
    def zeros(
        cls,
        time: FloatArray,
        n_elements: int,
    ) -> "LatinState":
        """Create a zero state on a specified time-space grid."""
        time_array = np.asarray(time, dtype=np.float64)

        if n_elements < 1:
            raise ValueError("n_elements must be at least one.")

        shape = (time_array.size, int(n_elements))

        def zero_field() -> FloatArray:
            return np.zeros(shape, dtype=np.float64)

        return cls(
            time=time_array.copy(),
            plastic_strain_rate=zero_field(),
            elastic_strain=zero_field(),
            alpha_rate=zero_field(),
            r_bar_rate=zero_field(),
            damage_rate=zero_field(),
            stress=zero_field(),
            beta=zero_field(),
            R_bar=zero_field(),
            energy_release_rate=zero_field(),
            plastic_strain=zero_field(),
            alpha=zero_field(),
            r_bar=zero_field(),
            damage=zero_field(),
        )

    @classmethod
    def from_elastic_initialization(
        cls,
        initialization: ElasticInitialization,
        materials: Sequence[MaterialParameters],
    ) -> "LatinState":
        """
        Convert the globally admissible elastic solution into s0 in A.

        Plastic strain, hardening variables and damage are zero. Stress and
        elastic strain come from the elastic solution. The damage energy
        release rate Y is evaluated from the elastic stress with D = 0.
        """
        _, n_elements = initialization.stress.shape

        if len(materials) != n_elements:
            raise ValueError(
                "One MaterialParameters object is required per element."
            )

        state = cls.zeros(
            time=initialization.time,
            n_elements=n_elements,
        )

        state.elastic_strain[:, :] = initialization.strain
        state.stress[:, :] = initialization.stress

        for element, material in enumerate(materials):
            stress = state.stress[:, element]
            tensile = stress >= 0.0

            state.energy_release_rate[:, element] = np.where(
                tensile,
                stress**2 / (2.0 * material.E),
                material.h * stress**2 / (2.0 * material.E),
            )

        return state

# -*- coding: utf-8 -*-
"""
Canonical material-point layout and immutable state representation for the
offshore-wind-turbine LATIN-PGD solver.

The tower material-point coordinate is

    q <-> (element, Gauss point, local fiber),

with element-major, Gauss-major, fiber-major ordering.  All LATIN material
histories use canonical shape

    (n_time, n_material_points).

This module is intentionally independent of FEM topology objects, constitutive
models, PGD bases, and equilibrium operators.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import ClassVar, Tuple

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


def _validated_positive_integer(value: int, name: str) -> int:
    """Return a validated strictly positive integer."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")

    result = int(value)
    if result < 1:
        raise ValueError(name + " must be at least 1.")
    return result


def _validated_index(
    value: int,
    upper_bound: int,
    name: str,
) -> int:
    """Return an integer index in [0, upper_bound)."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")

    result = int(value)
    if result < 0 or result >= upper_bound:
        raise IndexError(
            name
            + " must satisfy 0 <= "
            + name
            + " < "
            + str(upper_bound)
            + "."
        )
    return result


def _readonly_float_array(
    values: FloatArray,
    name: str,
    ndim: int,
) -> FloatArray:
    """Return a detached finite float64 array with write access disabled."""
    array = np.array(values, dtype=np.float64, copy=True)

    if array.ndim != ndim:
        raise ValueError(
            name
            + " must be "
            + ("one-dimensional." if ndim == 1 else "two-dimensional.")
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must contain only finite values.")

    array.setflags(write=False)
    return array


@dataclass(frozen=True)
class MaterialPointLayout:
    """
    Canonical mapping between q and (element, Gauss point, local fiber).

    Ordering is element-major, then Gauss-major, then fiber-major:

        q = (element_index * n_gauss + gauss_index) * n_fibers
            + fiber_index.
    """

    n_elements: int
    n_gauss: int
    n_fibers: int

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "n_elements",
            _validated_positive_integer(
                self.n_elements,
                "n_elements",
            ),
        )
        object.__setattr__(
            self,
            "n_gauss",
            _validated_positive_integer(
                self.n_gauss,
                "n_gauss",
            ),
        )
        object.__setattr__(
            self,
            "n_fibers",
            _validated_positive_integer(
                self.n_fibers,
                "n_fibers",
            ),
        )

    @property
    def n_material_points(self) -> int:
        """Return the total number of beam-Gauss-fiber material points."""
        return int(
            self.n_elements
            * self.n_gauss
            * self.n_fibers
        )

    def flatten(
        self,
        element_index: int,
        gauss_index: int,
        fiber_index: int,
    ) -> int:
        """Return canonical q for one (element, Gauss, fiber) triple."""
        element = _validated_index(
            element_index,
            self.n_elements,
            "element_index",
        )
        gauss = _validated_index(
            gauss_index,
            self.n_gauss,
            "gauss_index",
        )
        fiber = _validated_index(
            fiber_index,
            self.n_fibers,
            "fiber_index",
        )

        return int(
            (element * self.n_gauss + gauss)
            * self.n_fibers
            + fiber
        )

    def unflatten(
        self,
        material_point_index: int,
    ) -> Tuple[int, int, int]:
        """Return (element, Gauss, fiber) for canonical material point q."""
        q = _validated_index(
            material_point_index,
            self.n_material_points,
            "material_point_index",
        )

        points_per_element = self.n_gauss * self.n_fibers
        element = q // points_per_element
        remainder = q % points_per_element
        gauss = remainder // self.n_fibers
        fiber = remainder % self.n_fibers

        return int(element), int(gauss), int(fiber)

    @property
    def element_indices(self) -> IntArray:
        """Return the element index associated with each canonical q."""
        result = np.repeat(
            np.arange(self.n_elements, dtype=np.int64),
            self.n_gauss * self.n_fibers,
        )
        result.setflags(write=False)
        return result

    @property
    def gauss_indices(self) -> IntArray:
        """Return the local Gauss-point index associated with each q."""
        one_element = np.repeat(
            np.arange(self.n_gauss, dtype=np.int64),
            self.n_fibers,
        )
        result = np.tile(one_element, self.n_elements)
        result.setflags(write=False)
        return result

    @property
    def fiber_indices(self) -> IntArray:
        """Return the local fiber index associated with each canonical q."""
        result = np.tile(
            np.arange(self.n_fibers, dtype=np.int64),
            self.n_elements * self.n_gauss,
        )
        result.setflags(write=False)
        return result


@dataclass(frozen=True)
class LatinStateTower:
    """
    Immutable complete LATIN state on the tower time-material-point grid.

    The 13 material fields are stored with common shape
    (n_time, n_material_points).  The object owns detached read-only copies of
    all input arrays so that accepted states cannot be polluted through writable
    aliases during local/global trial construction.
    """

    time: FloatArray

    # Primary LATIN fields
    plastic_strain_rate: FloatArray
    elastic_strain: FloatArray
    alpha_rate: FloatArray
    r_bar_rate: FloatArray
    damage_rate: FloatArray
    stress: FloatArray
    beta: FloatArray
    R_bar: FloatArray
    energy_release_rate: FloatArray

    # Explicit integrated support histories
    plastic_strain: FloatArray
    alpha: FloatArray
    r_bar: FloatArray
    damage: FloatArray

    MATERIAL_FIELD_NAMES: ClassVar[Tuple[str, ...]] = (
        "plastic_strain_rate",
        "elastic_strain",
        "alpha_rate",
        "r_bar_rate",
        "damage_rate",
        "stress",
        "beta",
        "R_bar",
        "energy_release_rate",
        "plastic_strain",
        "alpha",
        "r_bar",
        "damage",
    )

    def __post_init__(self) -> None:
        time = _readonly_float_array(
            self.time,
            "time",
            ndim=1,
        )
        if time.size < 2:
            raise ValueError("At least two time points are required.")
        if np.any(np.diff(time) <= 0.0):
            raise ValueError("time must be strictly increasing.")

        expected_shape = None

        for field_name in self.MATERIAL_FIELD_NAMES:
            array = _readonly_float_array(
                getattr(self, field_name),
                field_name,
                ndim=2,
            )

            if expected_shape is None:
                expected_shape = array.shape
                if expected_shape[1] < 1:
                    raise ValueError(
                        "The tower LATIN state must contain at least "
                        "one material point."
                    )
            elif array.shape != expected_shape:
                raise ValueError(
                    field_name
                    + " has shape "
                    + str(array.shape)
                    + "; expected "
                    + str(expected_shape)
                    + "."
                )

            if array.shape[0] != time.size:
                raise ValueError(
                    field_name
                    + " must contain one row per time point."
                )

            object.__setattr__(self, field_name, array)

        damage = self.damage
        if np.any(damage < 0.0) or np.any(damage >= 1.0):
            raise ValueError(
                "damage must satisfy 0 <= damage < 1."
            )

        object.__setattr__(self, "time", time)

    @property
    def n_time(self) -> int:
        """Return the number of time points."""
        return int(self.time.size)

    @property
    def n_material_points(self) -> int:
        """Return the number of canonical beam-Gauss-fiber material points."""
        return int(self.stress.shape[1])

    @property
    def field_shape(self) -> Tuple[int, int]:
        """Return common (n_time, n_material_points) field shape."""
        return (
            int(self.stress.shape[0]),
            int(self.stress.shape[1]),
        )

    def copy(self) -> "LatinStateTower":
        """Return a fully detached state value with identical contents."""
        return LatinStateTower(
            time=self.time,
            plastic_strain_rate=self.plastic_strain_rate,
            elastic_strain=self.elastic_strain,
            alpha_rate=self.alpha_rate,
            r_bar_rate=self.r_bar_rate,
            damage_rate=self.damage_rate,
            stress=self.stress,
            beta=self.beta,
            R_bar=self.R_bar,
            energy_release_rate=self.energy_release_rate,
            plastic_strain=self.plastic_strain,
            alpha=self.alpha,
            r_bar=self.r_bar,
            damage=self.damage,
        )

    @classmethod
    def zeros(
        cls,
        time: FloatArray,
        n_material_points: int,
    ) -> "LatinStateTower":
        """Return a zero state on a specified time-material-point grid."""
        n_points = _validated_positive_integer(
            n_material_points,
            "n_material_points",
        )
        time_array = np.asarray(time, dtype=np.float64)

        shape = (time_array.size, n_points)

        def zero_field() -> FloatArray:
            return np.zeros(shape, dtype=np.float64)

        return cls(
            time=time_array,
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

# -*- coding: utf-8 -*-
"""
Fiber discretization utilities for thin-walled circular annular sections.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray


FloatArray = NDArray[np.float64]
IntArray = NDArray[np.int64]


@dataclass(frozen=True)
class AnnularFiberSection:
    """
    Polar fiber discretization of a circular annular cross-section.

    Each fiber represents one annular sector. Its area and centroid are
    evaluated analytically rather than approximated by a rectangular cell.
    """

    outer_radius: float
    inner_radius: float
    n_circumferential: int
    n_radial: int
    areas: FloatArray
    x_coordinates: FloatArray
    y_coordinates: FloatArray
    radial_indices: IntArray
    circumferential_indices: IntArray

    def __post_init__(self) -> None:
        outer_radius = float(self.outer_radius)
        inner_radius = float(self.inner_radius)
        n_circumferential = int(self.n_circumferential)
        n_radial = int(self.n_radial)

        areas = np.asarray(self.areas, dtype=np.float64)
        x_coordinates = np.asarray(
            self.x_coordinates,
            dtype=np.float64,
        )
        y_coordinates = np.asarray(
            self.y_coordinates,
            dtype=np.float64,
        )
        radial_indices = np.asarray(
            self.radial_indices,
            dtype=np.int64,
        )
        circumferential_indices = np.asarray(
            self.circumferential_indices,
            dtype=np.int64,
        )

        if not np.isfinite(outer_radius) or outer_radius <= 0.0:
            raise ValueError("outer_radius must be finite and positive.")
        if not np.isfinite(inner_radius) or inner_radius < 0.0:
            raise ValueError(
                "inner_radius must be finite and non-negative."
            )
        if inner_radius >= outer_radius:
            raise ValueError(
                "inner_radius must be smaller than outer_radius."
            )

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

        expected_size = n_circumferential * n_radial
        arrays = (
            areas,
            x_coordinates,
            y_coordinates,
            radial_indices,
            circumferential_indices,
        )

        for array in arrays:
            if array.ndim != 1:
                raise ValueError("Fiber arrays must be one-dimensional.")
            if array.size != expected_size:
                raise ValueError(
                    "Fiber arrays must contain "
                    "n_circumferential * n_radial entries."
                )

        if np.any(~np.isfinite(areas)) or np.any(areas <= 0.0):
            raise ValueError(
                "All fiber areas must be finite and positive."
            )
        if np.any(~np.isfinite(x_coordinates)):
            raise ValueError(
                "All fiber x-coordinates must be finite."
            )
        if np.any(~np.isfinite(y_coordinates)):
            raise ValueError(
                "All fiber y-coordinates must be finite."
            )

        if np.any(radial_indices < 0):
            raise ValueError("radial_indices must be non-negative.")
        if np.any(radial_indices >= n_radial):
            raise ValueError(
                "radial_indices exceed the radial layer count."
            )
        if np.any(circumferential_indices < 0):
            raise ValueError(
                "circumferential_indices must be non-negative."
            )
        if np.any(
            circumferential_indices >= n_circumferential
        ):
            raise ValueError(
                "circumferential_indices exceed the sector count."
            )

        object.__setattr__(self, "outer_radius", outer_radius)
        object.__setattr__(self, "inner_radius", inner_radius)
        object.__setattr__(
            self,
            "n_circumferential",
            n_circumferential,
        )
        object.__setattr__(self, "n_radial", n_radial)
        object.__setattr__(self, "areas", areas)
        object.__setattr__(
            self,
            "x_coordinates",
            x_coordinates,
        )
        object.__setattr__(
            self,
            "y_coordinates",
            y_coordinates,
        )
        object.__setattr__(
            self,
            "radial_indices",
            radial_indices,
        )
        object.__setattr__(
            self,
            "circumferential_indices",
            circumferential_indices,
        )

    @property
    def n_fibers(self) -> int:
        """Return the total number of section fibers."""
        return int(self.areas.size)

    @property
    def area(self) -> float:
        """Return the area represented by all fibers."""
        return float(np.sum(self.areas))

    @property
    def exact_area(self) -> float:
        """Return the analytical annular-section area."""
        return float(
            np.pi
            * (
                self.outer_radius ** 2
                - self.inner_radius ** 2
            )
        )

    @property
    def centroid_x(self) -> float:
        """Return the fiber-discretized x-coordinate of the centroid."""
        return float(
            np.dot(self.areas, self.x_coordinates) / self.area
        )

    @property
    def centroid_y(self) -> float:
        """Return the fiber-discretized y-coordinate of the centroid."""
        return float(
            np.dot(self.areas, self.y_coordinates) / self.area
        )

    @property
    def second_moment_x(self) -> float:
        """
        Return the fiber-discretized second moment about the x-axis.

        Fibers are represented as concentrated areas at their exact
        annular-sector centroids.
        """
        return float(
            np.dot(self.areas, self.y_coordinates ** 2)
        )

    @property
    def second_moment_y(self) -> float:
        """
        Return the fiber-discretized second moment about the y-axis.

        Fibers are represented as concentrated areas at their exact
        annular-sector centroids.
        """
        return float(
            np.dot(self.areas, self.x_coordinates ** 2)
        )

    @property
    def product_moment_xy(self) -> float:
        """Return the fiber-discretized product moment of area."""
        return float(
            np.dot(
                self.areas,
                self.x_coordinates * self.y_coordinates,
            )
        )

    @property
    def exact_second_moment(self) -> float:
        """
        Return the analytical second moment about either centroidal axis.
        """
        return float(
            0.25
            * np.pi
            * (
                self.outer_radius ** 4
                - self.inner_radius ** 4
            )
        )

    @property
    def coordinates(self) -> FloatArray:
        """Return fiber centroid coordinates with shape (n_fibers, 2)."""
        return np.column_stack(
            (
                self.x_coordinates,
                self.y_coordinates,
            )
        )


def _validated_integer(value: int, name: str) -> int:
    """Return a validated integer discretization parameter."""
    if isinstance(value, (bool, np.bool_)):
        raise TypeError(name + " must be an integer.")
    if not isinstance(value, (int, np.integer)):
        raise TypeError(name + " must be an integer.")
    return int(value)


def create_annular_fiber_section(
    outer_diameter: float,
    thickness: float,
    n_circumferential: int = 32,
    n_radial: int = 2,
) -> AnnularFiberSection:
    """
    Create an annular-sector fiber discretization.

    Circumferential sector centroids are placed at angles

        0, delta_theta, 2 * delta_theta, ...

    so that fibers are aligned with the positive and negative coordinate
    axes whenever n_circumferential is divisible by four. Radial layers
    have equal thickness.

    Parameters
    ----------
    outer_diameter
        Outside diameter of the circular section.
    thickness
        Wall thickness of the circular section.
    n_circumferential
        Number of equal circumferential sectors.
    n_radial
        Number of equal radial layers through the wall thickness.
    """
    outer_diameter = float(outer_diameter)
    thickness = float(thickness)
    n_circumferential = _validated_integer(
        n_circumferential,
        "n_circumferential",
    )
    n_radial = _validated_integer(
        n_radial,
        "n_radial",
    )

    if not np.isfinite(outer_diameter) or outer_diameter <= 0.0:
        raise ValueError(
            "outer_diameter must be finite and positive."
        )
    if not np.isfinite(thickness) or thickness <= 0.0:
        raise ValueError("thickness must be finite and positive.")

    outer_radius = 0.5 * outer_diameter
    inner_radius = outer_radius - thickness

    if inner_radius <= 0.0:
        raise ValueError(
            "thickness must be smaller than the outer radius."
        )
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

    delta_theta = 2.0 * np.pi / n_circumferential

    radial_edges = np.linspace(
        inner_radius,
        outer_radius,
        n_radial + 1,
        dtype=np.float64,
    )

    radial_inner = np.repeat(
        radial_edges[:-1],
        n_circumferential,
    )
    radial_outer = np.repeat(
        radial_edges[1:],
        n_circumferential,
    )

    circumferential_indices = np.tile(
        np.arange(n_circumferential, dtype=np.int64),
        n_radial,
    )
    radial_indices = np.repeat(
        np.arange(n_radial, dtype=np.int64),
        n_circumferential,
    )

    centroid_angles = (
        circumferential_indices.astype(np.float64)
        * delta_theta
    )

    areas = (
        0.5
        * (
            radial_outer ** 2
            - radial_inner ** 2
        )
        * delta_theta
    )

    centroid_radii = (
        4.0
        * np.sin(0.5 * delta_theta)
        / (3.0 * delta_theta)
        * (
            radial_outer ** 3
            - radial_inner ** 3
        )
        / (
            radial_outer ** 2
            - radial_inner ** 2
        )
    )

    x_coordinates = centroid_radii * np.cos(
        centroid_angles
    )
    y_coordinates = centroid_radii * np.sin(
        centroid_angles
    )

    return AnnularFiberSection(
        outer_radius=outer_radius,
        inner_radius=inner_radius,
        n_circumferential=n_circumferential,
        n_radial=n_radial,
        areas=areas,
        x_coordinates=x_coordinates,
        y_coordinates=y_coordinates,
        radial_indices=radial_indices,
        circumferential_indices=circumferential_indices,
    )

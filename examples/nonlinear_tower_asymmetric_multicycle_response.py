# -*- coding: utf-8 -*-
"""
Asymmetric multi-cycle nonlinear response of the NREL 5 MW tower.

The formal benchmark uses a positive-mean, sign-reversing force history with
R_F = F_min / F_max = -0.5. Its periodic checkpoints are

    F_mean -> F_max -> F_mean -> F_min -> F_mean.

A separate preload 0 -> F_mean is handled by the verified asymmetric driver
and is not counted as a cycle.

For ratcheting assessment, the direct global metric is the displacement change
over one complete cycle at the same reference force F_mean:

    Delta u_cycle = u_end - u_start.

The corresponding fixed-fiber material metric is the signed plastic-strain
increment over the same cycle.
"""

from __future__ import annotations

from typing import Tuple

import numpy as np
from numpy.typing import NDArray

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_response import (
    run_nonlinear_asymmetric_analysis,
)
from examples.nonlinear_tower_multicycle_diagnostics import (
    extract_multicycle_diagnostics,
)
from examples.nonlinear_tower_multicycle_response import (
    MulticycleTowerResult,
    evaluate_cycle_similarities,
)
from fem.tower_loading import (
    AsymmetricCyclicTopForceHistory,
    create_asymmetric_cyclic_top_force_history,
)
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


def run_asymmetric_multicycle_tower_analysis(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    maximum_force: float = 1.0e6,
    force_ratio: float = -0.5,
    period: float = 10.0,
    n_cycles: int = 5,
    increments_per_cycle: int = 40,
    similarity_tolerance: float = 1.0e-3,
    max_iterations: int = 40,
) -> MulticycleTowerResult:
    """Run asymmetric cyclic tower analysis and extract per-cycle metrics."""
    loading = create_asymmetric_cyclic_top_force_history(
        maximum_force=maximum_force,
        force_ratio=force_ratio,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
    )
    response = run_nonlinear_asymmetric_analysis(
        configuration=configuration,
        material=material,
        loading=loading,
        max_iterations=max_iterations,
    )
    diagnostics = extract_multicycle_diagnostics(response=response)
    similarities = evaluate_cycle_similarities(
        diagnostics=diagnostics,
        tolerance=similarity_tolerance,
    )

    return MulticycleTowerResult(
        response=response,
        diagnostics=diagnostics,
        similarities=similarities,
        similarity_tolerance=similarity_tolerance,
    )


def _validate_asymmetric_result(result: MulticycleTowerResult) -> None:
    """Validate that a multicycle result uses asymmetric cyclic loading."""
    if not isinstance(result, MulticycleTowerResult):
        raise TypeError("result must be a MulticycleTowerResult.")
    if not isinstance(
        result.response.loading,
        AsymmetricCyclicTopForceHistory,
    ):
        raise TypeError(
            "result must use AsymmetricCyclicTopForceHistory."
        )
    if result.diagnostics.n_cycles != result.response.loading.n_cycles:
        raise ValueError(
            "diagnostics must cover every asymmetric loading cycle."
        )


def cycle_displacement_drifts(
    result: MulticycleTowerResult,
) -> FloatArray:
    """Return u_end - u_start at the same F_mean for every complete cycle."""
    _validate_asymmetric_result(result)
    response = result.response
    values = np.empty(result.diagnostics.n_cycles, dtype=np.float64)

    for index, cycle in enumerate(result.diagnostics.cycles):
        values[index] = (
            response.top_displacements[cycle.indices.end]
            - response.top_displacements[cycle.indices.start]
        )

    return values


def normalized_cycle_displacement_drifts(
    result: MulticycleTowerResult,
    scale_floor: float = 1.0e-14,
) -> FloatArray:
    """Return |Delta u_cycle| normalized by each cycle displacement range."""
    drifts = cycle_displacement_drifts(result)
    ranges = result.diagnostics.displacement_ranges

    scale_floor = float(scale_floor)
    if not np.isfinite(scale_floor):
        raise ValueError("scale_floor must be finite.")
    if scale_floor <= 0.0:
        raise ValueError("scale_floor must be positive.")

    denominators = np.maximum(np.abs(ranges), scale_floor)
    return np.abs(drifts) / denominators


def critical_plastic_strain_drifts(
    result: MulticycleTowerResult,
) -> FloatArray:
    """Return signed fixed-fiber plastic-strain increment for every cycle."""
    _validate_asymmetric_result(result)
    return np.asarray(
        result.diagnostics.critical_plastic_strain_increments,
        dtype=np.float64,
    ).copy()


def ratcheting_metrics(
    result: MulticycleTowerResult,
) -> Tuple[FloatArray, FloatArray, FloatArray]:
    """Return global drift, normalized global drift, and plastic drift."""
    return (
        cycle_displacement_drifts(result),
        normalized_cycle_displacement_drifts(result),
        critical_plastic_strain_drifts(result),
    )

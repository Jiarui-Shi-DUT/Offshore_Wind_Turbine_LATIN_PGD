# -*- coding: utf-8 -*-
"""
Damage-mechanism separation for asymmetric multi-cycle tower loading.

The two analyses use the same tower discretisation, asymmetric loading,
plasticity, hardening, and time integration:

    1. coupled viscoplastic-damage response;
    2. damage-disabled viscoplastic response with k_damage = 0.

The formal asymmetric benchmark uses

    R_F = F_min / F_max = -0.5

with periodic checkpoints

    F_mean -> F_max -> F_mean -> F_min -> F_mean.

The common critical fiber is selected from the coupled analysis and imposed on
the damage-disabled response before cycle diagnostics are extracted.

In addition to the existing damage-mechanism comparison quantities, this
module compares two same-force ratcheting diagnostics:

    Delta u_n =
        u_end,n - u_start,n

and

    Delta eps_p,n =
        eps_p,end,n - eps_p,start,n.

Both cycle endpoints correspond to F_mean. Persistent non-zero drift over many
cycles is required before interpreting these metrics as sustained ratcheting.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from numpy.typing import NDArray

from examples.elastic_tapered_tower import TowerConfiguration
from examples.nonlinear_tower_asymmetric_multicycle_response import (
    critical_plastic_strain_drifts,
    cycle_displacement_drifts,
    run_asymmetric_multicycle_tower_analysis,
)
from examples.nonlinear_tower_damage_mechanism_probe import (
    DamageMechanismComparison,
    build_multicycle_result,
    compare_cycle_diagnostics,
    response_with_common_critical_location,
)
from material.viscoplastic_damage_1d import MaterialParameters


FloatArray = NDArray[np.float64]


def _validated_history(
    values: FloatArray,
    n_cycles: int,
    name: str,
) -> FloatArray:
    """Return one validated finite per-cycle history."""
    array = np.asarray(values, dtype=np.float64)
    if array.shape != (n_cycles,):
        raise ValueError(
            name + " must contain one value per cycle."
        )
    if np.any(~np.isfinite(array)):
        raise ValueError(name + " must be finite.")
    return array.copy()


@dataclass(frozen=True)
class AsymmetricDamageMechanismComparison:
    """Paired asymmetric analyses with ratcheting-drift diagnostics."""

    mechanism: DamageMechanismComparison

    coupled_cycle_displacement_drifts: FloatArray
    damage_disabled_cycle_displacement_drifts: FloatArray
    cycle_displacement_drift_differences: FloatArray

    coupled_plastic_strain_drifts: FloatArray
    damage_disabled_plastic_strain_drifts: FloatArray
    plastic_strain_drift_differences: FloatArray

    def __post_init__(self) -> None:
        if not isinstance(
            self.mechanism,
            DamageMechanismComparison,
        ):
            raise TypeError(
                "mechanism must be a DamageMechanismComparison."
            )

        n_cycles = self.mechanism.n_cycles
        history_names = (
            "coupled_cycle_displacement_drifts",
            "damage_disabled_cycle_displacement_drifts",
            "cycle_displacement_drift_differences",
            "coupled_plastic_strain_drifts",
            "damage_disabled_plastic_strain_drifts",
            "plastic_strain_drift_differences",
        )
        for name in history_names:
            values = _validated_history(
                getattr(self, name),
                n_cycles,
                name,
            )
            object.__setattr__(self, name, values)

        expected_displacement_difference = (
            self.coupled_cycle_displacement_drifts
            - self.damage_disabled_cycle_displacement_drifts
        )
        if not np.allclose(
            self.cycle_displacement_drift_differences,
            expected_displacement_difference,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "cycle_displacement_drift_differences "
                "is inconsistent."
            )

        expected_plastic_difference = (
            self.coupled_plastic_strain_drifts
            - self.damage_disabled_plastic_strain_drifts
        )
        if not np.allclose(
            self.plastic_strain_drift_differences,
            expected_plastic_difference,
            rtol=1.0e-12,
            atol=1.0e-15,
        ):
            raise ValueError(
                "plastic_strain_drift_differences is inconsistent."
            )

        coupled_loading = self.mechanism.coupled.response.loading
        disabled_loading = (
            self.mechanism.damage_disabled.response.loading
        )
        if coupled_loading.force_ratio != disabled_loading.force_ratio:
            raise ValueError(
                "Paired asymmetric analyses must use "
                "the same force ratio."
            )
        if coupled_loading.maximum_force != disabled_loading.maximum_force:
            raise ValueError(
                "Paired asymmetric analyses must use "
                "the same maximum force."
            )

    @property
    def n_cycles(self) -> int:
        """Return number of compared cycles."""
        return self.mechanism.n_cycles

    @property
    def cycle_numbers(self) -> NDArray[np.int64]:
        """Return one-based cycle numbers."""
        return self.mechanism.cycle_numbers


def run_asymmetric_damage_mechanism_comparison(
    configuration: TowerConfiguration,
    material: MaterialParameters,
    maximum_force: float = 1.0e6,
    force_ratio: float = -0.5,
    period: float = 10.0,
    n_cycles: int = 5,
    increments_per_cycle: int = 40,
    similarity_tolerance: float = 1.0e-3,
    max_iterations: int = 40,
) -> AsymmetricDamageMechanismComparison:
    """
    Run paired asymmetric coupled and damage-disabled tower analyses.

    All settings other than k_damage are kept identical.
    """
    if not isinstance(
        configuration,
        TowerConfiguration,
    ):
        raise TypeError(
            "configuration must be a TowerConfiguration."
        )
    if not isinstance(
        material,
        MaterialParameters,
    ):
        raise TypeError(
            "material must be a MaterialParameters object."
        )

    coupled = run_asymmetric_multicycle_tower_analysis(
        configuration=configuration,
        material=material,
        maximum_force=maximum_force,
        force_ratio=force_ratio,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
        similarity_tolerance=similarity_tolerance,
        max_iterations=max_iterations,
    )

    disabled_material = replace(
        material,
        k_damage=0.0,
    )
    raw_disabled = run_asymmetric_multicycle_tower_analysis(
        configuration=configuration,
        material=disabled_material,
        maximum_force=maximum_force,
        force_ratio=force_ratio,
        period=period,
        n_cycles=n_cycles,
        increments_per_cycle=increments_per_cycle,
        similarity_tolerance=similarity_tolerance,
        max_iterations=max_iterations,
    )

    common_disabled_response = (
        response_with_common_critical_location(
            response=raw_disabled.response,
            reference=coupled.response,
        )
    )
    damage_disabled = build_multicycle_result(
        response=common_disabled_response,
        similarity_tolerance=similarity_tolerance,
    )

    mechanism = DamageMechanismComparison(
        coupled=coupled,
        damage_disabled=damage_disabled,
        cycles=compare_cycle_diagnostics(
            coupled=coupled.diagnostics,
            damage_disabled=damage_disabled.diagnostics,
        ),
    )

    coupled_displacement_drifts = (
        cycle_displacement_drifts(coupled)
    )
    disabled_displacement_drifts = (
        cycle_displacement_drifts(damage_disabled)
    )
    coupled_plastic_drifts = (
        critical_plastic_strain_drifts(coupled)
    )
    disabled_plastic_drifts = (
        critical_plastic_strain_drifts(damage_disabled)
    )

    return AsymmetricDamageMechanismComparison(
        mechanism=mechanism,
        coupled_cycle_displacement_drifts=(
            coupled_displacement_drifts
        ),
        damage_disabled_cycle_displacement_drifts=(
            disabled_displacement_drifts
        ),
        cycle_displacement_drift_differences=(
            coupled_displacement_drifts
            - disabled_displacement_drifts
        ),
        coupled_plastic_strain_drifts=(
            coupled_plastic_drifts
        ),
        damage_disabled_plastic_strain_drifts=(
            disabled_plastic_drifts
        ),
        plastic_strain_drift_differences=(
            coupled_plastic_drifts
            - disabled_plastic_drifts
        ),
    )

# -*- coding: utf-8 -*-
"""
Tower-compatible nonlinear LATIN local stage on canonical material points.

The constitutive algorithm is intentionally the same as the validated
one-dimensional local stage.  The only structural generalisation is the space
index:

    element  ->  canonical material point q.

For the ascent-direction choice used by the project,

    (B_plus)^(-1) = 0,
    (b_plus)^(-1) = 0,

the thermodynamic-force histories remain fixed during the local projection:

    sigma_hat = sigma_i,
    beta_hat = beta_i,
    R_bar_hat = R_bar_i,
    Y_hat = Y_i.

The four integrated internal histories [eps_p, alpha, r_bar, D] are advanced
sequentially in time and independently for every q-point using the same RK4
material-point integrator as latin.local_stage.  The returned LatinStateTower
is a new immutable value; the accepted input state is never used as scratch
storage.
"""

from __future__ import annotations

from typing import Sequence, Tuple, Union

import numpy as np

from latin.local_stage import (
    _integrate_one_local_step,
    local_rates_from_forces,
    unilateral_elastic_strain,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


MaterialInput = Union[
    MaterialParameters,
    Sequence[MaterialParameters],
]


def _material_sequence(
    materials: MaterialInput,
    n_material_points: int,
) -> Tuple[MaterialParameters, ...]:
    """Broadcast one material or validate one material object per q-point."""
    if isinstance(materials, MaterialParameters):
        material_tuple = (materials,) * n_material_points
    else:
        material_tuple = tuple(materials)
        if len(material_tuple) == 1:
            material_tuple = material_tuple * n_material_points

    if len(material_tuple) != n_material_points:
        raise ValueError(
            "materials must contain one MaterialParameters object or one "
            "object per material point."
        )
    if any(
        not isinstance(material, MaterialParameters)
        for material in material_tuple
    ):
        raise TypeError(
            "materials must contain only MaterialParameters objects."
        )

    return material_tuple


def solve_tower_local_stage(
    global_state: LatinStateTower,
    materials: MaterialInput,
) -> LatinStateTower:
    """
    Project one immutable tower global state from A onto local manifold Gamma.

    Time remains sequential inside each material point because the integrated
    constitutive variables carry history.  Distinct q-points are independent.
    """
    if not isinstance(global_state, LatinStateTower):
        raise TypeError(
            "global_state must be a LatinStateTower."
        )

    material_tuple = _material_sequence(
        materials=materials,
        n_material_points=global_state.n_material_points,
    )

    fields = {
        name: np.array(
            getattr(global_state, name),
            dtype=np.float64,
            copy=True,
        )
        for name in LatinStateTower.MATERIAL_FIELD_NAMES
    }

    # The ascent search direction keeps the thermodynamic forces fixed.
    fields["stress"][:, :] = global_state.stress
    fields["beta"][:, :] = global_state.beta
    fields["R_bar"][:, :] = global_state.R_bar
    fields["energy_release_rate"][:, :] = (
        global_state.energy_release_rate
    )

    for q, material in enumerate(material_tuple):
        internal_state = np.array(
            [
                global_state.plastic_strain[0, q],
                global_state.alpha[0, q],
                global_state.r_bar[0, q],
                global_state.damage[0, q],
            ],
            dtype=np.float64,
        )

        fields["plastic_strain"][0, q] = internal_state[0]
        fields["alpha"][0, q] = internal_state[1]
        fields["r_bar"][0, q] = internal_state[2]
        fields["damage"][0, q] = internal_state[3]

        initial_rates = local_rates_from_forces(
            stress=float(fields["stress"][0, q]),
            beta=float(fields["beta"][0, q]),
            transformed_force=float(fields["R_bar"][0, q]),
            energy_release_rate=float(
                fields["energy_release_rate"][0, q]
            ),
            damage=float(fields["damage"][0, q]),
            material=material,
        )

        (
            fields["plastic_strain_rate"][0, q],
            fields["alpha_rate"][0, q],
            fields["r_bar_rate"][0, q],
            fields["damage_rate"][0, q],
        ) = initial_rates

        fields["elastic_strain"][0, q] = unilateral_elastic_strain(
            stress=float(fields["stress"][0, q]),
            damage=float(fields["damage"][0, q]),
            material=material,
        )

        for step in range(1, global_state.n_time):
            time_step = float(
                global_state.time[step]
                - global_state.time[step - 1]
            )

            internal_state = _integrate_one_local_step(
                state_old=internal_state,
                time_step=time_step,
                stress_old=float(
                    fields["stress"][step - 1, q]
                ),
                stress_new=float(
                    fields["stress"][step, q]
                ),
                beta_old=float(
                    fields["beta"][step - 1, q]
                ),
                beta_new=float(
                    fields["beta"][step, q]
                ),
                transformed_force_old=float(
                    fields["R_bar"][step - 1, q]
                ),
                transformed_force_new=float(
                    fields["R_bar"][step, q]
                ),
                energy_old=float(
                    fields["energy_release_rate"][
                        step - 1,
                        q,
                    ]
                ),
                energy_new=float(
                    fields["energy_release_rate"][step, q]
                ),
                material=material,
            )

            fields["plastic_strain"][step, q] = internal_state[0]
            fields["alpha"][step, q] = internal_state[1]
            fields["r_bar"][step, q] = internal_state[2]
            fields["damage"][step, q] = internal_state[3]

            current_rates = local_rates_from_forces(
                stress=float(fields["stress"][step, q]),
                beta=float(fields["beta"][step, q]),
                transformed_force=float(
                    fields["R_bar"][step, q]
                ),
                energy_release_rate=float(
                    fields["energy_release_rate"][step, q]
                ),
                damage=float(fields["damage"][step, q]),
                material=material,
            )

            (
                fields["plastic_strain_rate"][step, q],
                fields["alpha_rate"][step, q],
                fields["r_bar_rate"][step, q],
                fields["damage_rate"][step, q],
            ) = current_rates

            fields["elastic_strain"][step, q] = (
                unilateral_elastic_strain(
                    stress=float(fields["stress"][step, q]),
                    damage=float(fields["damage"][step, q]),
                    material=material,
                )
            )

    return LatinStateTower(
        time=global_state.time,
        **fields,
    )

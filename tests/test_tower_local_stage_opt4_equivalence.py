# -*- coding: utf-8 -*-
"""OPT-4 equivalence test for homogeneous tower material points."""

import numpy as np

from latin.local_stage import solve_local_stage
from latin.state import LatinState
from latin.tower_local_stage import solve_tower_local_stage
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


def test_homogeneous_tower_matches_validated_1d_pointwise():
    """
    Verify the homogeneous tower local stage against the validated
    pointwise one-dimensional implementation.

    The case includes nonuniform time steps, tension-compression reversal,
    and nonzero q-dependent initial internal variables.
    """
    time = np.array(
        [0.0, 0.03, 0.09, 0.17, 0.28],
        dtype=np.float64,
    )

    n_points = 4
    shape = (time.size, n_points)

    stress = np.array(
        [
            [95.0, -70.0, 110.0, -85.0],
            [120.0, -92.0, 82.0, -55.0],
            [-60.0, 105.0, -75.0, 98.0],
            [88.0, -65.0, 125.0, -100.0],
            [-72.0, 115.0, -90.0, 108.0],
        ],
        dtype=np.float64,
    )

    beta = np.array(
        [
            [1.20, -0.80, 0.60, -1.10],
            [1.35, -0.65, 0.75, -0.95],
            [0.90, -0.40, 0.30, -0.70],
            [1.10, -0.55, 0.50, -0.85],
            [0.75, -0.25, 0.20, -0.60],
        ],
        dtype=np.float64,
    )

    R_bar = np.array(
        [
            [1.40, 1.55, 1.70, 1.85],
            [1.48, 1.63, 1.78, 1.93],
            [1.56, 1.71, 1.86, 2.01],
            [1.64, 1.79, 1.94, 2.09],
            [1.72, 1.87, 2.02, 2.17],
        ],
        dtype=np.float64,
    )

    energy = np.array(
        [
            [0.0260, 0.0270, 0.0280, 0.0290],
            [0.0300, 0.0310, 0.0320, 0.0330],
            [0.0285, 0.0295, 0.0305, 0.0315],
            [0.0340, 0.0350, 0.0360, 0.0370],
            [0.0320, 0.0330, 0.0340, 0.0350],
        ],
        dtype=np.float64,
    )

    zeros = np.zeros(shape, dtype=np.float64)

    plastic_strain = zeros.copy()
    alpha = zeros.copy()
    r_bar = zeros.copy()
    damage = zeros.copy()

    plastic_strain[0, :] = np.array(
        [1.0e-5, -2.0e-5, 3.0e-5, -4.0e-5],
        dtype=np.float64,
    )
    alpha[0, :] = np.array(
        [2.0e-5, -1.0e-5, 4.0e-5, -3.0e-5],
        dtype=np.float64,
    )
    r_bar[0, :] = np.array(
        [0.010, 0.015, 0.020, 0.025],
        dtype=np.float64,
    )
    damage[0, :] = np.array(
        [0.0010, 0.0020, 0.0030, 0.0040],
        dtype=np.float64,
    )

    fields = {
        "plastic_strain_rate": zeros.copy(),
        "elastic_strain": zeros.copy(),
        "alpha_rate": zeros.copy(),
        "r_bar_rate": zeros.copy(),
        "damage_rate": zeros.copy(),
        "stress": stress,
        "beta": beta,
        "R_bar": R_bar,
        "energy_release_rate": energy,
        "plastic_strain": plastic_strain,
        "alpha": alpha,
        "r_bar": r_bar,
        "damage": damage,
    }

    one_d_state = LatinState(
        time=time,
        **{
            name: np.array(values, copy=True)
            for name, values in fields.items()
        },
    )

    tower_state = LatinStateTower(
        time=time,
        **{
            name: np.array(values, copy=True)
            for name, values in fields.items()
        },
    )

    material = MaterialParameters()

    one_d_local = solve_local_stage(
        global_state=one_d_state,
        materials=(material,) * n_points,
    )

    tower_local = solve_tower_local_stage(
        global_state=tower_state,
        materials=material,
    )

    for name in LatinStateTower.MATERIAL_FIELD_NAMES:
        np.testing.assert_allclose(
            getattr(tower_local, name),
            getattr(one_d_local, name),
            rtol=0.0,
            atol=1.0e-14,
        )

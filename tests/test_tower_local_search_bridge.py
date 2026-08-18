# -*- coding: utf-8 -*-
"""Unit tests for the tower local-stage/search-direction bridge."""

import unittest

import numpy as np

from latin.local_stage import solve_local_stage
from latin.search_directions import (
    compute_descent_search_directions,
)
from latin.state import LatinState
from latin.tower_local_stage import (
    solve_tower_local_stage,
)
from latin.tower_search_directions import (
    compute_tower_descent_search_directions,
)
from latin.tower_state import LatinStateTower
from material.viscoplastic_damage_1d import MaterialParameters


def _common_fields(time, n_points):
    t = np.arange(time.size, dtype=np.float64)[:, None]
    q = np.arange(n_points, dtype=np.float64)[None, :]
    stress = 90.0 + 18.0 * t + 3.0 * q
    beta = 1.0 + 0.15 * t + 0.05 * q
    R_bar = 1.5 + 0.20 * t + 0.04 * q
    energy = 0.025 + 0.004 * t + 0.001 * q

    shape = stress.shape
    return {
        "plastic_strain_rate": np.zeros(shape),
        "elastic_strain": np.zeros(shape),
        "alpha_rate": np.zeros(shape),
        "r_bar_rate": np.zeros(shape),
        "damage_rate": np.zeros(shape),
        "stress": stress,
        "beta": beta,
        "R_bar": R_bar,
        "energy_release_rate": energy,
        "plastic_strain": np.zeros(shape),
        "alpha": np.zeros(shape),
        "r_bar": np.zeros(shape),
        "damage": np.zeros(shape),
    }


def _latin_state(time, fields):
    return LatinState(
        time=time,
        **{
            name: np.array(values, copy=True)
            for name, values in fields.items()
        },
    )


def _tower_state(time, fields):
    return LatinStateTower(
        time=time,
        **{
            name: np.array(values, copy=True)
            for name, values in fields.items()
        },
    )


class TestTowerLocalSearchBridge(unittest.TestCase):
    """Verify q-generalisation against the validated 1D implementation."""

    def test_tower_local_stage_matches_current_1d_pointwise(self) -> None:
        time = np.array(
            [0.0, 0.05, 0.12, 0.20],
            dtype=np.float64,
        )
        n_points = 3
        fields = _common_fields(time, n_points)
        materials = (
            MaterialParameters(E=134_000.0, C=5_500.0),
            MaterialParameters(E=150_000.0, C=6_000.0),
            MaterialParameters(E=175_000.0, C=6_500.0),
        )

        one_d_global = _latin_state(time, fields)
        tower_global = _tower_state(time, fields)

        one_d_local = solve_local_stage(
            global_state=one_d_global,
            materials=materials,
        )
        tower_local = solve_tower_local_stage(
            global_state=tower_global,
            materials=materials,
        )

        for name in LatinStateTower.MATERIAL_FIELD_NAMES:
            np.testing.assert_allclose(
                getattr(tower_local, name),
                getattr(one_d_local, name),
                rtol=0.0,
                atol=1.0e-14,
            )

    def test_tower_local_stage_preserves_input_and_fixed_forces(self) -> None:
        time = np.array(
            [0.0, 0.08, 0.16],
            dtype=np.float64,
        )
        fields = _common_fields(time, 2)
        global_state = _tower_state(time, fields)
        snapshots = {
            name: np.array(
                getattr(global_state, name),
                copy=True,
            )
            for name in LatinStateTower.MATERIAL_FIELD_NAMES
        }

        local_state = solve_tower_local_stage(
            global_state=global_state,
            materials=MaterialParameters(),
        )

        for name in LatinStateTower.MATERIAL_FIELD_NAMES:
            np.testing.assert_allclose(
                getattr(global_state, name),
                snapshots[name],
                rtol=0.0,
                atol=0.0,
            )

        for force_name in (
            "stress",
            "beta",
            "R_bar",
            "energy_release_rate",
        ):
            np.testing.assert_allclose(
                getattr(local_state, force_name),
                getattr(global_state, force_name),
                rtol=0.0,
                atol=0.0,
            )

        self.assertFalse(local_state.stress.flags.writeable)
        self.assertFalse(local_state.damage.flags.writeable)

    def test_tower_search_directions_match_current_1d_pointwise(self) -> None:
        time = np.array(
            [0.0, 0.05, 0.12, 0.20],
            dtype=np.float64,
        )
        n_points = 3
        fields = _common_fields(time, n_points)
        materials = (
            MaterialParameters(E=134_000.0, C=5_500.0),
            MaterialParameters(E=150_000.0, C=6_000.0),
            MaterialParameters(E=175_000.0, C=6_500.0),
        )

        one_d_local = solve_local_stage(
            global_state=_latin_state(time, fields),
            materials=materials,
        )
        tower_local = solve_tower_local_stage(
            global_state=_tower_state(time, fields),
            materials=materials,
        )

        one_d_directions = compute_descent_search_directions(
            local_state=one_d_local,
            materials=materials,
            regularization=0.15,
        )
        tower_directions = (
            compute_tower_descent_search_directions(
                local_state=tower_local,
                materials=materials,
                regularization=0.15,
            )
        )

        for name in (
            "H_sigma",
            "H_beta",
            "H_R_bar",
            "b_damage",
        ):
            np.testing.assert_allclose(
                getattr(tower_directions, name),
                getattr(one_d_directions, name),
                rtol=0.0,
                atol=1.0e-14,
            )

    def test_single_material_broadcast_and_invalid_material_count(self) -> None:
        time = np.array(
            [0.0, 0.1, 0.2],
            dtype=np.float64,
        )
        fields = _common_fields(time, 4)
        global_state = _tower_state(time, fields)
        material = MaterialParameters()

        broadcast_local = solve_tower_local_stage(
            global_state=global_state,
            materials=material,
        )
        explicit_local = solve_tower_local_stage(
            global_state=global_state,
            materials=(material,) * 4,
        )

        for name in LatinStateTower.MATERIAL_FIELD_NAMES:
            np.testing.assert_allclose(
                getattr(broadcast_local, name),
                getattr(explicit_local, name),
                rtol=0.0,
                atol=0.0,
            )

        with self.assertRaises(ValueError):
            solve_tower_local_stage(
                global_state=global_state,
                materials=(material, material),
            )

        with self.assertRaises(ValueError):
            compute_tower_descent_search_directions(
                local_state=broadcast_local,
                materials=(material, material),
            )


if __name__ == "__main__":
    unittest.main()

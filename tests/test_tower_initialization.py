# -*- coding: utf-8 -*-
"""Unit contracts for tower elastic LATIN initialization."""

import unittest
import numpy as np

from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
)
from latin.tower_initialization import (
    compute_tower_elastic_initialization,
)
from latin.tower_state import (
    LatinStateTower,
    MaterialPointLayout,
)
from material.viscoplastic_damage_1d import MaterialParameters


class TestTowerElasticInitialization(unittest.TestCase):
    def setUp(self) -> None:
        layout = MaterialPointLayout(
            n_elements=1,
            n_gauss=1,
            n_fibers=1,
        )
        metric = MaterialPointMetric(
            weights=np.array([2.0], dtype=np.float64)
        )
        self.operator = TowerEquilibriumOperator(
            layout=layout,
            metric=metric,
            reference_modulus=np.array([100.0], dtype=np.float64),
            compatibility_matrix=np.array([[0.5]], dtype=np.float64),
            free_dofs=np.array([1], dtype=np.int64),
            n_dof=2,
        )
        self.material = MaterialParameters(E=100.0, h=0.25)
        self.time = np.array([0.0, 0.5, 1.0], dtype=np.float64)
        self.loads = np.array(
            [[0.0, 0.0], [0.0, 20.0], [0.0, -10.0]],
            dtype=np.float64,
        )

    def _run(self):
        return compute_tower_elastic_initialization(
            time=self.time,
            load_vectors=self.loads,
            materials=self.material,
            equilibrium_operator=self.operator,
            stress_to_force_factor=10.0,
        )

    def test_exact_reference_history(self) -> None:
        result = self._run()
        np.testing.assert_allclose(
            result.displacements,
            np.array(
                [[0.0, 0.0], [0.0, 0.04], [0.0, -0.02]]
            ),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            result.state.elastic_strain,
            np.array([[0.0], [0.02], [-0.01]]),
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            result.state.stress,
            np.array([[0.0], [2.0], [-1.0]]),
            rtol=0.0,
            atol=1.0e-14,
        )
        np.testing.assert_allclose(
            result.state.energy_release_rate,
            np.array([[0.0], [0.02], [0.00125]]),
            rtol=0.0,
            atol=1.0e-15,
        )

    def test_inelastic_fields_are_zero(self) -> None:
        state = self._run().state
        for name in (
            "plastic_strain_rate",
            "alpha_rate",
            "r_bar_rate",
            "damage_rate",
            "beta",
            "R_bar",
            "plastic_strain",
            "alpha",
            "r_bar",
            "damage",
        ):
            np.testing.assert_allclose(
                getattr(state, name),
                np.zeros(state.field_shape),
                rtol=0.0,
                atol=0.0,
            )

    def test_equilibrium_and_compatibility_hold(self) -> None:
        result = self._run()
        np.testing.assert_allclose(
            result.state.elastic_strain,
            result.displacements[:, self.operator.free_dofs]
            @ self.operator.compatibility_matrix.T,
            rtol=0.0,
            atol=1.0e-15,
        )
        np.testing.assert_allclose(
            result.free_equilibrium_residual,
            np.zeros((self.time.size, 1)),
            rtol=0.0,
            atol=1.0e-13,
        )

    def test_result_is_immutable(self) -> None:
        result = self._run()
        self.assertIsInstance(result.state, LatinStateTower)
        self.assertFalse(result.displacements.flags.writeable)
        self.assertFalse(result.load_vectors.flags.writeable)
        self.assertFalse(
            result.free_equilibrium_residual.flags.writeable
        )
        with self.assertRaises(ValueError):
            result.displacements[0, 0] = 1.0
        with self.assertRaises(ValueError):
            result.state.stress[0, 0] = 1.0

    def test_rejects_modulus_mismatch(self) -> None:
        with self.assertRaisesRegex(
            ValueError, "same reference elastic modulus"
        ):
            compute_tower_elastic_initialization(
                time=self.time,
                load_vectors=self.loads,
                materials=MaterialParameters(E=101.0),
                equilibrium_operator=self.operator,
                stress_to_force_factor=10.0,
            )

    def test_rejects_bad_force_conversion(self) -> None:
        for factor in (0.0, -1.0):
            with self.assertRaises(ValueError):
                compute_tower_elastic_initialization(
                    time=self.time,
                    load_vectors=self.loads,
                    materials=self.material,
                    equilibrium_operator=self.operator,
                    stress_to_force_factor=factor,
                )

    def test_rejects_wrong_load_shape(self) -> None:
        with self.assertRaisesRegex(
            ValueError, "load_vectors must have shape"
        ):
            compute_tower_elastic_initialization(
                time=self.time,
                load_vectors=np.zeros((3, 3)),
                materials=self.material,
                equilibrium_operator=self.operator,
                stress_to_force_factor=10.0,
            )


if __name__ == "__main__":
    unittest.main()

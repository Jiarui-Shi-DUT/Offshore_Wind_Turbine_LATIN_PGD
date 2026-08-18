# -*- coding: utf-8 -*-
"""Tests for the tower material-point metric and reference projection."""

import unittest
from types import SimpleNamespace

import numpy as np

from latin.tower_equilibrium_operator import (
    MaterialPointMetric,
    TowerEquilibriumOperator,
    build_tower_equilibrium_operator,
)
from latin.tower_state import MaterialPointLayout


class TestMaterialPointMetric(unittest.TestCase):
    """Verify geometric material-point metric operations."""

    def test_metric_operations(self) -> None:
        source_weights = np.array([1.0, 2.0, 3.0])
        metric = MaterialPointMetric(source_weights)

        source_weights[:] = 99.0
        np.testing.assert_allclose(
            metric.weights,
            np.array([1.0, 2.0, 3.0]),
        )
        self.assertFalse(metric.weights.flags.writeable)
        self.assertEqual(metric.n_material_points, 3)
        self.assertAlmostEqual(metric.total_measure, 6.0)

        values = np.array([2.0, -1.0, 4.0])
        np.testing.assert_allclose(
            metric.apply(values),
            np.array([2.0, -2.0, 12.0]),
        )
        self.assertAlmostEqual(
            metric.integrate(values),
            12.0,
        )

        history = np.vstack(
            (
                values,
                np.ones(3),
            )
        )
        np.testing.assert_allclose(
            metric.integrate(history),
            np.array([12.0, 6.0]),
        )

        other = np.array([1.0, 3.0, -2.0])
        expected_inner = float(
            np.sum(
                np.array([1.0, 2.0, 3.0])
                * values
                * other
            )
        )
        self.assertAlmostEqual(
            metric.inner_product(values, other),
            expected_inner,
        )
        self.assertAlmostEqual(
            metric.norm(values),
            np.sqrt(
                np.sum(
                    np.array([1.0, 2.0, 3.0])
                    * values**2
                )
            ),
        )

    def test_metric_rejects_nonpositive_weights(self) -> None:
        with self.assertRaises(ValueError):
            MaterialPointMetric(
                np.array([1.0, 0.0]),
            )
        with self.assertRaises(ValueError):
            MaterialPointMetric(
                np.array([1.0, -1.0]),
            )


class TestTowerEquilibriumOperator(unittest.TestCase):
    """Verify E = H (H^T M C0 H)^-1 H^T M C0."""

    def setUp(self) -> None:
        self.layout = MaterialPointLayout(
            n_elements=1,
            n_gauss=2,
            n_fibers=2,
        )
        self.metric = MaterialPointMetric(
            np.array([1.0, 2.0, 1.5, 0.5]),
        )
        self.h_matrix = np.array(
            [
                [1.0, 0.0],
                [0.0, 1.0],
                [1.0, 1.0],
                [2.0, -1.0],
            ],
            dtype=np.float64,
        )
        self.operator = TowerEquilibriumOperator(
            layout=self.layout,
            metric=self.metric,
            reference_modulus=np.array(
                [10.0, 20.0, 15.0, 25.0],
            ),
            compatibility_matrix=self.h_matrix,
            free_dofs=np.array([1, 3]),
            n_dof=4,
        )

    def test_compatible_source_is_fixed_point(self) -> None:
        displacement = np.array([0.25, -0.4])
        source = self.h_matrix @ displacement

        result = self.operator.apply_spatial(source)

        np.testing.assert_allclose(
            result.compatible_strain,
            source,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            result.stress,
            np.zeros(4),
            rtol=0.0,
            atol=1.0e-11,
        )
        np.testing.assert_allclose(
            result.displacement_free,
            displacement,
            rtol=1.0e-12,
            atol=1.0e-12,
        )

    def test_arbitrary_source_returns_equilibrated_stress(self) -> None:
        source = np.array([0.2, -0.1, 0.5, -0.3])
        result = self.operator.apply_spatial(source)

        residual = self.operator.equilibrium_residual(
            result.stress
        )
        np.testing.assert_allclose(
            residual,
            np.zeros(self.operator.n_free_dofs),
            rtol=0.0,
            atol=1.0e-11,
        )

        expected_stress = (
            self.operator.reference_modulus
            * (
                result.compatible_strain
                - source
            )
        )
        np.testing.assert_allclose(
            result.stress,
            expected_stress,
            rtol=0.0,
            atol=0.0,
        )

    def test_history_projection_matches_rowwise_spatial_projection(self) -> None:
        history = np.array(
            [
                [0.2, -0.1, 0.5, -0.3],
                [0.0, 0.1, -0.2, 0.4],
                [0.3, 0.2, 0.1, -0.1],
            ],
            dtype=np.float64,
        )

        history_result = self.operator.apply_history(history)

        expected_compatible = []
        expected_stress = []
        expected_displacement = []
        for row in history:
            spatial = self.operator.apply_spatial(row)
            expected_compatible.append(
                spatial.compatible_strain
            )
            expected_stress.append(spatial.stress)
            expected_displacement.append(
                spatial.displacement_free
            )

        np.testing.assert_allclose(
            history_result.compatible_strain,
            np.vstack(expected_compatible),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            history_result.stress,
            np.vstack(expected_stress),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            history_result.displacement_free,
            np.vstack(expected_displacement),
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            self.operator.equilibrium_residual(
                history_result.stress
            ),
            np.zeros(
                (
                    history.shape[0],
                    self.operator.n_free_dofs,
                )
            ),
            rtol=0.0,
            atol=1.0e-11,
        )


class TestTowerOperatorBuilder(unittest.TestCase):
    """Verify q ordering, metric construction, and H assembly from tower data."""

    def _fake_system(self):
        section_0 = SimpleNamespace(
            section=SimpleNamespace(
                n_fibers=2,
                areas=np.array([0.5, 0.5]),
                y_coordinates=np.array([-1.0, 1.0]),
            )
        )
        section_1 = SimpleNamespace(
            section=SimpleNamespace(
                n_fibers=2,
                areas=np.array([0.5, 0.5]),
                y_coordinates=np.array([-1.0, 1.0]),
            )
        )

        b_matrices = np.zeros((2, 2, 6))
        b_matrices[0, 0] = np.array(
            [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        )
        b_matrices[0, 1] = np.array(
            [0.0, 0.0, 0.0, 0.0, 1.0, 0.0]
        )
        b_matrices[1, 0] = np.array(
            [-1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        )
        b_matrices[1, 1] = np.array(
            [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]
        )

        element = SimpleNamespace(
            sections=(section_0, section_1),
            gauss_weights=np.array([1.0, 1.0]),
            b_matrices=b_matrices,
            transformation=np.eye(6),
            jacobian=2.0,
        )

        return SimpleNamespace(
            mesh=SimpleNamespace(
                n_elements=1,
                n_dof=6,
            ),
            elements=(element,),
            element_dofs=np.arange(6).reshape(1, 6),
            material=SimpleNamespace(E=10.0),
        )

    def test_builder_constructs_expected_metric_and_compatibility(self) -> None:
        operator = build_tower_equilibrium_operator(
            self._fake_system(),
        )

        self.assertEqual(
            operator.layout,
            MaterialPointLayout(1, 2, 2),
        )
        np.testing.assert_allclose(
            operator.metric.weights,
            np.ones(4),
            rtol=0.0,
            atol=0.0,
        )
        self.assertAlmostEqual(
            operator.metric.total_measure,
            4.0,
        )
        np.testing.assert_array_equal(
            operator.free_dofs,
            np.array([3, 4, 5]),
        )

        expected_h = np.array(
            [
                [1.0, 1.0, 0.0],
                [1.0, -1.0, 0.0],
                [1.0, 0.0, 1.0],
                [1.0, 0.0, -1.0],
            ]
        )
        np.testing.assert_allclose(
            operator.compatibility_matrix,
            expected_h,
            rtol=0.0,
            atol=0.0,
        )

        displacement = np.array([0.1, -0.2, 0.3])
        compatible = expected_h @ displacement
        projection = operator.apply_spatial(compatible)
        np.testing.assert_allclose(
            projection.compatible_strain,
            compatible,
            rtol=1.0e-12,
            atol=1.0e-12,
        )
        np.testing.assert_allclose(
            projection.stress,
            np.zeros(4),
            rtol=0.0,
            atol=1.0e-11,
        )


if __name__ == "__main__":
    unittest.main()

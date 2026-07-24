# -*- coding:utf-8 -*-
"""
作者：Shi Jiarui
日期：2026年07月24日
"""
import unittest

import numpy as np

from material.viscoplastic_damage_1d import (
    MaterialParameters,
    run_material_point,
)


class TestMaterialModel1D(unittest.TestCase):
    """一维黏塑性—损伤材料模型的回归测试。"""

    @classmethod
    def setUpClass(cls) -> None:
        material = MaterialParameters(sigma_y=80.0)

        cls.response = run_material_point(
            material=material,
            total_time=200.0,
            time_step=0.1,
        )

    def test_results_are_finite(self) -> None:
        """所有主要结果都必须是有限数值。"""
        self.assertTrue(np.all(np.isfinite(self.response.stress)))
        self.assertTrue(np.all(np.isfinite(self.response.damage)))
        self.assertTrue(
            np.all(
                np.isfinite(
                    self.response.energy_release_rate
                )
            )
        )

    def test_damage_is_admissible(self) -> None:
        """损伤应位于允许范围内，并保持单调不减。"""
        self.assertGreaterEqual(
            float(np.min(self.response.damage)),
            0.0,
        )
        self.assertLess(
            float(np.max(self.response.damage)),
            1.0,
        )
        self.assertTrue(
            np.all(
                np.diff(self.response.damage) >= -1.0e-12
            )
        )

    def test_reference_results(self) -> None:
        """重构后的结果应与已验证基准一致。"""
        self.assertAlmostEqual(
            float(self.response.damage[-1]),
            0.1869603995378383,
            places=10,
        )
        self.assertAlmostEqual(
            float(np.max(self.response.stress)),
            134.948339,
            places=5,
        )
        self.assertAlmostEqual(
            float(np.min(self.response.stress)),
            -135.677151,
            places=5,
        )


if __name__ == "__main__":
    unittest.main()
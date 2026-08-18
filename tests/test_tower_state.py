# -*- coding: utf-8 -*-
"""Unit tests for tower LATIN material-point layout and state ownership."""

import unittest

import numpy as np

from latin.tower_state import LatinStateTower, MaterialPointLayout


class TestMaterialPointLayout(unittest.TestCase):
    """Verify canonical q <-> (element, Gauss, fiber) topology."""

    def test_canonical_order_and_round_trip(self) -> None:
        layout = MaterialPointLayout(
            n_elements=2,
            n_gauss=2,
            n_fibers=3,
        )

        self.assertEqual(layout.n_material_points, 12)

        expected = [
            (0, 0, 0),
            (0, 0, 1),
            (0, 0, 2),
            (0, 1, 0),
            (0, 1, 1),
            (0, 1, 2),
            (1, 0, 0),
            (1, 0, 1),
            (1, 0, 2),
            (1, 1, 0),
            (1, 1, 1),
            (1, 1, 2),
        ]

        for q, triple in enumerate(expected):
            self.assertEqual(layout.unflatten(q), triple)
            self.assertEqual(layout.flatten(*triple), q)

        np.testing.assert_array_equal(
            layout.element_indices,
            np.array(
                [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1],
                dtype=np.int64,
            ),
        )
        np.testing.assert_array_equal(
            layout.gauss_indices,
            np.array(
                [0, 0, 0, 1, 1, 1, 0, 0, 0, 1, 1, 1],
                dtype=np.int64,
            ),
        )
        np.testing.assert_array_equal(
            layout.fiber_indices,
            np.array(
                [0, 1, 2, 0, 1, 2, 0, 1, 2, 0, 1, 2],
                dtype=np.int64,
            ),
        )

    def test_layout_rejects_invalid_sizes_and_indices(self) -> None:
        with self.assertRaises(TypeError):
            MaterialPointLayout(True, 2, 3)
        with self.assertRaises(ValueError):
            MaterialPointLayout(0, 2, 3)
        with self.assertRaises(IndexError):
            MaterialPointLayout(2, 2, 3).flatten(2, 0, 0)
        with self.assertRaises(IndexError):
            MaterialPointLayout(2, 2, 3).unflatten(12)


class TestLatinStateTower(unittest.TestCase):
    """Verify shape, validation, immutability, and no-alias ownership."""

    def test_zero_state_has_canonical_shape_and_read_only_arrays(self) -> None:
        time = np.array([0.2, 0.4, 0.7], dtype=np.float64)
        state = LatinStateTower.zeros(
            time=time,
            n_material_points=12,
        )

        self.assertEqual(state.n_time, 3)
        self.assertEqual(state.n_material_points, 12)
        self.assertEqual(state.field_shape, (3, 12))
        np.testing.assert_array_equal(state.time, time)

        self.assertFalse(state.time.flags.writeable)
        for field_name in LatinStateTower.MATERIAL_FIELD_NAMES:
            field = getattr(state, field_name)
            self.assertEqual(field.shape, (3, 12))
            self.assertTrue(np.all(field == 0.0))
            self.assertFalse(field.flags.writeable)

        with self.assertRaises(ValueError):
            state.stress[0, 0] = 1.0

    def test_constructor_owns_detached_copies(self) -> None:
        time = np.array([0.0, 1.0], dtype=np.float64)
        source = np.arange(6, dtype=np.float64).reshape(2, 3)
        damage = np.zeros((2, 3), dtype=np.float64)

        kwargs = {
            field_name: source
            for field_name in LatinStateTower.MATERIAL_FIELD_NAMES
        }
        kwargs["damage"] = damage

        state = LatinStateTower(
            time=time,
            **kwargs,
        )

        time[:] = 99.0
        source[:] = -99.0
        damage[:] = 0.9

        np.testing.assert_array_equal(
            state.time,
            np.array([0.0, 1.0], dtype=np.float64),
        )
        np.testing.assert_array_equal(
            state.stress,
            np.arange(6, dtype=np.float64).reshape(2, 3),
        )
        np.testing.assert_array_equal(
            state.damage,
            np.zeros((2, 3), dtype=np.float64),
        )

    def test_copy_is_detached_and_content_identical(self) -> None:
        original = LatinStateTower.zeros(
            time=np.array([0.0, 0.5, 1.0]),
            n_material_points=4,
        )
        copied = original.copy()

        self.assertIsNot(original, copied)
        np.testing.assert_array_equal(original.time, copied.time)

        for field_name in LatinStateTower.MATERIAL_FIELD_NAMES:
            original_field = getattr(original, field_name)
            copied_field = getattr(copied, field_name)
            np.testing.assert_array_equal(
                original_field,
                copied_field,
            )
            self.assertFalse(
                np.shares_memory(
                    original_field,
                    copied_field,
                )
            )

    def test_state_rejects_invalid_time_shape_and_damage(self) -> None:
        with self.assertRaises(ValueError):
            LatinStateTower.zeros(
                time=np.array([0.0, 0.0]),
                n_material_points=2,
            )

        valid = LatinStateTower.zeros(
            time=np.array([0.0, 1.0]),
            n_material_points=2,
        )

        fields = {
            name: np.array(getattr(valid, name), copy=True)
            for name in LatinStateTower.MATERIAL_FIELD_NAMES
        }
        fields["damage"][1, 1] = 1.0

        with self.assertRaises(ValueError):
            LatinStateTower(
                time=valid.time,
                **fields,
            )


if __name__ == "__main__":
    unittest.main()

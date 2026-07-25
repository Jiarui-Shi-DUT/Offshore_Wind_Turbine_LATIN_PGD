# -*- coding: utf-8 -*-
"""Regression tests for the LATIN-PGD saturation criterion."""

import unittest

from latin.pgd_saturation import (
    SaturationAction,
    decide_pgd_saturation,
    saturation_indicator,
)


class TestPGDSaturation(unittest.TestCase):
    """Verify Eq. (60) and the three adaptive decisions."""

    def test_saturation_indicator_formula(self) -> None:
        value = saturation_indicator(
            previous_indicator=0.2,
            current_indicator=0.1,
        )

        self.assertAlmostEqual(
            value,
            1.0 / 3.0,
            places=15,
        )

    def test_advance_latin_when_reduction_is_sufficient(self) -> None:
        decision = decide_pgd_saturation(
            previous_indicator=0.2,
            current_indicator=0.1,
        )

        self.assertEqual(
            decision.action,
            SaturationAction.ADVANCE_LATIN,
        )
        self.assertTrue(decision.should_advance_latin)
        self.assertFalse(decision.should_enrich_basis)
        self.assertFalse(decision.should_stop)

    def test_enrich_basis_when_reduction_is_insufficient(self) -> None:
        decision = decide_pgd_saturation(
            previous_indicator=0.2,
            current_indicator=0.18,
        )

        self.assertEqual(
            decision.action,
            SaturationAction.ENRICH_BASIS,
        )
        self.assertFalse(decision.should_advance_latin)
        self.assertTrue(decision.should_enrich_basis)
        self.assertFalse(decision.should_stop)
        self.assertGreater(
            decision.value,
            decision.stopping_tolerance,
        )
        self.assertLessEqual(
            decision.value,
            decision.enrichment_tolerance,
        )

    def test_stop_when_indicator_is_saturated(self) -> None:
        decision = decide_pgd_saturation(
            previous_indicator=0.2,
            current_indicator=0.19998,
        )

        self.assertEqual(
            decision.action,
            SaturationAction.STOP_SATURATED,
        )
        self.assertFalse(decision.should_advance_latin)
        self.assertFalse(decision.should_enrich_basis)
        self.assertTrue(decision.should_stop)
        self.assertLessEqual(
            decision.value,
            decision.stopping_tolerance,
        )

    def test_increasing_indicator_is_not_accepted(self) -> None:
        decision = decide_pgd_saturation(
            previous_indicator=0.1,
            current_indicator=0.12,
        )

        self.assertLess(decision.value, 0.0)
        self.assertEqual(
            decision.action,
            SaturationAction.STOP_SATURATED,
        )


if __name__ == "__main__":
    unittest.main()

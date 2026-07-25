# -*- coding: utf-8 -*-
"""
Saturation control for adaptive LATIN-PGD basis enrichment.

The paper defines the saturation parameter

    zeta = (xi_previous - xi_current)
           / (xi_previous + xi_current),

where xi is the relative LATIN indicator.  Three cases are distinguished:

1. zeta > zeta_enrich:
   updating the existing PGD basis has reduced the LATIN error sufficiently,
   so the algorithm proceeds to the next LATIN iteration;

2. zeta_stop < zeta <= zeta_enrich:
   the LATIN error reduction has become insufficient, so a new space-time
   PGD pair should be added;

3. zeta <= zeta_stop:
   the LATIN indicator is considered saturated and the outer algorithm may
   stop because further basis updates are no longer effective.

For the one-dimensional bar example, the paper uses zeta_enrich = 0.1 and
zeta_stop = 1.0e-4.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

import numpy as np


class SaturationAction(Enum):
    """Action selected from the PGD saturation parameter."""

    ADVANCE_LATIN = "advance_latin"
    ENRICH_BASIS = "enrich_basis"
    STOP_SATURATED = "stop_saturated"


@dataclass(frozen=True)
class PGDSaturationDecision:
    """Saturation value and the corresponding adaptive action."""

    value: float
    action: SaturationAction
    previous_indicator: float
    current_indicator: float
    enrichment_tolerance: float
    stopping_tolerance: float

    @property
    def should_advance_latin(self) -> bool:
        """Whether the current PGD approximation is satisfactory."""
        return self.action is SaturationAction.ADVANCE_LATIN

    @property
    def should_enrich_basis(self) -> bool:
        """Whether one additional space-time PGD pair is required."""
        return self.action is SaturationAction.ENRICH_BASIS

    @property
    def should_stop(self) -> bool:
        """Whether the LATIN indicator has numerically saturated."""
        return self.action is SaturationAction.STOP_SATURATED


def saturation_indicator(
    previous_indicator: float,
    current_indicator: float,
) -> float:
    """
    Evaluate the paper's saturation parameter from Eq. (60).

    Parameters
    ----------
    previous_indicator:
        Relative LATIN indicator xi_i before the current PGD update.
    current_indicator:
        Relative LATIN indicator xi_(i+1) after the current PGD update.

    Returns
    -------
    float
        zeta = (xi_i - xi_(i+1)) / (xi_i + xi_(i+1)).

    Notes
    -----
    A negative value is retained rather than clipped.  It indicates that the
    current LATIN indicator increased and therefore cannot be regarded as a
    satisfactory reduction.
    """
    previous = float(previous_indicator)
    current = float(current_indicator)

    if not np.isfinite(previous) or not np.isfinite(current):
        raise ValueError("LATIN indicators must be finite.")
    if previous < 0.0 or current < 0.0:
        raise ValueError("LATIN indicators must be non-negative.")

    denominator = previous + current
    numerical_zero = float(np.finfo(np.float64).eps)

    if denominator <= numerical_zero:
        return 0.0

    value = (previous - current) / denominator

    if not np.isfinite(value):
        raise FloatingPointError(
            "The PGD saturation parameter is non-finite."
        )

    return float(value)


def decide_pgd_saturation(
    previous_indicator: float,
    current_indicator: float,
    *,
    enrichment_tolerance: float = 0.1,
    stopping_tolerance: float = 1.0e-4,
) -> PGDSaturationDecision:
    """
    Select the adaptive PGD action associated with the saturation parameter.

    The stopping test has priority over enrichment.  Therefore a nearly zero
    improvement terminates the iteration instead of repeatedly adding modes.
    """
    enrichment = float(enrichment_tolerance)
    stopping = float(stopping_tolerance)

    if not np.isfinite(enrichment) or enrichment <= 0.0:
        raise ValueError(
            "enrichment_tolerance must be positive and finite."
        )
    if not np.isfinite(stopping) or stopping < 0.0:
        raise ValueError(
            "stopping_tolerance must be non-negative and finite."
        )
    if stopping >= enrichment:
        raise ValueError(
            "stopping_tolerance must be smaller than "
            "enrichment_tolerance."
        )

    value = saturation_indicator(
        previous_indicator=previous_indicator,
        current_indicator=current_indicator,
    )

    if value <= stopping:
        action = SaturationAction.STOP_SATURATED
    elif value <= enrichment:
        action = SaturationAction.ENRICH_BASIS
    else:
        action = SaturationAction.ADVANCE_LATIN

    return PGDSaturationDecision(
        value=value,
        action=action,
        previous_indicator=float(previous_indicator),
        current_indicator=float(current_indicator),
        enrichment_tolerance=enrichment,
        stopping_tolerance=stopping,
    )

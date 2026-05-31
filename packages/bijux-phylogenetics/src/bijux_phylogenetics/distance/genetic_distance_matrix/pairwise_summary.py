from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class _PairSummary:
    distance: float | None
    comparable_sites: int
    mismatch_sites: float
    transition_sites: float
    ag_transition_sites: float
    ct_transition_sites: float
    transversion_sites: float
    ambiguity_sites: int
    skipped_sites: int
    saturated: bool
    saturation_reason: str | None

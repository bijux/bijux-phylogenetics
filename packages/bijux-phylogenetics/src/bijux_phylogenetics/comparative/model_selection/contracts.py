from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class ComparativeModelComparisonRow:
    """Likelihood-comparison row for one comparative trait model."""

    model: str
    parameter_count: int
    log_likelihood: float
    aic: float
    aicc: float
    delta_aic: float | None = None
    delta_aicc: float | None = None
    rank: int | None = None
    comparable: bool = True
    comparability_note: str | None = None
    selected: bool = False
    likelihood_constant_policy: str | None = None
    akaike_weight: float | None = None
    within_delta_aic_threshold: bool | None = None
    within_delta_aicc_threshold: bool | None = None


@dataclass(slots=True)
class ComparativeModelComparisonReport:
    """Likelihood comparison between BM and OU trait models."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    rows: list[ComparativeModelComparisonRow]
    better_model: str
    warnings: list[str]

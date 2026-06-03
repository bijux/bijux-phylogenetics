from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.common import ComparativeReadinessReport
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeParameterInterval,
    ComparativeResidualSummary,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)


@dataclass(slots=True)
class BrownianRegimeExclusion:
    """One taxon excluded before multi-rate Brownian fitting."""

    taxon: str
    reason: str


@dataclass(slots=True)
class BrownianRegimeBranchRow:
    """One non-root branch assigned to a regime for multi-rate Brownian fitting."""

    branch_id: str
    regime: str
    branch_length: float
    descendant_taxa: list[str]
    analyzed_descendant_taxa: list[str]
    contributes_to_analysis: bool


@dataclass(slots=True)
class BrownianRegimeRateRow:
    """One regime-specific Brownian rate estimate."""

    regime: str
    branch_count: int
    contributing_branch_count: int
    total_branch_length: float
    contributing_branch_length: float
    sigma_squared: float
    lower_95: float
    upper_95: float
    interval_method: str


@dataclass(slots=True)
class BrownianRegimeProfileRow:
    """One conditional likelihood-profile row for one regime rate."""

    regime: str
    sigma_squared: float
    log_likelihood: float
    delta_log_likelihood: float
    in_support_interval: bool
    selected: bool


@dataclass(slots=True)
class BrownianRegimeIdentifiabilityWarning:
    """Warning that one regime-specific Brownian rate is weakly identified."""

    regime: str
    kind: str
    message: str


@dataclass(slots=True)
class BrownianRegimeFitSummaryReport:
    """Reviewer-facing multi-rate Brownian fit driven by a branch regime map."""

    tree_path: Path
    traits_path: Path
    regime_map_path: Path
    taxon_column: str
    branch_id_column: str
    regime_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[BrownianRegimeExclusion]
    branch_rows: list[BrownianRegimeBranchRow]
    regime_rows: list[BrownianRegimeRateRow]
    profile_rows: list[BrownianRegimeProfileRow]
    root_state: float
    root_state_interval: ComparativeParameterInterval
    log_likelihood: float
    aic: float
    aicc: float
    comparison_rows: list[ComparativeModelComparisonRow]
    better_model: str
    likelihood_ratio_statistic: float
    likelihood_ratio_degrees_of_freedom: int
    likelihood_ratio_p_value: float
    likelihood_ratio_p_value_method: str
    identifiability_warnings: list[BrownianRegimeIdentifiabilityWarning]
    residual_diagnostics: ComparativeResidualSummary
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport

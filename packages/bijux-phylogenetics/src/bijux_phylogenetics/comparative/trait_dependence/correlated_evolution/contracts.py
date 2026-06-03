from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class CorrelatedTraitExclusion:
    """One taxon excluded before pairwise trait-evolution analysis."""

    taxon: str
    reason: str
    missing_traits: list[str]


@dataclass(slots=True)
class CorrelatedTraitComparisonRow:
    """One model-comparison row for independent versus correlated evolution."""

    model_kind: str
    model_description: str
    parameter_count: int
    log_likelihood: float
    aic: float
    delta_aic: float
    selected: bool


@dataclass(slots=True)
class CorrelatedTraitObservationRow:
    """One reviewer-facing observation row used by a coupling analysis."""

    row_kind: str
    label: str
    taxon: str | None
    left_taxa: list[str]
    right_taxa: list[str]
    left_numeric_value: float | None
    right_numeric_value: float | None
    expected_variance: float | None
    left_state: str | None
    right_state: str | None
    joint_state: str | None


@dataclass(slots=True)
class CorrelatedTraitEvolutionReport:
    """Reviewer-facing pairwise trait-evolution coupling report."""

    tree_path: Path
    traits_path: Path
    left_trait: str
    right_trait: str
    taxon_column: str
    analysis_kind: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    excluded_taxa: list[CorrelatedTraitExclusion]
    observation_rows: list[CorrelatedTraitObservationRow]
    comparison_rows: list[CorrelatedTraitComparisonRow]
    association_measure_name: str
    association_measure_value: float
    evolutionary_covariance: float | None
    evolutionary_correlation: float | None
    lower_95_confidence_interval: float | None
    upper_95_confidence_interval: float | None
    independent_parameter_count: int
    independent_log_likelihood: float
    independent_aic: float
    correlated_parameter_count: int
    correlated_log_likelihood: float
    correlated_aic: float
    better_model: str
    likelihood_ratio_statistic: float
    likelihood_ratio_degrees_of_freedom: int
    likelihood_ratio_p_value: float
    likelihood_ratio_p_value_method: str
    left_root_estimate: float | None
    right_root_estimate: float | None
    left_state_order: list[str]
    right_state_order: list[str]
    joint_state_counts: dict[str, int]
    warnings: list[str]

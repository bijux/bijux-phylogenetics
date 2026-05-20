from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.assessment import (
    ComparativeMethodMaturityReport,
    ComparativeSensitivityReport,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeReadinessReport,
    NumericTraitSummary,
)
from bijux_phylogenetics.comparative.continuous import (
    BrownianMotionFitReport,
    OUTraitEvolutionSummaryReport,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonReport,
)
from bijux_phylogenetics.comparative.pgls import (
    ComparativeFormulaSpecification,
    PGLSInputReport,
    PGLSResult,
)
from bijux_phylogenetics.comparative.signal import (
    BlombergKReport,
    IndependentContrastReport,
    PagelLambdaReport,
)


@dataclass(slots=True)
class ComparativeAuditRow:
    """Auditable record for one comparative analysis surface."""

    analysis: str
    taxa_used: list[str]
    traits_used: list[str]
    excluded_taxa: list[str]
    assumptions: list[str]
    warnings: list[str]


@dataclass(slots=True)
class ComparativePredictorInfluenceRow:
    """One predictor term ranked by effect and test strength."""

    term: str
    estimate: float
    test_statistic: float
    p_value: float
    significant: bool


@dataclass(slots=True)
class ComparativeTaxonInfluenceRow:
    """Combined taxon influence summary across PGLS and leave-one-out sensitivity."""

    taxon: str
    leverage: float
    standardized_residual: float
    sensitivity_delta_log_likelihood: float
    sensitivity_delta_primary_parameter: float
    influence_score: float


@dataclass(slots=True)
class ComparativeInfluenceReport:
    """Predictor and taxon drivers for one comparative analysis."""

    tree_path: Path
    traits_path: Path
    response: str
    selected_model: str
    predictor_rows: list[ComparativePredictorInfluenceRow]
    taxon_rows: list[ComparativeTaxonInfluenceRow]
    top_predictor_terms: list[str]
    top_taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class ComparativeModelSnapshot:
    """Shared comparative snapshot used by reports and comparisons."""

    tree_path: Path
    traits_path: Path
    response: str
    formula: ComparativeFormulaSpecification
    readiness: ComparativeReadinessReport
    summary: NumericTraitSummary
    signal_k: BlombergKReport
    signal_lambda: PagelLambdaReport
    contrasts: IndependentContrastReport
    brownian: BrownianMotionFitReport
    ou: OUTraitEvolutionSummaryReport
    model_comparison: ComparativeModelComparisonReport
    pgls_inputs: PGLSInputReport
    pgls_model: PGLSResult
    sensitivity: ComparativeSensitivityReport
    maturity: ComparativeMethodMaturityReport
    audit_rows: list[ComparativeAuditRow]
    limitations: list[str]


@dataclass(slots=True)
class ComparativeMethodReport:
    """Integrated comparative-method report for one response trait."""

    snapshot: ComparativeModelSnapshot
    influence: ComparativeInfluenceReport


@dataclass(slots=True)
class ComparativeMethodsSummaryTextResult:
    """Reviewer-facing Markdown methods summary for one comparative analysis."""

    output_path: Path
    title: str
    selected_model: str
    predictor_count: int
    analysis_taxa: int
    excluded_taxa: int
    warning_count: int
    text: str
    report: ComparativeMethodReport


@dataclass(slots=True)
class ComparativeCoefficientDeltaRow:
    """Difference in one encoded coefficient between two comparative analyses."""

    term: str
    left_estimate: float | None
    right_estimate: float | None
    delta: float | None
    sign_changed: bool


@dataclass(slots=True)
class ComparativeTreeComparisonReport:
    """Compare comparative results across two alternative trees."""

    left_tree_path: Path
    right_tree_path: Path
    response: str
    predictors: list[str]
    left_selected_model: str
    right_selected_model: str
    delta_blombergs_k: float
    delta_pagels_lambda: float
    delta_brownian_rate: float
    delta_ou_alpha: float
    coefficient_deltas: list[ComparativeCoefficientDeltaRow]
    sign_changed_terms: list[str]
    conclusion_changed: bool
    warnings: list[str]


@dataclass(slots=True)
class ComparativePruningComparisonReport:
    """Compare comparative results before and after explicit taxon pruning."""

    tree_path: Path
    response: str
    predictors: list[str]
    baseline_taxa: list[str]
    pruned_taxa: list[str]
    dropped_taxa: list[str]
    delta_blombergs_k: float
    delta_pagels_lambda: float
    coefficient_deltas: list[ComparativeCoefficientDeltaRow]
    baseline_selected_model: str
    pruned_selected_model: str
    sign_changed_terms: list[str]
    conclusion_changed: bool
    warnings: list[str]


__all__ = [
    "ComparativeAuditRow",
    "ComparativeCoefficientDeltaRow",
    "ComparativeInfluenceReport",
    "ComparativeMethodReport",
    "ComparativeMethodsSummaryTextResult",
    "ComparativeModelSnapshot",
    "ComparativePredictorInfluenceRow",
    "ComparativePruningComparisonReport",
    "ComparativeTaxonInfluenceRow",
    "ComparativeTreeComparisonReport",
]

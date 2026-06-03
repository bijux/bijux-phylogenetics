from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, replace
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import stable_value, write_ancestral_rows
from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.ancestral.tree_set import (
    ContinuousAncestralTreeSetReport,
    DiscreteAncestralTreeSetReport,
)


@dataclass(frozen=True, slots=True)
class AncestralConfidenceSummary:
    """Reviewer-facing confidence summary for one ancestral evidence surface."""

    trait: str
    taxon_column: str
    source_kind: str
    reconstruction_kind: str
    target_kind: str
    model: str
    state_ordering: str | None
    alpha: float | None
    analyzed_taxon_count: int
    kept_tree_count: int | None
    confidence_row_count: int
    low_confidence_count: int
    unstable_count: int
    high_entropy_count: int
    top_uncertain_id: str | None
    top_uncertain_label: str | None
    top_uncertain_score: float | None
    warning_count: int


@dataclass(frozen=True, slots=True)
class ContinuousAncestralConfidenceRow:
    """One ranked internal-node confidence row for a continuous reconstruction."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    estimate: float
    standard_error: float
    lower_95_interval: float
    upper_95_interval: float
    uncertainty_width: float
    relative_uncertainty: float
    confidence: float
    uncertainty_score: float
    uncertainty_rank: int
    confidence_class: str
    unstable: bool


@dataclass(frozen=True, slots=True)
class DiscreteAncestralConfidenceRow:
    """One ranked internal-node confidence row for a discrete reconstruction."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_state: str
    state_set: list[str]
    state_probabilities: dict[str, float]
    max_posterior_probability: float
    runner_up_probability: float
    probability_margin: float
    entropy: float
    normalized_entropy: float
    uncertainty_score: float
    uncertainty_rank: int
    confidence_class: str
    ambiguous: bool
    unstable: bool


@dataclass(frozen=True, slots=True)
class ContinuousAncestralTreeSetConfidenceRow:
    """One ranked comparable-clade confidence row for a continuous tree set."""

    clade_id: str
    clade_taxa: list[str]
    tree_presence_count: int
    tree_presence_fraction: float
    mean_confidence: float
    mean_standard_error: float
    empirical_interval_width: float
    normalized_empirical_interval_width: float
    unstable_tree_count: int
    unstable_tree_fraction: float
    instability_score: float
    uncertainty_score: float
    uncertainty_rank: int
    confidence_class: str
    stability_class: str


@dataclass(frozen=True, slots=True)
class DiscreteAncestralTreeSetConfidenceRow:
    """One ranked comparable-clade confidence row for a discrete tree set."""

    clade_id: str
    clade_taxa: list[str]
    tree_presence_count: int
    tree_presence_fraction: float
    dominant_state: str
    dominant_state_fraction: float
    unique_state_count: int
    ambiguous_tree_fraction: float
    unstable_tree_fraction: float
    state_distribution: dict[str, int]
    entropy: float
    normalized_entropy: float
    instability_score: float
    uncertainty_score: float
    uncertainty_rank: int
    confidence_class: str
    stability_class: str


def build_continuous_ancestral_confidence_rows(
    report: ContinuousAncestralReport,
) -> list[ContinuousAncestralConfidenceRow]:
    """Build ranked internal-node confidence rows for one continuous reconstruction."""
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    if not internal_estimates:
        return []
    trait_range = _continuous_trait_range(report)
    rows = [
        ContinuousAncestralConfidenceRow(
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=estimate.descendant_taxa,
            estimate=estimate.estimate,
            standard_error=estimate.standard_error,
            lower_95_interval=estimate.lower_95_interval,
            upper_95_interval=estimate.upper_95_interval,
            uncertainty_width=estimate.uncertainty_width,
            relative_uncertainty=stable_value(
                estimate.uncertainty_width / _safe_scale(trait_range)
            ),
            confidence=estimate.confidence,
            uncertainty_score=stable_value(
                (1.0 - estimate.confidence)
                + min(estimate.uncertainty_width / _safe_scale(trait_range), 1.0)
            ),
            uncertainty_rank=0,
            confidence_class="",
            unstable=estimate.unstable,
        )
        for estimate in internal_estimates
    ]
    return _rank_continuous_rows(rows)


def summarize_continuous_ancestral_confidence(
    report: ContinuousAncestralReport,
) -> AncestralConfidenceSummary:
    """Summarize confidence signals for one continuous reconstruction."""
    rows = build_continuous_ancestral_confidence_rows(report)
    top_row = rows[0] if rows else None
    return AncestralConfidenceSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        source_kind="tree",
        reconstruction_kind="continuous",
        target_kind="internal_node",
        model=report.model,
        state_ordering=None,
        alpha=report.alpha,
        analyzed_taxon_count=report.taxon_count,
        kept_tree_count=None,
        confidence_row_count=len(rows),
        low_confidence_count=sum(row.confidence_class == "low" for row in rows),
        unstable_count=sum(row.unstable for row in rows),
        high_entropy_count=0,
        top_uncertain_id=top_row.node if top_row is not None else None,
        top_uncertain_label=top_row.node_name or top_row.node
        if top_row is not None
        else None,
        top_uncertain_score=top_row.uncertainty_score if top_row is not None else None,
        warning_count=len(report.warnings),
    )


def build_discrete_ancestral_confidence_rows(
    report: DiscreteAncestralReport,
) -> list[DiscreteAncestralConfidenceRow]:
    """Build ranked internal-node confidence rows for one discrete reconstruction."""
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    rows = [
        _build_discrete_confidence_row(
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=estimate.descendant_taxa,
            most_likely_state=estimate.most_likely_state,
            state_set=estimate.state_set,
            state_probabilities=estimate.state_probabilities,
            ambiguous=estimate.ambiguous,
            unstable=estimate.unstable,
        )
        for estimate in internal_estimates
    ]
    return _rank_discrete_rows(rows)


def summarize_discrete_ancestral_confidence(
    report: DiscreteAncestralReport,
) -> AncestralConfidenceSummary:
    """Summarize confidence signals for one discrete reconstruction."""
    rows = build_discrete_ancestral_confidence_rows(report)
    top_row = rows[0] if rows else None
    return AncestralConfidenceSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        source_kind="tree",
        reconstruction_kind="discrete",
        target_kind="internal_node",
        model=report.model,
        state_ordering=report.state_ordering,
        alpha=None,
        analyzed_taxon_count=report.taxon_count,
        kept_tree_count=None,
        confidence_row_count=len(rows),
        low_confidence_count=sum(row.confidence_class == "low" for row in rows),
        unstable_count=sum(row.unstable for row in rows),
        high_entropy_count=sum(row.normalized_entropy >= 0.5 for row in rows),
        top_uncertain_id=top_row.node if top_row is not None else None,
        top_uncertain_label=top_row.node_name or top_row.node
        if top_row is not None
        else None,
        top_uncertain_score=top_row.uncertainty_score if top_row is not None else None,
        warning_count=len(report.warnings),
    )


def build_continuous_ancestral_tree_set_confidence_rows(
    report: ContinuousAncestralTreeSetReport,
) -> list[ContinuousAncestralTreeSetConfidenceRow]:
    """Build ranked comparable-clade confidence rows for one continuous tree set."""
    confidence_by_clade: dict[str, list[float]] = defaultdict(list)
    for row in report.node_rows:
        confidence_by_clade[row.clade_id].append(row.confidence)
    maximum_interval_width = max(
        (row.empirical_interval_width for row in report.clade_summaries),
        default=0.0,
    )
    rows: list[ContinuousAncestralTreeSetConfidenceRow] = [
        ContinuousAncestralTreeSetConfidenceRow(
            clade_id=row.clade_id,
            clade_taxa=row.clade_taxa,
            tree_presence_count=row.tree_presence_count,
            tree_presence_fraction=row.tree_presence_fraction,
            mean_confidence=stable_value(
                sum(confidence_by_clade[row.clade_id])
                / len(confidence_by_clade[row.clade_id])
            ),
            mean_standard_error=row.mean_standard_error,
            empirical_interval_width=row.empirical_interval_width,
            normalized_empirical_interval_width=stable_value(
                row.empirical_interval_width / _safe_scale(maximum_interval_width)
            ),
            unstable_tree_count=row.unstable_tree_count,
            unstable_tree_fraction=row.unstable_tree_fraction,
            instability_score=row.instability_score,
            uncertainty_score=0.0,
            uncertainty_rank=0,
            confidence_class="",
            stability_class=row.stability_class,
        )
        for row in report.clade_summaries
    ]
    scored_rows = []
    for row in rows:
        uncertainty_score = stable_value(
            (1.0 - row.tree_presence_fraction)
            + (1.0 - row.mean_confidence)
            + row.unstable_tree_fraction
            + row.normalized_empirical_interval_width
        )
        scored_rows.append(
            ContinuousAncestralTreeSetConfidenceRow(
                clade_id=row.clade_id,
                clade_taxa=row.clade_taxa,
                tree_presence_count=row.tree_presence_count,
                tree_presence_fraction=row.tree_presence_fraction,
                mean_confidence=row.mean_confidence,
                mean_standard_error=row.mean_standard_error,
                empirical_interval_width=row.empirical_interval_width,
                normalized_empirical_interval_width=row.normalized_empirical_interval_width,
                unstable_tree_count=row.unstable_tree_count,
                unstable_tree_fraction=row.unstable_tree_fraction,
                instability_score=row.instability_score,
                uncertainty_score=uncertainty_score,
                uncertainty_rank=row.uncertainty_rank,
                confidence_class=row.confidence_class,
                stability_class=row.stability_class,
            )
        )
    return _rank_continuous_tree_set_rows(scored_rows)


def summarize_continuous_ancestral_tree_set_confidence(
    report: ContinuousAncestralTreeSetReport,
) -> AncestralConfidenceSummary:
    """Summarize confidence signals for one continuous ancestral tree-set review."""
    rows = build_continuous_ancestral_tree_set_confidence_rows(report)
    top_row = rows[0] if rows else None
    return AncestralConfidenceSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        source_kind="tree_set",
        reconstruction_kind="continuous",
        target_kind="comparable_clade",
        model=report.model,
        state_ordering=None,
        alpha=report.alpha,
        analyzed_taxon_count=len(report.analysis_taxa),
        kept_tree_count=report.kept_tree_count,
        confidence_row_count=len(rows),
        low_confidence_count=sum(row.confidence_class == "low" for row in rows),
        unstable_count=sum(row.stability_class != "stable" for row in rows),
        high_entropy_count=0,
        top_uncertain_id=top_row.clade_id if top_row is not None else None,
        top_uncertain_label=top_row.clade_id if top_row is not None else None,
        top_uncertain_score=top_row.uncertainty_score if top_row is not None else None,
        warning_count=len(report.warnings),
    )


def build_discrete_ancestral_tree_set_confidence_rows(
    report: DiscreteAncestralTreeSetReport,
) -> list[DiscreteAncestralTreeSetConfidenceRow]:
    """Build ranked comparable-clade confidence rows for one discrete tree set."""
    rows = []
    for row in report.clade_summaries:
        entropy, normalized_entropy = _distribution_entropy(
            row.state_distribution,
            total_count=row.tree_presence_count,
        )
        uncertainty_score = stable_value(
            (1.0 - row.tree_presence_fraction)
            + (1.0 - row.dominant_state_fraction)
            + row.ambiguous_tree_fraction
            + row.unstable_tree_fraction
        )
        rows.append(
            DiscreteAncestralTreeSetConfidenceRow(
                clade_id=row.clade_id,
                clade_taxa=row.clade_taxa,
                tree_presence_count=row.tree_presence_count,
                tree_presence_fraction=row.tree_presence_fraction,
                dominant_state=row.dominant_state,
                dominant_state_fraction=row.dominant_state_fraction,
                unique_state_count=row.unique_state_count,
                ambiguous_tree_fraction=row.ambiguous_tree_fraction,
                unstable_tree_fraction=row.unstable_tree_fraction,
                state_distribution=row.state_distribution,
                entropy=entropy,
                normalized_entropy=normalized_entropy,
                instability_score=row.instability_score,
                uncertainty_score=uncertainty_score,
                uncertainty_rank=0,
                confidence_class="",
                stability_class=row.stability_class,
            )
        )
    return _rank_discrete_tree_set_rows(rows)


def summarize_discrete_ancestral_tree_set_confidence(
    report: DiscreteAncestralTreeSetReport,
) -> AncestralConfidenceSummary:
    """Summarize confidence signals for one discrete ancestral tree-set review."""
    rows = build_discrete_ancestral_tree_set_confidence_rows(report)
    top_row = rows[0] if rows else None
    return AncestralConfidenceSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        source_kind="tree_set",
        reconstruction_kind="discrete",
        target_kind="comparable_clade",
        model=report.model,
        state_ordering=report.state_ordering,
        alpha=None,
        analyzed_taxon_count=len(report.analysis_taxa),
        kept_tree_count=report.kept_tree_count,
        confidence_row_count=len(rows),
        low_confidence_count=sum(row.confidence_class == "low" for row in rows),
        unstable_count=sum(row.stability_class != "stable" for row in rows),
        high_entropy_count=sum(row.normalized_entropy >= 0.5 for row in rows),
        top_uncertain_id=top_row.clade_id if top_row is not None else None,
        top_uncertain_label=top_row.clade_id if top_row is not None else None,
        top_uncertain_score=top_row.uncertainty_score if top_row is not None else None,
        warning_count=len(report.warnings),
    )


def write_ancestral_confidence_summary_table(
    path: Path,
    summary: AncestralConfidenceSummary,
) -> Path:
    """Write one summary ledger for an ancestral confidence review surface."""
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "source_kind",
            "reconstruction_kind",
            "target_kind",
            "model",
            "state_ordering",
            "alpha",
            "analyzed_taxon_count",
            "kept_tree_count",
            "confidence_row_count",
            "low_confidence_count",
            "unstable_count",
            "high_entropy_count",
            "top_uncertain_id",
            "top_uncertain_label",
            "top_uncertain_score",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "source_kind": summary.source_kind,
                "reconstruction_kind": summary.reconstruction_kind,
                "target_kind": summary.target_kind,
                "model": summary.model,
                "state_ordering": summary.state_ordering or "",
                "alpha": _format_optional_float(summary.alpha),
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "kept_tree_count": _format_optional_int(summary.kept_tree_count),
                "confidence_row_count": str(summary.confidence_row_count),
                "low_confidence_count": str(summary.low_confidence_count),
                "unstable_count": str(summary.unstable_count),
                "high_entropy_count": str(summary.high_entropy_count),
                "top_uncertain_id": summary.top_uncertain_id or "",
                "top_uncertain_label": summary.top_uncertain_label or "",
                "top_uncertain_score": _format_optional_float(
                    summary.top_uncertain_score
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_continuous_ancestral_confidence_table(
    path: Path,
    report: ContinuousAncestralReport,
) -> Path:
    """Write one ranked internal-node confidence ledger for a continuous reconstruction."""
    rows = build_continuous_ancestral_confidence_rows(report)
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "estimate",
            "standard_error",
            "lower_95_interval",
            "upper_95_interval",
            "uncertainty_width",
            "relative_uncertainty",
            "confidence",
            "uncertainty_score",
            "uncertainty_rank",
            "confidence_class",
            "unstable",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "estimate": str(row.estimate),
                "standard_error": str(row.standard_error),
                "lower_95_interval": str(row.lower_95_interval),
                "upper_95_interval": str(row.upper_95_interval),
                "uncertainty_width": str(row.uncertainty_width),
                "relative_uncertainty": str(row.relative_uncertainty),
                "confidence": str(row.confidence),
                "uncertainty_score": str(row.uncertainty_score),
                "uncertainty_rank": str(row.uncertainty_rank),
                "confidence_class": row.confidence_class,
                "unstable": str(row.unstable).lower(),
            }
            for row in rows
        ],
    )


def write_discrete_ancestral_confidence_table(
    path: Path,
    report: DiscreteAncestralReport,
) -> Path:
    """Write one ranked internal-node confidence ledger for a discrete reconstruction."""
    rows = build_discrete_ancestral_confidence_rows(report)
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_state",
            "state_set",
            "state_probabilities",
            "max_posterior_probability",
            "runner_up_probability",
            "probability_margin",
            "entropy",
            "normalized_entropy",
            "uncertainty_score",
            "uncertainty_rank",
            "confidence_class",
            "ambiguous",
            "unstable",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "most_likely_state": row.most_likely_state,
                "state_set": ",".join(row.state_set),
                "state_probabilities": json.dumps(
                    row.state_probabilities, sort_keys=True
                ),
                "max_posterior_probability": str(row.max_posterior_probability),
                "runner_up_probability": str(row.runner_up_probability),
                "probability_margin": str(row.probability_margin),
                "entropy": str(row.entropy),
                "normalized_entropy": str(row.normalized_entropy),
                "uncertainty_score": str(row.uncertainty_score),
                "uncertainty_rank": str(row.uncertainty_rank),
                "confidence_class": row.confidence_class,
                "ambiguous": str(row.ambiguous).lower(),
                "unstable": str(row.unstable).lower(),
            }
            for row in rows
        ],
    )


def write_continuous_ancestral_tree_set_confidence_table(
    path: Path,
    report: ContinuousAncestralTreeSetReport,
) -> Path:
    """Write one ranked comparable-clade confidence ledger for a continuous tree set."""
    rows = build_continuous_ancestral_tree_set_confidence_rows(report)
    return write_ancestral_rows(
        path,
        columns=[
            "clade_id",
            "clade_taxa",
            "tree_presence_count",
            "tree_presence_fraction",
            "mean_confidence",
            "mean_standard_error",
            "empirical_interval_width",
            "normalized_empirical_interval_width",
            "unstable_tree_count",
            "unstable_tree_fraction",
            "instability_score",
            "uncertainty_score",
            "uncertainty_rank",
            "confidence_class",
            "stability_class",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "mean_confidence": str(row.mean_confidence),
                "mean_standard_error": str(row.mean_standard_error),
                "empirical_interval_width": str(row.empirical_interval_width),
                "normalized_empirical_interval_width": str(
                    row.normalized_empirical_interval_width
                ),
                "unstable_tree_count": str(row.unstable_tree_count),
                "unstable_tree_fraction": str(row.unstable_tree_fraction),
                "instability_score": str(row.instability_score),
                "uncertainty_score": str(row.uncertainty_score),
                "uncertainty_rank": str(row.uncertainty_rank),
                "confidence_class": row.confidence_class,
                "stability_class": row.stability_class,
            }
            for row in rows
        ],
    )


def write_discrete_ancestral_tree_set_confidence_table(
    path: Path,
    report: DiscreteAncestralTreeSetReport,
) -> Path:
    """Write one ranked comparable-clade confidence ledger for a discrete tree set."""
    rows = build_discrete_ancestral_tree_set_confidence_rows(report)
    return write_ancestral_rows(
        path,
        columns=[
            "clade_id",
            "clade_taxa",
            "tree_presence_count",
            "tree_presence_fraction",
            "dominant_state",
            "dominant_state_fraction",
            "unique_state_count",
            "ambiguous_tree_fraction",
            "unstable_tree_fraction",
            "state_distribution",
            "entropy",
            "normalized_entropy",
            "instability_score",
            "uncertainty_score",
            "uncertainty_rank",
            "confidence_class",
            "stability_class",
        ],
        rows=[
            {
                "clade_id": row.clade_id,
                "clade_taxa": ",".join(row.clade_taxa),
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "dominant_state": row.dominant_state,
                "dominant_state_fraction": str(row.dominant_state_fraction),
                "unique_state_count": str(row.unique_state_count),
                "ambiguous_tree_fraction": str(row.ambiguous_tree_fraction),
                "unstable_tree_fraction": str(row.unstable_tree_fraction),
                "state_distribution": json.dumps(
                    row.state_distribution, sort_keys=True
                ),
                "entropy": str(row.entropy),
                "normalized_entropy": str(row.normalized_entropy),
                "instability_score": str(row.instability_score),
                "uncertainty_score": str(row.uncertainty_score),
                "uncertainty_rank": str(row.uncertainty_rank),
                "confidence_class": row.confidence_class,
                "stability_class": row.stability_class,
            }
            for row in rows
        ],
    )


def _rank_continuous_rows(
    rows: list[ContinuousAncestralConfidenceRow],
) -> list[ContinuousAncestralConfidenceRow]:
    ranked = sorted(
        rows,
        key=lambda row: (-row.uncertainty_score, row.node),
    )
    return [
        replace(
            row,
            uncertainty_rank=index,
            confidence_class=_score_to_confidence_class(row.uncertainty_score),
        )
        for index, row in enumerate(ranked, start=1)
    ]


def _rank_discrete_rows(
    rows: list[DiscreteAncestralConfidenceRow],
) -> list[DiscreteAncestralConfidenceRow]:
    ranked = sorted(
        rows,
        key=lambda row: (-row.uncertainty_score, row.node),
    )
    return [
        replace(
            row,
            uncertainty_rank=index,
            confidence_class=_score_to_confidence_class(row.uncertainty_score),
        )
        for index, row in enumerate(ranked, start=1)
    ]


def _rank_continuous_tree_set_rows(
    rows: list[ContinuousAncestralTreeSetConfidenceRow],
) -> list[ContinuousAncestralTreeSetConfidenceRow]:
    ranked = sorted(
        rows,
        key=lambda row: (-row.uncertainty_score, row.clade_id),
    )
    return [
        replace(
            row,
            uncertainty_rank=index,
            confidence_class=_score_to_confidence_class(row.uncertainty_score),
        )
        for index, row in enumerate(ranked, start=1)
    ]


def _rank_discrete_tree_set_rows(
    rows: list[DiscreteAncestralTreeSetConfidenceRow],
) -> list[DiscreteAncestralTreeSetConfidenceRow]:
    ranked = sorted(
        rows,
        key=lambda row: (-row.uncertainty_score, row.clade_id),
    )
    return [
        replace(
            row,
            uncertainty_rank=index,
            confidence_class=_score_to_confidence_class(row.uncertainty_score),
        )
        for index, row in enumerate(ranked, start=1)
    ]


def _build_discrete_confidence_row(
    *,
    node: str,
    node_name: str | None,
    descendant_taxa: list[str],
    most_likely_state: str,
    state_set: list[str],
    state_probabilities: dict[str, float],
    ambiguous: bool,
    unstable: bool,
) -> DiscreteAncestralConfidenceRow:
    ordered_probabilities = sorted(state_probabilities.values(), reverse=True)
    max_probability = ordered_probabilities[0] if ordered_probabilities else 0.0
    runner_up = ordered_probabilities[1] if len(ordered_probabilities) > 1 else 0.0
    entropy, normalized_entropy = _distribution_entropy(
        {
            state: 1
            for state, probability in state_probabilities.items()
            if probability > 0.0
        },
        probability_distribution=state_probabilities,
    )
    uncertainty_score = stable_value((1.0 - max_probability) + normalized_entropy)
    return DiscreteAncestralConfidenceRow(
        node=node,
        node_name=node_name,
        descendant_taxa=descendant_taxa,
        most_likely_state=most_likely_state,
        state_set=state_set,
        state_probabilities=state_probabilities,
        max_posterior_probability=stable_value(max_probability),
        runner_up_probability=stable_value(runner_up),
        probability_margin=stable_value(max_probability - runner_up),
        entropy=entropy,
        normalized_entropy=normalized_entropy,
        uncertainty_score=uncertainty_score,
        uncertainty_rank=0,
        confidence_class="",
        ambiguous=ambiguous,
        unstable=unstable,
    )


def _distribution_entropy(
    distribution: dict[str, int],
    *,
    total_count: int | None = None,
    probability_distribution: dict[str, float] | None = None,
) -> tuple[float, float]:
    if probability_distribution is not None:
        probabilities = [
            value for value in probability_distribution.values() if value > 0.0
        ]
    else:
        resolved_total = (
            total_count if total_count is not None else sum(distribution.values())
        )
        if resolved_total <= 0:
            return 0.0, 0.0
        probabilities = [
            count / resolved_total for count in distribution.values() if count > 0
        ]
    if len(probabilities) <= 1:
        return 0.0, 0.0
    entropy = -sum(
        probability * math.log2(probability) for probability in probabilities
    )
    maximum_entropy = math.log2(len(probabilities))
    normalized_entropy = entropy / maximum_entropy if maximum_entropy > 0.0 else 0.0
    return stable_value(entropy), stable_value(normalized_entropy)


def _continuous_trait_range(report: ContinuousAncestralReport) -> float:
    tip_values = [estimate.estimate for estimate in report.estimates if estimate.is_tip]
    if not tip_values:
        return 0.0
    return max(tip_values) - min(tip_values)


def _safe_scale(value: float) -> float:
    return max(value, 1e-12)


def _score_to_confidence_class(score: float) -> str:
    if score <= 0.4:
        return "high"
    if score <= 0.9:
        return "moderate"
    return "low"


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)


def _format_optional_int(value: int | None) -> str:
    if value is None:
        return ""
    return str(value)

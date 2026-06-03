from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)


@dataclass(frozen=True, slots=True)
class OrderedDiscreteFitRow:
    """One fitted discrete ancestral model under one state-ordering assumption."""

    ordering_mode: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    analyzed_taxon_count: int
    log_likelihood: float
    parameter_count: int
    aic: float
    root_most_likely_state: str
    root_confidence: float


@dataclass(frozen=True, slots=True)
class OrderedDiscreteNodeComparisonRow:
    """One node-wise comparison between ordered and unordered reconstructions."""

    node: str
    descendant_taxa: list[str]
    ordered_state: str
    unordered_state: str
    ordered_confidence: float
    unordered_confidence: float
    confidence_delta: float
    differs: bool
    ambiguity_changed: bool


@dataclass(frozen=True, slots=True)
class OrderedDiscreteTransitionComparisonRow:
    """One directed transition under ordered and unordered fitted models."""

    source_state: str
    target_state: str
    step_distance: int
    ordered_transition_allowed: bool
    unordered_transition_allowed: bool
    ordered_rate: float
    unordered_rate: float


@dataclass(frozen=True, slots=True)
class OrderedDiscreteSummary:
    """Reviewer-facing summary for ordered discrete ancestral reconstruction."""

    trait: str
    taxon_column: str
    model: str
    analyzed_taxon_count: int
    state_count: int
    ordered_log_likelihood: float
    unordered_log_likelihood: float
    ordered_parameter_count: int
    unordered_parameter_count: int
    ordered_aic: float
    unordered_aic: float
    delta_aic: float
    preferred_ordering: str
    differing_node_count: int
    ambiguity_change_count: int
    restricted_transition_count: int
    warning_count: int


@dataclass(slots=True)
class OrderedDiscreteReport:
    """Comparison between ordered and unordered discrete ancestral likelihood fits."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    ordered_states: list[str]
    fit_rows: list[OrderedDiscreteFitRow]
    node_rows: list[OrderedDiscreteNodeComparisonRow]
    transition_rows: list[OrderedDiscreteTransitionComparisonRow]
    warnings: list[str]


def summarize_ordered_discrete_reconstruction(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    ordered_states: list[str],
    taxon_column: str | None = None,
    model: str = "equal-rates",
) -> OrderedDiscreteReport:
    """Compare one ordered discrete likelihood fit against the unordered baseline."""
    if model == "meristic":
        _resolve_discrete_model_name(model)
    if model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ValueError(
            "ordered discrete ancestral reconstruction requires a discrete likelihood model"
        )
    ordered_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering="ordered",
        ordered_states=ordered_states,
    )
    unordered_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering="unordered",
    )
    warnings = list(ordered_report.warnings)
    ordered_fit = _build_fit_row(ordered_report, ordering_mode="ordered")
    unordered_fit = _build_fit_row(unordered_report, ordering_mode="unordered")
    node_rows = _build_node_rows(ordered_report, unordered_report)
    transition_rows = _build_transition_rows(ordered_report, unordered_report)
    if any(row.differs for row in node_rows):
        warnings.append(
            "one or more internal-node ancestral conclusions differ between ordered and unordered state-ordering assumptions"
        )
    return OrderedDiscreteReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=ordered_report.taxon_column,
        model=model,
        ordered_states=list(ordered_states),
        fit_rows=[ordered_fit, unordered_fit],
        node_rows=node_rows,
        transition_rows=transition_rows,
        warnings=warnings,
    )


def summarize_ordered_discrete_report(
    report: OrderedDiscreteReport,
) -> OrderedDiscreteSummary:
    """Summarize the main review facts for ordered discrete reconstruction."""
    ordered_fit = next(row for row in report.fit_rows if row.ordering_mode == "ordered")
    unordered_fit = next(
        row for row in report.fit_rows if row.ordering_mode == "unordered"
    )
    delta_aic = ordered_fit.aic - unordered_fit.aic
    if abs(delta_aic) < 1e-9:
        preferred_ordering = "tie"
    elif delta_aic < 0.0:
        preferred_ordering = "ordered"
    else:
        preferred_ordering = "unordered"
    return OrderedDiscreteSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        analyzed_taxon_count=ordered_fit.analyzed_taxon_count,
        state_count=len(report.ordered_states),
        ordered_log_likelihood=ordered_fit.log_likelihood,
        unordered_log_likelihood=unordered_fit.log_likelihood,
        ordered_parameter_count=ordered_fit.parameter_count,
        unordered_parameter_count=unordered_fit.parameter_count,
        ordered_aic=ordered_fit.aic,
        unordered_aic=unordered_fit.aic,
        delta_aic=delta_aic,
        preferred_ordering=preferred_ordering,
        differing_node_count=sum(row.differs for row in report.node_rows),
        ambiguity_change_count=sum(row.ambiguity_changed for row in report.node_rows),
        restricted_transition_count=sum(
            not row.ordered_transition_allowed for row in report.transition_rows
        ),
        warning_count=len(report.warnings),
    )


def write_ordered_discrete_summary_table(
    path: Path,
    report: OrderedDiscreteReport,
) -> Path:
    """Write one overall summary ledger for ordered discrete reconstruction."""
    summary = summarize_ordered_discrete_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "analyzed_taxon_count",
            "state_count",
            "ordered_log_likelihood",
            "unordered_log_likelihood",
            "ordered_parameter_count",
            "unordered_parameter_count",
            "ordered_aic",
            "unordered_aic",
            "delta_aic",
            "preferred_ordering",
            "differing_node_count",
            "ambiguity_change_count",
            "restricted_transition_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "state_count": str(summary.state_count),
                "ordered_log_likelihood": str(summary.ordered_log_likelihood),
                "unordered_log_likelihood": str(summary.unordered_log_likelihood),
                "ordered_parameter_count": str(summary.ordered_parameter_count),
                "unordered_parameter_count": str(summary.unordered_parameter_count),
                "ordered_aic": str(summary.ordered_aic),
                "unordered_aic": str(summary.unordered_aic),
                "delta_aic": str(summary.delta_aic),
                "preferred_ordering": summary.preferred_ordering,
                "differing_node_count": str(summary.differing_node_count),
                "ambiguity_change_count": str(summary.ambiguity_change_count),
                "restricted_transition_count": str(summary.restricted_transition_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_ordered_discrete_fit_table(
    path: Path,
    report: OrderedDiscreteReport,
) -> Path:
    """Write one per-fit ledger for ordered discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "ordering_mode",
            "model",
            "state_ordering",
            "ordered_states",
            "analyzed_taxon_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "root_most_likely_state",
            "root_confidence",
        ],
        rows=[
            {
                "ordering_mode": row.ordering_mode,
                "model": row.model,
                "state_ordering": row.state_ordering,
                "ordered_states": ",".join(row.ordered_states),
                "analyzed_taxon_count": str(row.analyzed_taxon_count),
                "log_likelihood": str(row.log_likelihood),
                "parameter_count": str(row.parameter_count),
                "aic": str(row.aic),
                "root_most_likely_state": row.root_most_likely_state,
                "root_confidence": str(row.root_confidence),
            }
            for row in report.fit_rows
        ],
    )


def write_ordered_discrete_node_table(
    path: Path,
    report: OrderedDiscreteReport,
) -> Path:
    """Write one node-wise comparison ledger for ordered discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "ordered_state",
            "unordered_state",
            "ordered_confidence",
            "unordered_confidence",
            "confidence_delta",
            "differs",
            "ambiguity_changed",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "ordered_state": row.ordered_state,
                "unordered_state": row.unordered_state,
                "ordered_confidence": str(row.ordered_confidence),
                "unordered_confidence": str(row.unordered_confidence),
                "confidence_delta": str(row.confidence_delta),
                "differs": str(row.differs).lower(),
                "ambiguity_changed": str(row.ambiguity_changed).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_ordered_discrete_transition_table(
    path: Path,
    report: OrderedDiscreteReport,
) -> Path:
    """Write one directed-transition ledger for ordered discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_state",
            "target_state",
            "step_distance",
            "ordered_transition_allowed",
            "unordered_transition_allowed",
            "ordered_rate",
            "unordered_rate",
        ],
        rows=[
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "step_distance": str(row.step_distance),
                "ordered_transition_allowed": str(
                    row.ordered_transition_allowed
                ).lower(),
                "unordered_transition_allowed": str(
                    row.unordered_transition_allowed
                ).lower(),
                "ordered_rate": str(row.ordered_rate),
                "unordered_rate": str(row.unordered_rate),
            }
            for row in report.transition_rows
        ],
    )


def _build_fit_row(report, *, ordering_mode: str) -> OrderedDiscreteFitRow:
    root_estimate = max(
        (estimate for estimate in report.estimates if not estimate.is_tip),
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    )
    return OrderedDiscreteFitRow(
        ordering_mode=ordering_mode,
        model=report.model,
        state_ordering=report.state_ordering,
        ordered_states=list(report.ordered_states),
        analyzed_taxon_count=report.taxon_count,
        log_likelihood=_required_fit_value(report.log_likelihood),
        parameter_count=_required_parameter_count(report.parameter_count),
        aic=_required_fit_value(report.aic),
        root_most_likely_state=root_estimate.most_likely_state,
        root_confidence=root_estimate.confidence,
    )


def _build_node_rows(
    ordered_report, unordered_report
) -> list[OrderedDiscreteNodeComparisonRow]:
    unordered_by_node = {
        estimate.node: estimate
        for estimate in unordered_report.estimates
        if not estimate.is_tip
    }
    rows: list[OrderedDiscreteNodeComparisonRow] = []
    for ordered_estimate in ordered_report.estimates:
        if ordered_estimate.is_tip or ordered_estimate.node not in unordered_by_node:
            continue
        unordered_estimate = unordered_by_node[ordered_estimate.node]
        rows.append(
            OrderedDiscreteNodeComparisonRow(
                node=ordered_estimate.node,
                descendant_taxa=ordered_estimate.descendant_taxa,
                ordered_state=ordered_estimate.most_likely_state,
                unordered_state=unordered_estimate.most_likely_state,
                ordered_confidence=ordered_estimate.confidence,
                unordered_confidence=unordered_estimate.confidence,
                confidence_delta=ordered_estimate.confidence
                - unordered_estimate.confidence,
                differs=ordered_estimate.most_likely_state
                != unordered_estimate.most_likely_state,
                ambiguity_changed=ordered_estimate.ambiguous
                != unordered_estimate.ambiguous,
            )
        )
    return rows


def _build_transition_rows(
    ordered_report,
    unordered_report,
) -> list[OrderedDiscreteTransitionComparisonRow]:
    ordered_by_pair = {
        (row.source_state, row.target_state): row
        for row in ordered_report.transition_rate_rows
    }
    unordered_by_pair = {
        (row.source_state, row.target_state): row
        for row in unordered_report.transition_rate_rows
    }
    rows: list[OrderedDiscreteTransitionComparisonRow] = []
    for pair in sorted(ordered_by_pair):
        ordered_row = ordered_by_pair[pair]
        unordered_row = unordered_by_pair[pair]
        rows.append(
            OrderedDiscreteTransitionComparisonRow(
                source_state=ordered_row.source_state,
                target_state=ordered_row.target_state,
                step_distance=ordered_row.step_distance,
                ordered_transition_allowed=ordered_row.transition_allowed,
                unordered_transition_allowed=unordered_row.transition_allowed,
                ordered_rate=ordered_row.rate,
                unordered_rate=unordered_row.rate,
            )
        )
    return rows


def _required_fit_value(value: float | None) -> float:
    if value is None:
        raise ValueError(
            "ordered discrete reconstruction requires likelihood fit details"
        )
    return value


def _required_parameter_count(value: int | None) -> int:
    if value is None:
        raise ValueError(
            "ordered discrete reconstruction requires likelihood fit details"
        )
    return value

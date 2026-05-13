from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.discrete import (
    reconstruct_discrete_ancestral_states,
)


@dataclass(frozen=True, slots=True)
class IrreversibleDiscreteFitRow:
    """One fitted discrete ancestral model under one transition-constraint regime."""

    constraint_mode: str
    model: str
    analyzed_taxon_count: int
    log_likelihood: float
    parameter_count: int
    aic: float
    root_most_likely_state: str
    root_confidence: float


@dataclass(frozen=True, slots=True)
class IrreversibleDiscreteNodeComparisonRow:
    """One node-wise comparison between constrained and unconstrained fits."""

    node: str
    descendant_taxa: list[str]
    constrained_state: str
    unconstrained_state: str
    constrained_confidence: float
    unconstrained_confidence: float
    confidence_delta: float
    differs: bool
    ambiguity_changed: bool


@dataclass(frozen=True, slots=True)
class IrreversibleDiscreteTransitionComparisonRow:
    """One directed transition under constrained and unconstrained fits."""

    source_state: str
    target_state: str
    constrained_transition_allowed: bool
    unconstrained_transition_allowed: bool
    constrained_rate: float
    unconstrained_rate: float


@dataclass(frozen=True, slots=True)
class IrreversibleDiscreteSummary:
    """Reviewer-facing summary for irreversible discrete reconstruction."""

    trait: str
    taxon_column: str
    model: str
    analyzed_taxon_count: int
    constrained_log_likelihood: float
    unconstrained_log_likelihood: float
    likelihood_difference: float
    constrained_parameter_count: int
    unconstrained_parameter_count: int
    constrained_aic: float
    unconstrained_aic: float
    delta_aic: float
    preferred_constraint: str
    differing_node_count: int
    ambiguity_change_count: int
    forbidden_transition_count: int
    warning_count: int


@dataclass(slots=True)
class IrreversibleDiscreteReport:
    """Comparison between constrained and unconstrained discrete ancestral fits."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    allowed_transition_pairs: list[tuple[str, str]]
    fit_rows: list[IrreversibleDiscreteFitRow]
    node_rows: list[IrreversibleDiscreteNodeComparisonRow]
    transition_rows: list[IrreversibleDiscreteTransitionComparisonRow]
    warnings: list[str]


def summarize_irreversible_discrete_reconstruction(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    allowed_transition_pairs: list[tuple[str, str]],
    taxon_column: str | None = None,
    model: str = "all-rates-different",
) -> IrreversibleDiscreteReport:
    """Compare constrained and unconstrained discrete ancestral likelihood fits."""
    if model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ValueError(
            "irreversible discrete ancestral reconstruction requires a discrete likelihood model"
        )
    constrained_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    unconstrained_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    warnings = list(constrained_report.warnings)
    node_rows = _build_node_rows(constrained_report, unconstrained_report)
    transition_rows = _build_transition_rows(constrained_report, unconstrained_report)
    if any(row.differs for row in node_rows):
        warnings.append(
            "one or more internal-node ancestral conclusions differ between constrained and unconstrained transition graphs"
        )
    return IrreversibleDiscreteReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=constrained_report.taxon_column,
        model=model,
        allowed_transition_pairs=list(allowed_transition_pairs),
        fit_rows=[
            _build_fit_row(constrained_report, constraint_mode="constrained"),
            _build_fit_row(unconstrained_report, constraint_mode="unconstrained"),
        ],
        node_rows=node_rows,
        transition_rows=transition_rows,
        warnings=warnings,
    )


def summarize_irreversible_discrete_report(
    report: IrreversibleDiscreteReport,
) -> IrreversibleDiscreteSummary:
    """Summarize the main review facts for irreversible discrete reconstruction."""
    constrained_fit = next(
        row for row in report.fit_rows if row.constraint_mode == "constrained"
    )
    unconstrained_fit = next(
        row for row in report.fit_rows if row.constraint_mode == "unconstrained"
    )
    delta_aic = constrained_fit.aic - unconstrained_fit.aic
    if abs(delta_aic) < 1e-9:
        preferred_constraint = "tie"
    elif delta_aic < 0.0:
        preferred_constraint = "constrained"
    else:
        preferred_constraint = "unconstrained"
    return IrreversibleDiscreteSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        analyzed_taxon_count=constrained_fit.analyzed_taxon_count,
        constrained_log_likelihood=constrained_fit.log_likelihood,
        unconstrained_log_likelihood=unconstrained_fit.log_likelihood,
        likelihood_difference=(
            constrained_fit.log_likelihood - unconstrained_fit.log_likelihood
        ),
        constrained_parameter_count=constrained_fit.parameter_count,
        unconstrained_parameter_count=unconstrained_fit.parameter_count,
        constrained_aic=constrained_fit.aic,
        unconstrained_aic=unconstrained_fit.aic,
        delta_aic=delta_aic,
        preferred_constraint=preferred_constraint,
        differing_node_count=sum(row.differs for row in report.node_rows),
        ambiguity_change_count=sum(row.ambiguity_changed for row in report.node_rows),
        forbidden_transition_count=sum(
            not row.constrained_transition_allowed for row in report.transition_rows
        ),
        warning_count=len(report.warnings),
    )


def write_irreversible_discrete_summary_table(
    path: Path,
    report: IrreversibleDiscreteReport,
) -> Path:
    """Write one overall summary ledger for irreversible discrete reconstruction."""
    summary = summarize_irreversible_discrete_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "analyzed_taxon_count",
            "constrained_log_likelihood",
            "unconstrained_log_likelihood",
            "likelihood_difference",
            "constrained_parameter_count",
            "unconstrained_parameter_count",
            "constrained_aic",
            "unconstrained_aic",
            "delta_aic",
            "preferred_constraint",
            "differing_node_count",
            "ambiguity_change_count",
            "forbidden_transition_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "constrained_log_likelihood": str(summary.constrained_log_likelihood),
                "unconstrained_log_likelihood": str(
                    summary.unconstrained_log_likelihood
                ),
                "likelihood_difference": str(summary.likelihood_difference),
                "constrained_parameter_count": str(summary.constrained_parameter_count),
                "unconstrained_parameter_count": str(
                    summary.unconstrained_parameter_count
                ),
                "constrained_aic": str(summary.constrained_aic),
                "unconstrained_aic": str(summary.unconstrained_aic),
                "delta_aic": str(summary.delta_aic),
                "preferred_constraint": summary.preferred_constraint,
                "differing_node_count": str(summary.differing_node_count),
                "ambiguity_change_count": str(summary.ambiguity_change_count),
                "forbidden_transition_count": str(summary.forbidden_transition_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_irreversible_discrete_fit_table(
    path: Path,
    report: IrreversibleDiscreteReport,
) -> Path:
    """Write one per-fit ledger for irreversible discrete reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "constraint_mode",
            "model",
            "analyzed_taxon_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "root_most_likely_state",
            "root_confidence",
        ],
        rows=[
            {
                "constraint_mode": row.constraint_mode,
                "model": row.model,
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


def write_irreversible_discrete_node_table(
    path: Path,
    report: IrreversibleDiscreteReport,
) -> Path:
    """Write one node-wise comparison ledger for irreversible reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "constrained_state",
            "unconstrained_state",
            "constrained_confidence",
            "unconstrained_confidence",
            "confidence_delta",
            "differs",
            "ambiguity_changed",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "constrained_state": row.constrained_state,
                "unconstrained_state": row.unconstrained_state,
                "constrained_confidence": str(row.constrained_confidence),
                "unconstrained_confidence": str(row.unconstrained_confidence),
                "confidence_delta": str(row.confidence_delta),
                "differs": str(row.differs).lower(),
                "ambiguity_changed": str(row.ambiguity_changed).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_irreversible_discrete_transition_table(
    path: Path,
    report: IrreversibleDiscreteReport,
) -> Path:
    """Write one directed-transition ledger for irreversible reconstruction."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_state",
            "target_state",
            "constrained_transition_allowed",
            "unconstrained_transition_allowed",
            "constrained_rate",
            "unconstrained_rate",
        ],
        rows=[
            {
                "source_state": row.source_state,
                "target_state": row.target_state,
                "constrained_transition_allowed": str(
                    row.constrained_transition_allowed
                ).lower(),
                "unconstrained_transition_allowed": str(
                    row.unconstrained_transition_allowed
                ).lower(),
                "constrained_rate": str(row.constrained_rate),
                "unconstrained_rate": str(row.unconstrained_rate),
            }
            for row in report.transition_rows
        ],
    )


def _build_fit_row(report, *, constraint_mode: str) -> IrreversibleDiscreteFitRow:
    root_estimate = max(
        (estimate for estimate in report.estimates if not estimate.is_tip),
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    )
    return IrreversibleDiscreteFitRow(
        constraint_mode=constraint_mode,
        model=report.model,
        analyzed_taxon_count=report.taxon_count,
        log_likelihood=_required_fit_value(report.log_likelihood),
        parameter_count=_required_parameter_count(report.parameter_count),
        aic=_required_fit_value(report.aic),
        root_most_likely_state=root_estimate.most_likely_state,
        root_confidence=root_estimate.confidence,
    )


def _build_node_rows(
    constrained_report,
    unconstrained_report,
) -> list[IrreversibleDiscreteNodeComparisonRow]:
    unconstrained_by_node = {
        estimate.node: estimate
        for estimate in unconstrained_report.estimates
        if not estimate.is_tip
    }
    rows: list[IrreversibleDiscreteNodeComparisonRow] = []
    for constrained_estimate in constrained_report.estimates:
        if (
            constrained_estimate.is_tip
            or constrained_estimate.node not in unconstrained_by_node
        ):
            continue
        unconstrained_estimate = unconstrained_by_node[constrained_estimate.node]
        rows.append(
            IrreversibleDiscreteNodeComparisonRow(
                node=constrained_estimate.node,
                descendant_taxa=constrained_estimate.descendant_taxa,
                constrained_state=constrained_estimate.most_likely_state,
                unconstrained_state=unconstrained_estimate.most_likely_state,
                constrained_confidence=constrained_estimate.confidence,
                unconstrained_confidence=unconstrained_estimate.confidence,
                confidence_delta=constrained_estimate.confidence
                - unconstrained_estimate.confidence,
                differs=constrained_estimate.most_likely_state
                != unconstrained_estimate.most_likely_state,
                ambiguity_changed=constrained_estimate.ambiguous
                != unconstrained_estimate.ambiguous,
            )
        )
    return rows


def _build_transition_rows(
    constrained_report,
    unconstrained_report,
) -> list[IrreversibleDiscreteTransitionComparisonRow]:
    constrained_by_pair = {
        (row.source_state, row.target_state): row
        for row in constrained_report.transition_rate_rows
    }
    unconstrained_by_pair = {
        (row.source_state, row.target_state): row
        for row in unconstrained_report.transition_rate_rows
    }
    rows: list[IrreversibleDiscreteTransitionComparisonRow] = []
    for pair in sorted(constrained_by_pair):
        constrained_row = constrained_by_pair[pair]
        unconstrained_row = unconstrained_by_pair[pair]
        rows.append(
            IrreversibleDiscreteTransitionComparisonRow(
                source_state=constrained_row.source_state,
                target_state=constrained_row.target_state,
                constrained_transition_allowed=constrained_row.transition_allowed,
                unconstrained_transition_allowed=(unconstrained_row.transition_allowed),
                constrained_rate=constrained_row.rate,
                unconstrained_rate=unconstrained_row.rate,
            )
        )
    return rows


def _required_fit_value(value: float | None) -> float:
    if value is None:
        raise ValueError(
            "irreversible discrete reconstruction requires likelihood fit details"
        )
    return value


def _required_parameter_count(value: int | None) -> int:
    if value is None:
        raise ValueError(
            "irreversible discrete reconstruction requires likelihood fit details"
        )
    return value

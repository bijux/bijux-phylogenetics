from __future__ import annotations

from dataclasses import dataclass, replace
import json
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    load_discrete_dataset,
    write_ancestral_rows,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_root_prior as _resolve_root_prior,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_state_order as _resolve_state_order,
)


@dataclass(frozen=True, slots=True)
class RootSensitivityAssumptionRow:
    """One reconstruction run under one explicit root assumption."""

    assumption_id: str
    root_prior_mode: str
    fixed_root_state: str | None
    root_prior_distribution: dict[str, float]
    root_most_likely_state: str
    root_confidence: float
    root_entropy: float
    unstable_node_count: int
    weak_support_node_count: int


@dataclass(frozen=True, slots=True)
class RootSensitivityNodeRow:
    """One node-level comparison across explicit root assumptions."""

    node: str
    descendant_taxa: list[str]
    assumption_states: dict[str, str]
    assumption_confidences: dict[str, float]
    assumption_entropies: dict[str, float]
    unique_state_count: int
    state_changed: bool
    max_confidence_delta: float
    max_entropy_delta: float
    sensitivity_score: float
    sensitivity_rank: int
    stability_class: str


@dataclass(frozen=True, slots=True)
class RootSensitivitySummary:
    """Reviewer-facing summary for ancestral root-assumption sensitivity."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    analyzed_taxon_count: int
    assumption_count: int
    compared_node_count: int
    state_changed_node_count: int
    support_changed_node_count: int
    top_sensitive_node: str | None
    top_sensitive_score: float | None
    warning_count: int


@dataclass(slots=True)
class RootSensitivityReport:
    """Discrete ancestral sensitivity to explicit root assumptions."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    analyzed_taxon_count: int
    fixed_root_state: str | None
    assumption_rows: list[RootSensitivityAssumptionRow]
    node_rows: list[RootSensitivityNodeRow]
    warnings: list[str]


def summarize_ancestral_root_sensitivity(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    fixed_root_state: str | None = None,
) -> RootSensitivityReport:
    """Compare discrete ancestral reconstructions across explicit root assumptions."""
    if model == "meristic":
        _resolve_discrete_model_name(model)
    if model not in {"equal-rates", "symmetric", "all-rates-different"}:
        raise ValueError(
            "ancestral root sensitivity requires a discrete likelihood model"
        )
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    assumption_specs: list[tuple[str, str, str | None]] = [
        ("equal_root_prior", "equal", None),
        ("empirical_root_prior", "empirical", None),
    ]
    if fixed_root_state is not None:
        assumption_specs.append(("fixed_root_state", "fixed", fixed_root_state))
    assumption_reports: list[
        tuple[RootSensitivityAssumptionRow, DiscreteAncestralReport]
    ] = []
    warnings = list(dataset.warnings)
    if fixed_root_state is not None:
        warnings.append(
            "fixed-root-state sensitivity assumes the supplied root state is biologically admissible and should be treated as a scenario test"
        )
    for assumption_id, root_prior_mode, resolved_fixed_root_state in assumption_specs:
        report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            state_ordering=state_ordering,
            ordered_states=ordered_states,
            root_prior_mode=root_prior_mode,
            fixed_root_state=resolved_fixed_root_state,
        )
        root_prior_distribution = _root_prior_distribution(
            state_order,
            state_counts=dataset.state_counts,
            root_prior_mode=root_prior_mode,
            fixed_root_state=resolved_fixed_root_state,
        )
        root_estimate = _root_estimate(report)
        assumption_reports.append(
            (
                RootSensitivityAssumptionRow(
                    assumption_id=assumption_id,
                    root_prior_mode=root_prior_mode,
                    fixed_root_state=resolved_fixed_root_state,
                    root_prior_distribution=root_prior_distribution,
                    root_most_likely_state=root_estimate.most_likely_state,
                    root_confidence=root_estimate.confidence,
                    root_entropy=_estimate_entropy(root_estimate.state_probabilities),
                    unstable_node_count=len(report.unstable_nodes),
                    weak_support_node_count=len(report.weak_support_nodes),
                ),
                report,
            )
        )
    node_rows = _build_root_sensitivity_node_rows(
        [row for row, _report in assumption_reports],
        [report for _row, report in assumption_reports],
    )
    if any(row.state_changed for row in node_rows):
        warnings.append(
            "one or more internal-node ancestral conclusions change state across explicit root assumptions"
        )
    if any(row.stability_class == "root_sensitive_support" for row in node_rows):
        warnings.append(
            "one or more internal-node ancestral conclusions remain state-stable but change support materially across explicit root assumptions"
        )
    return RootSensitivityReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=dataset.taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=list(ordered_states or []),
        analyzed_taxon_count=len(dataset.taxa),
        fixed_root_state=fixed_root_state,
        assumption_rows=[row for row, _report in assumption_reports],
        node_rows=node_rows,
        warnings=warnings,
    )


def summarize_ancestral_root_sensitivity_report(
    report: RootSensitivityReport,
) -> RootSensitivitySummary:
    """Summarize the main review facts for ancestral root-assumption sensitivity."""
    top_row = report.node_rows[0] if report.node_rows else None
    return RootSensitivitySummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        analyzed_taxon_count=report.analyzed_taxon_count,
        assumption_count=len(report.assumption_rows),
        compared_node_count=len(report.node_rows),
        state_changed_node_count=sum(row.state_changed for row in report.node_rows),
        support_changed_node_count=sum(
            row.stability_class == "root_sensitive_support" for row in report.node_rows
        ),
        top_sensitive_node=top_row.node if top_row is not None else None,
        top_sensitive_score=top_row.sensitivity_score if top_row is not None else None,
        warning_count=len(report.warnings),
    )


def write_ancestral_root_sensitivity_summary_table(
    path: Path,
    report: RootSensitivityReport,
) -> Path:
    """Write one overall summary ledger for root-assumption sensitivity."""
    summary = summarize_ancestral_root_sensitivity_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "analyzed_taxon_count",
            "assumption_count",
            "compared_node_count",
            "state_changed_node_count",
            "support_changed_node_count",
            "top_sensitive_node",
            "top_sensitive_score",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "assumption_count": str(summary.assumption_count),
                "compared_node_count": str(summary.compared_node_count),
                "state_changed_node_count": str(summary.state_changed_node_count),
                "support_changed_node_count": str(summary.support_changed_node_count),
                "top_sensitive_node": summary.top_sensitive_node or "",
                "top_sensitive_score": _format_optional_float(
                    summary.top_sensitive_score
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_ancestral_root_assumption_table(
    path: Path,
    report: RootSensitivityReport,
) -> Path:
    """Write one assumption-level ledger for root-assumption sensitivity."""
    return write_ancestral_rows(
        path,
        columns=[
            "assumption_id",
            "root_prior_mode",
            "fixed_root_state",
            "root_prior_distribution",
            "root_most_likely_state",
            "root_confidence",
            "root_entropy",
            "unstable_node_count",
            "weak_support_node_count",
        ],
        rows=[
            {
                "assumption_id": row.assumption_id,
                "root_prior_mode": row.root_prior_mode,
                "fixed_root_state": row.fixed_root_state or "",
                "root_prior_distribution": json.dumps(
                    row.root_prior_distribution,
                    sort_keys=True,
                ),
                "root_most_likely_state": row.root_most_likely_state,
                "root_confidence": str(row.root_confidence),
                "root_entropy": str(row.root_entropy),
                "unstable_node_count": str(row.unstable_node_count),
                "weak_support_node_count": str(row.weak_support_node_count),
            }
            for row in report.assumption_rows
        ],
    )


def write_ancestral_root_sensitivity_node_table(
    path: Path,
    report: RootSensitivityReport,
) -> Path:
    """Write one node-wise comparison ledger for root-assumption sensitivity."""
    return write_ancestral_rows(
        path,
        columns=[
            "node",
            "descendant_taxa",
            "assumption_states",
            "assumption_confidences",
            "assumption_entropies",
            "unique_state_count",
            "state_changed",
            "max_confidence_delta",
            "max_entropy_delta",
            "sensitivity_score",
            "sensitivity_rank",
            "stability_class",
        ],
        rows=[
            {
                "node": row.node,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "assumption_states": json.dumps(row.assumption_states, sort_keys=True),
                "assumption_confidences": json.dumps(
                    row.assumption_confidences,
                    sort_keys=True,
                ),
                "assumption_entropies": json.dumps(
                    row.assumption_entropies,
                    sort_keys=True,
                ),
                "unique_state_count": str(row.unique_state_count),
                "state_changed": str(row.state_changed).lower(),
                "max_confidence_delta": str(row.max_confidence_delta),
                "max_entropy_delta": str(row.max_entropy_delta),
                "sensitivity_score": str(row.sensitivity_score),
                "sensitivity_rank": str(row.sensitivity_rank),
                "stability_class": row.stability_class,
            }
            for row in report.node_rows
        ],
    )


def _build_root_sensitivity_node_rows(
    assumption_rows: list[RootSensitivityAssumptionRow],
    reports: list[DiscreteAncestralReport],
) -> list[RootSensitivityNodeRow]:
    node_to_estimates: dict[str, dict[str, object]] = {}
    descendant_taxa: dict[str, list[str]] = {}
    for assumption_row, report in zip(assumption_rows, reports, strict=True):
        for estimate in report.estimates:
            if estimate.is_tip:
                continue
            node_to_estimates.setdefault(estimate.node, {})[
                assumption_row.assumption_id
            ] = estimate
            descendant_taxa[estimate.node] = estimate.descendant_taxa
    rows: list[RootSensitivityNodeRow] = []
    for node, estimates_by_assumption in sorted(node_to_estimates.items()):
        assumption_states = {
            assumption_id: estimate.most_likely_state
            for assumption_id, estimate in sorted(estimates_by_assumption.items())
        }
        assumption_confidences = {
            assumption_id: estimate.confidence
            for assumption_id, estimate in sorted(estimates_by_assumption.items())
        }
        assumption_entropies = {
            assumption_id: _estimate_entropy(estimate.state_probabilities)
            for assumption_id, estimate in sorted(estimates_by_assumption.items())
        }
        confidence_values = list(assumption_confidences.values())
        entropy_values = list(assumption_entropies.values())
        unique_state_count = len(set(assumption_states.values()))
        state_changed = unique_state_count > 1
        max_confidence_delta = max(confidence_values) - min(confidence_values)
        max_entropy_delta = max(entropy_values) - min(entropy_values)
        sensitivity_score = (
            float(state_changed) + max_confidence_delta + max_entropy_delta
        )
        if state_changed:
            stability_class = "root_sensitive_state"
        elif max_confidence_delta >= 0.15 or max_entropy_delta >= 0.15:
            stability_class = "root_sensitive_support"
        else:
            stability_class = "stable"
        rows.append(
            RootSensitivityNodeRow(
                node=node,
                descendant_taxa=descendant_taxa[node],
                assumption_states=assumption_states,
                assumption_confidences=assumption_confidences,
                assumption_entropies=assumption_entropies,
                unique_state_count=unique_state_count,
                state_changed=state_changed,
                max_confidence_delta=max_confidence_delta,
                max_entropy_delta=max_entropy_delta,
                sensitivity_score=sensitivity_score,
                sensitivity_rank=0,
                stability_class=stability_class,
            )
        )
    ranked_rows = sorted(
        rows,
        key=lambda row: (-row.sensitivity_score, row.node),
    )
    return [
        replace(row, sensitivity_rank=index)
        for index, row in enumerate(ranked_rows, start=1)
    ]


def _root_estimate(report: DiscreteAncestralReport):
    internal_estimates = [
        estimate for estimate in report.estimates if not estimate.is_tip
    ]
    return max(
        internal_estimates,
        key=lambda estimate: (len(estimate.descendant_taxa), estimate.node),
    )


def _root_prior_distribution(
    state_order: list[str],
    *,
    state_counts: dict[str, int],
    root_prior_mode: str,
    fixed_root_state: str | None,
) -> dict[str, float]:
    root_prior = _resolve_root_prior(
        state_order,
        state_counts=state_counts,
        mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
    )
    return {state: float(root_prior[index]) for index, state in enumerate(state_order)}


def _estimate_entropy(state_probabilities: dict[str, float]) -> float:
    positive_probabilities = [
        probability for probability in state_probabilities.values() if probability > 0.0
    ]
    if len(positive_probabilities) <= 1:
        return 0.0
    return -sum(
        probability * math.log2(probability) for probability in positive_probabilities
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)

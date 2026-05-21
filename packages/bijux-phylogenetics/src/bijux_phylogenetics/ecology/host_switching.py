from __future__ import annotations

import csv
from dataclasses import dataclass
import json
from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError


@dataclass(frozen=True, slots=True)
class HostStateNodeRow:
    """One internal-node host-state reconstruction row."""

    node: str
    node_name: str | None
    descendant_taxa: list[str]
    most_likely_host: str
    host_probabilities: dict[str, float]
    confidence: float
    ambiguous: bool
    is_root: bool


@dataclass(frozen=True, slots=True)
class HostSwitchBranchRow:
    """One branchwise host-switch review row."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_host: str
    child_most_likely_host: str
    parent_host_set: list[str]
    child_host_set: list[str]
    overlapping_hosts: list[str]
    changed: bool
    transition: str
    certainty_class: str
    parent_confidence: float
    child_confidence: float
    transition_allowed: bool


@dataclass(frozen=True, slots=True)
class HostSwitchCountRow:
    """One aggregated host-switch count row."""

    transition: str
    source_host: str
    target_host: str
    transition_allowed: bool
    certain_switch_count: int
    uncertain_switch_count: int
    total_switch_count: int


@dataclass(frozen=True, slots=True)
class HostSwitchFitRow:
    """One fitted host-transition regime row."""

    constraint_mode: str
    model: str
    analyzed_taxon_count: int
    log_likelihood: float
    parameter_count: int
    aic: float
    root_host: str
    root_confidence: float


@dataclass(frozen=True, slots=True)
class UnsupportedHostSwitchClaimRow:
    """One unconstrained host-switch claim forbidden by the supplied transition policy."""

    branch_id: str
    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    unconstrained_source_host: str
    unconstrained_target_host: str
    unconstrained_certainty_class: str
    constrained_source_host: str
    constrained_target_host: str
    constrained_certainty_class: str
    claim_resolved: bool


@dataclass(frozen=True, slots=True)
class HostSwitchExclusionRow:
    """One excluded row from the host-switching workflow."""

    taxon: str
    raw_host: str
    normalized_host: str | None
    reason: str
    note: str


@dataclass(frozen=True, slots=True)
class HostSwitchSummary:
    """Reviewer-facing summary for host-switching analysis."""

    trait: str
    taxon_column: str
    model: str
    internal_model: str
    analysis_constraint_mode: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    observed_host_count: int
    internal_node_count: int
    ambiguous_internal_node_count: int
    host_switch_count: int
    certain_host_switch_count: int
    uncertain_host_switch_count: int
    allowed_transition_count: int
    forbidden_transition_count: int
    constrained_log_likelihood: float | None
    unconstrained_log_likelihood: float
    constrained_aic: float | None
    unconstrained_aic: float
    preferred_constraint: str
    unsupported_switch_claim_count: int
    root_host: str
    root_confidence: float
    warning_count: int


@dataclass(slots=True)
class HostSwitchingReport:
    """Owned host-switching review surface on one parasite or pathogen tree."""

    tree_path: Path
    traits_path: Path
    constraint_path: Path | None
    trait: str
    taxon_column: str
    model: str
    internal_model: str
    analysis_constraint_mode: str
    summary: HostSwitchSummary
    node_rows: list[HostStateNodeRow]
    branch_rows: list[HostSwitchBranchRow]
    count_rows: list[HostSwitchCountRow]
    fit_rows: list[HostSwitchFitRow]
    unsupported_claim_rows: list[UnsupportedHostSwitchClaimRow]
    exclusion_rows: list[HostSwitchExclusionRow]
    warnings: list[str]


def summarize_host_switching(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "ard",
    constraint_path: Path | None = None,
) -> HostSwitchingReport:
    """Model host-state evolution and summarize inferred host switches."""
    unconstrained_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    exclusion_rows = _build_exclusion_rows(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    allowed_transition_pairs: list[tuple[str, str]] | None = None
    constrained_report: DiscreteAncestralReport | None = None
    if constraint_path is not None:
        allowed_transition_pairs = _load_allowed_host_transitions(
            constraint_path,
            observed_hosts=unconstrained_report.observed_states,
        )
        constrained_report = reconstruct_discrete_ancestral_states(
            tree_path,
            traits_path,
            trait=trait,
            taxon_column=taxon_column,
            model=model,
            allowed_transition_pairs=allowed_transition_pairs,
        )
    active_report = constrained_report or unconstrained_report
    active_allowed_transition_pairs = (
        set(allowed_transition_pairs)
        if allowed_transition_pairs is not None
        else {
            (left, right)
            for left in active_report.observed_states
            for right in active_report.observed_states
            if left != right
        }
    )
    node_rows = _build_node_rows(active_report)
    branch_rows = _build_branch_rows(
        active_report,
        allowed_transition_pairs=active_allowed_transition_pairs,
    )
    count_rows = _build_count_rows(branch_rows)
    fit_rows = _build_fit_rows(
        unconstrained_report=unconstrained_report,
        constrained_report=constrained_report,
    )
    unsupported_claim_rows = _build_unsupported_claim_rows(
        unconstrained_report=unconstrained_report,
        constrained_report=constrained_report,
        allowed_transition_pairs=active_allowed_transition_pairs,
    )
    warnings = list(
        dict.fromkeys(
            [
                *unconstrained_report.warnings,
                *(
                    constrained_report.warnings
                    if constrained_report is not None
                    else []
                ),
            ]
        )
    )
    summary = _build_summary(
        active_report=active_report,
        unconstrained_report=unconstrained_report,
        constrained_report=constrained_report,
        branch_rows=branch_rows,
        count_rows=count_rows,
        unsupported_claim_rows=unsupported_claim_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
        analysis_constraint_mode=(
            "constrained" if constrained_report is not None else "unconstrained"
        ),
    )
    return HostSwitchingReport(
        tree_path=tree_path,
        traits_path=traits_path,
        constraint_path=constraint_path,
        trait=active_report.trait,
        taxon_column=active_report.taxon_column,
        model=model,
        internal_model=active_report.model,
        analysis_constraint_mode=summary.analysis_constraint_mode,
        summary=summary,
        node_rows=node_rows,
        branch_rows=branch_rows,
        count_rows=count_rows,
        fit_rows=fit_rows,
        unsupported_claim_rows=unsupported_claim_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )


def write_host_switch_summary_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one overall summary ledger for host-switching analysis."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "internal_model",
            "analysis_constraint_mode",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "observed_host_count",
            "internal_node_count",
            "ambiguous_internal_node_count",
            "host_switch_count",
            "certain_host_switch_count",
            "uncertain_host_switch_count",
            "allowed_transition_count",
            "forbidden_transition_count",
            "constrained_log_likelihood",
            "unconstrained_log_likelihood",
            "constrained_aic",
            "unconstrained_aic",
            "preferred_constraint",
            "unsupported_switch_claim_count",
            "root_host",
            "root_confidence",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "internal_model": summary.internal_model,
                "analysis_constraint_mode": summary.analysis_constraint_mode,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "observed_host_count": str(summary.observed_host_count),
                "internal_node_count": str(summary.internal_node_count),
                "ambiguous_internal_node_count": str(
                    summary.ambiguous_internal_node_count
                ),
                "host_switch_count": str(summary.host_switch_count),
                "certain_host_switch_count": str(summary.certain_host_switch_count),
                "uncertain_host_switch_count": str(summary.uncertain_host_switch_count),
                "allowed_transition_count": str(summary.allowed_transition_count),
                "forbidden_transition_count": str(summary.forbidden_transition_count),
                "constrained_log_likelihood": _format_optional_float(
                    summary.constrained_log_likelihood
                ),
                "unconstrained_log_likelihood": str(
                    summary.unconstrained_log_likelihood
                ),
                "constrained_aic": _format_optional_float(summary.constrained_aic),
                "unconstrained_aic": str(summary.unconstrained_aic),
                "preferred_constraint": summary.preferred_constraint,
                "unsupported_switch_claim_count": str(
                    summary.unsupported_switch_claim_count
                ),
                "root_host": summary.root_host,
                "root_confidence": str(summary.root_confidence),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_host_state_node_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one internal-node host-state ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "descendant_taxa",
            "most_likely_host",
            "host_probabilities",
            "confidence",
            "ambiguous",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "descendant_taxa": ",".join(row.descendant_taxa),
                "most_likely_host": row.most_likely_host,
                "host_probabilities": json.dumps(
                    row.host_probabilities,
                    sort_keys=True,
                ),
                "confidence": str(row.confidence),
                "ambiguous": str(row.ambiguous).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.node_rows
        ],
    )


def write_host_switch_branch_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one branchwise host-switch ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_host",
            "child_most_likely_host",
            "parent_host_set",
            "child_host_set",
            "overlapping_hosts",
            "changed",
            "transition",
            "certainty_class",
            "parent_confidence",
            "child_confidence",
            "transition_allowed",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": _format_optional_float(row.branch_length),
                "parent_most_likely_host": row.parent_most_likely_host,
                "child_most_likely_host": row.child_most_likely_host,
                "parent_host_set": ",".join(row.parent_host_set),
                "child_host_set": ",".join(row.child_host_set),
                "overlapping_hosts": ",".join(row.overlapping_hosts),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
                "parent_confidence": str(row.parent_confidence),
                "child_confidence": str(row.child_confidence),
                "transition_allowed": str(row.transition_allowed).lower(),
            }
            for row in report.branch_rows
        ],
    )


def write_host_switch_count_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one aggregated host-switch count ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "transition",
            "source_host",
            "target_host",
            "transition_allowed",
            "certain_switch_count",
            "uncertain_switch_count",
            "total_switch_count",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_host": row.source_host,
                "target_host": row.target_host,
                "transition_allowed": str(row.transition_allowed).lower(),
                "certain_switch_count": str(row.certain_switch_count),
                "uncertain_switch_count": str(row.uncertain_switch_count),
                "total_switch_count": str(row.total_switch_count),
            }
            for row in report.count_rows
        ],
    )


def write_host_switch_fit_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one fit-comparison ledger for unconstrained and constrained host models."""
    return write_taxon_rows(
        path,
        columns=[
            "constraint_mode",
            "model",
            "analyzed_taxon_count",
            "log_likelihood",
            "parameter_count",
            "aic",
            "root_host",
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
                "root_host": row.root_host,
                "root_confidence": str(row.root_confidence),
            }
            for row in report.fit_rows
        ],
    )


def write_unsupported_host_switch_claim_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one forbidden unconstrained host-switch claim ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "unconstrained_source_host",
            "unconstrained_target_host",
            "unconstrained_certainty_class",
            "constrained_source_host",
            "constrained_target_host",
            "constrained_certainty_class",
            "claim_resolved",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "unconstrained_source_host": row.unconstrained_source_host,
                "unconstrained_target_host": row.unconstrained_target_host,
                "unconstrained_certainty_class": row.unconstrained_certainty_class,
                "constrained_source_host": row.constrained_source_host,
                "constrained_target_host": row.constrained_target_host,
                "constrained_certainty_class": row.constrained_certainty_class,
                "claim_resolved": str(row.claim_resolved).lower(),
            }
            for row in report.unsupported_claim_rows
        ],
    )


def write_host_switch_exclusion_table(
    path: Path,
    report: HostSwitchingReport,
) -> Path:
    """Write one excluded-row ledger for host-switching analysis."""
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_host",
            "normalized_host",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_host": row.raw_host,
                "normalized_host": row.normalized_host or "",
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )


def _build_node_rows(report: DiscreteAncestralReport) -> list[HostStateNodeRow]:
    root_node = report.estimates[0].node
    return [
        HostStateNodeRow(
            node=estimate.node,
            node_name=estimate.node_name,
            descendant_taxa=list(estimate.descendant_taxa),
            most_likely_host=estimate.most_likely_state,
            host_probabilities=dict(sorted(estimate.state_probabilities.items())),
            confidence=_stable_float(estimate.confidence),
            ambiguous=estimate.ambiguous,
            is_root=estimate.node == root_node,
        )
        for estimate in report.estimates
        if not estimate.is_tip
    ]


def _build_branch_rows(
    report: DiscreteAncestralReport,
    *,
    allowed_transition_pairs: set[tuple[str, str]],
) -> list[HostSwitchBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[HostSwitchBranchRow] = []

    def visit(node) -> None:
        parent_estimate = estimate_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[_node_signature(child)]
            overlapping_hosts = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = (
                parent_estimate.most_likely_state != child_estimate.most_likely_state
            )
            transition = (
                f"{parent_estimate.most_likely_state}->{child_estimate.most_likely_state}"
                if changed
                else ""
            )
            branch_rows.append(
                HostSwitchBranchRow(
                    branch_id=child_estimate.node,
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=list(child_estimate.descendant_taxa),
                    branch_length=child.branch_length,
                    parent_most_likely_host=parent_estimate.most_likely_state,
                    child_most_likely_host=child_estimate.most_likely_state,
                    parent_host_set=list(parent_estimate.state_set),
                    child_host_set=list(child_estimate.state_set),
                    overlapping_hosts=overlapping_hosts,
                    changed=changed,
                    transition=transition,
                    certainty_class=_transition_certainty_class(
                        changed=changed,
                        overlapping_hosts=overlapping_hosts,
                        parent_host_set=parent_estimate.state_set,
                        child_host_set=child_estimate.state_set,
                    ),
                    parent_confidence=_stable_float(parent_estimate.confidence),
                    child_confidence=_stable_float(child_estimate.confidence),
                    transition_allowed=(
                        not changed
                        or (
                            parent_estimate.most_likely_state,
                            child_estimate.most_likely_state,
                        )
                        in allowed_transition_pairs
                    ),
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def _build_count_rows(
    branch_rows: list[HostSwitchBranchRow],
) -> list[HostSwitchCountRow]:
    grouped: dict[str, list[HostSwitchBranchRow]] = {}
    for row in branch_rows:
        if row.changed:
            grouped.setdefault(row.transition, []).append(row)
    count_rows: list[HostSwitchCountRow] = []
    for transition in sorted(grouped):
        rows = grouped[transition]
        source_host, target_host = transition.split("->", maxsplit=1)
        certain_switch_count = sum(
            row.certainty_class == "certain_switch" for row in rows
        )
        uncertain_switch_count = sum(
            row.certainty_class == "uncertain_switch" for row in rows
        )
        count_rows.append(
            HostSwitchCountRow(
                transition=transition,
                source_host=source_host,
                target_host=target_host,
                transition_allowed=rows[0].transition_allowed,
                certain_switch_count=certain_switch_count,
                uncertain_switch_count=uncertain_switch_count,
                total_switch_count=len(rows),
            )
        )
    return count_rows


def _build_fit_rows(
    *,
    unconstrained_report: DiscreteAncestralReport,
    constrained_report: DiscreteAncestralReport | None,
) -> list[HostSwitchFitRow]:
    rows = [_build_fit_row("unconstrained", unconstrained_report)]
    if constrained_report is not None:
        rows.append(_build_fit_row("constrained", constrained_report))
    return rows


def _build_fit_row(
    constraint_mode: str,
    report: DiscreteAncestralReport,
) -> HostSwitchFitRow:
    root_estimate = next(
        estimate
        for estimate in report.estimates
        if estimate.node == report.estimates[0].node
    )
    if (
        report.log_likelihood is None
        or report.parameter_count is None
        or report.aic is None
    ):
        raise AncestralReconstructionError(
            "host-switching fit comparison requires a likelihood discrete ancestral model"
        )
    return HostSwitchFitRow(
        constraint_mode=constraint_mode,
        model=report.model,
        analyzed_taxon_count=report.taxon_count,
        log_likelihood=_stable_float(report.log_likelihood),
        parameter_count=report.parameter_count,
        aic=_stable_float(report.aic),
        root_host=root_estimate.most_likely_state,
        root_confidence=_stable_float(root_estimate.confidence),
    )


def _build_unsupported_claim_rows(
    *,
    unconstrained_report: DiscreteAncestralReport,
    constrained_report: DiscreteAncestralReport | None,
    allowed_transition_pairs: set[tuple[str, str]],
) -> list[UnsupportedHostSwitchClaimRow]:
    if constrained_report is None:
        return []
    unconstrained_rows = {
        row.branch_id: row
        for row in _build_branch_rows(
            unconstrained_report,
            allowed_transition_pairs=allowed_transition_pairs,
        )
        if row.changed and not row.transition_allowed
    }
    constrained_rows = {
        row.branch_id: row
        for row in _build_branch_rows(
            constrained_report,
            allowed_transition_pairs=allowed_transition_pairs,
        )
    }
    rows: list[UnsupportedHostSwitchClaimRow] = []
    for branch_id, unconstrained_row in sorted(unconstrained_rows.items()):
        constrained_row = constrained_rows[branch_id]
        rows.append(
            UnsupportedHostSwitchClaimRow(
                branch_id=branch_id,
                parent_node=unconstrained_row.parent_node,
                child_node=unconstrained_row.child_node,
                child_descendant_taxa=list(unconstrained_row.child_descendant_taxa),
                unconstrained_source_host=unconstrained_row.parent_most_likely_host,
                unconstrained_target_host=unconstrained_row.child_most_likely_host,
                unconstrained_certainty_class=unconstrained_row.certainty_class,
                constrained_source_host=constrained_row.parent_most_likely_host,
                constrained_target_host=constrained_row.child_most_likely_host,
                constrained_certainty_class=constrained_row.certainty_class,
                claim_resolved=(
                    unconstrained_row.transition != constrained_row.transition
                ),
            )
        )
    return rows


def _build_exclusion_rows(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
) -> list[HostSwitchExclusionRow]:
    audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return [
        HostSwitchExclusionRow(
            taxon=row.taxon,
            raw_host=row.raw_state,
            normalized_host=row.normalized_state,
            reason=row.issue_code or "excluded",
            note=row.note,
        )
        for row in audit.rows
        if not row.included
    ]


def _build_summary(
    *,
    active_report: DiscreteAncestralReport,
    unconstrained_report: DiscreteAncestralReport,
    constrained_report: DiscreteAncestralReport | None,
    branch_rows: list[HostSwitchBranchRow],
    count_rows: list[HostSwitchCountRow],
    unsupported_claim_rows: list[UnsupportedHostSwitchClaimRow],
    exclusion_rows: list[HostSwitchExclusionRow],
    warnings: list[str],
    analysis_constraint_mode: str,
) -> HostSwitchSummary:
    root_estimate = next(
        estimate
        for estimate in active_report.estimates
        if estimate.node == active_report.estimates[0].node
    )
    certain_host_switch_count = sum(
        row.certainty_class == "certain_switch" for row in branch_rows if row.changed
    )
    uncertain_host_switch_count = sum(
        row.certainty_class == "uncertain_switch" for row in branch_rows if row.changed
    )
    observed_host_count = len(active_report.observed_states)
    allowed_transition_count = len(active_report.allowed_transition_pairs)
    all_transition_count = observed_host_count * max(observed_host_count - 1, 0)
    forbidden_transition_count = max(all_transition_count - allowed_transition_count, 0)
    preferred_constraint = "unconstrained"
    if (
        constrained_report is not None
        and constrained_report.aic is not None
        and constrained_report.aic
        <= (unconstrained_report.aic or constrained_report.aic)
    ):
        preferred_constraint = "constrained"
    return HostSwitchSummary(
        trait=active_report.trait,
        taxon_column=active_report.taxon_column,
        model=unconstrained_report.model,
        internal_model=active_report.model,
        analysis_constraint_mode=analysis_constraint_mode,
        analyzed_taxon_count=active_report.taxon_count,
        excluded_taxon_count=len(exclusion_rows),
        observed_host_count=observed_host_count,
        internal_node_count=sum(
            not estimate.is_tip for estimate in active_report.estimates
        ),
        ambiguous_internal_node_count=sum(
            row.ambiguous for row in _build_node_rows(active_report)
        ),
        host_switch_count=sum(row.total_switch_count for row in count_rows),
        certain_host_switch_count=certain_host_switch_count,
        uncertain_host_switch_count=uncertain_host_switch_count,
        allowed_transition_count=allowed_transition_count,
        forbidden_transition_count=forbidden_transition_count,
        constrained_log_likelihood=(
            _stable_float(constrained_report.log_likelihood)
            if constrained_report is not None
            and constrained_report.log_likelihood is not None
            else None
        ),
        unconstrained_log_likelihood=_stable_float(
            unconstrained_report.log_likelihood or 0.0
        ),
        constrained_aic=(
            _stable_float(constrained_report.aic)
            if constrained_report is not None and constrained_report.aic is not None
            else None
        ),
        unconstrained_aic=_stable_float(unconstrained_report.aic or 0.0),
        preferred_constraint=preferred_constraint,
        unsupported_switch_claim_count=len(unsupported_claim_rows),
        root_host=root_estimate.most_likely_state,
        root_confidence=_stable_float(root_estimate.confidence),
        warning_count=len(warnings),
    )


def _load_allowed_host_transitions(
    path: Path,
    *,
    observed_hosts: list[str],
) -> list[tuple[str, str]]:
    if not path.exists():
        raise FileNotFoundError(f"host-transition constraint file not found: {path}")
    raw_text = path.read_text(encoding="utf-8")
    sample = raw_text[:1024]
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        dialect = csv.excel_tab if "\t" in sample else csv.excel
    reader = csv.DictReader(raw_text.splitlines(), dialect=dialect)
    if reader.fieldnames is None:
        raise ValueError("host-transition constraint file must contain a header row")
    required = {"source_host", "target_host"}
    if not required.issubset(reader.fieldnames):
        raise ValueError(
            "host-transition constraint file must contain source_host and target_host columns"
        )
    allowed_field = (
        "transition_allowed" if "transition_allowed" in reader.fieldnames else None
    )
    observed_host_set = set(observed_hosts)
    allowed_pairs: list[tuple[str, str]] = []
    for row in reader:
        source_host = (row.get("source_host") or "").strip()
        target_host = (row.get("target_host") or "").strip()
        if not source_host or not target_host:
            raise ValueError(
                "host-transition constraint rows must name both source_host and target_host"
            )
        if source_host == target_host:
            raise ValueError(
                "host-transition constraint rows must connect distinct hosts"
            )
        if source_host not in observed_host_set:
            raise ValueError(
                "host-transition source host is not present in the analyzed host vocabulary: "
                f"{source_host}"
            )
        if target_host not in observed_host_set:
            raise ValueError(
                "host-transition target host is not present in the analyzed host vocabulary: "
                f"{target_host}"
            )
        if allowed_field is not None and not _parse_truthy_cell(
            row.get(allowed_field, "")
        ):
            continue
        allowed_pairs.append((source_host, target_host))
    if not allowed_pairs:
        raise ValueError(
            "host-transition constraint file must allow at least one directed host transition"
        )
    return sorted(set(allowed_pairs))


def _parse_truthy_cell(raw: str) -> bool:
    normalized = raw.strip().lower()
    if normalized in {"", "0", "false", "no", "forbidden"}:
        return False
    if normalized in {"1", "true", "yes", "allowed", "x"}:
        return True
    raise ValueError(
        "host-transition constraint transition_allowed cells must be one of "
        "0,1,false,true,no,yes,forbidden,allowed,x"
    )


def _transition_certainty_class(
    *,
    changed: bool,
    overlapping_hosts: list[str],
    parent_host_set: list[str],
    child_host_set: list[str],
) -> str:
    if not changed:
        return "no_switch"
    if overlapping_hosts:
        return "uncertain_switch"
    if len(parent_host_set) == 1 and len(child_host_set) == 1:
        return "certain_switch"
    return "uncertain_switch"


def _node_signature(node) -> str:
    taxa = sorted(_node_descendant_taxa(node))
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


def _node_descendant_taxa(node) -> list[str]:
    if not node.children:
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_node_descendant_taxa(child))
    return taxa


def _stable_float(value: float) -> float:
    return float(format(round(value, 15), ".15g"))


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(value)

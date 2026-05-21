from __future__ import annotations

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
from .contracts import (
    HostStateNodeRow,
    HostSwitchBranchRow,
    HostSwitchCountRow,
    HostSwitchExclusionRow,
    HostSwitchFitRow,
    HostSwitchSummary,
    HostSwitchingReport,
    UnsupportedHostSwitchClaimRow,
)
from .shared import (
    format_optional_float,
    load_allowed_host_transitions,
    node_signature,
    resolve_host_state,
    stable_float,
    transition_certainty_class,
)


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
    coding_audit = audit_discrete_state_coding(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    unconstrained_report = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    observed_hosts_by_taxon = {
        row.taxon: row.normalized_state
        for row in coding_audit.rows
        if row.included and row.normalized_state is not None
    }
    exclusion_rows = _build_exclusion_rows(coding_audit)
    allowed_transition_pairs: list[tuple[str, str]] | None = None
    constrained_report: DiscreteAncestralReport | None = None
    if constraint_path is not None:
        allowed_transition_pairs = load_allowed_host_transitions(
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
        observed_hosts_by_taxon=observed_hosts_by_taxon,
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
        observed_hosts_by_taxon=observed_hosts_by_taxon,
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
                "constrained_log_likelihood": format_optional_float(
                    summary.constrained_log_likelihood
                ),
                "unconstrained_log_likelihood": str(
                    summary.unconstrained_log_likelihood
                ),
                "constrained_aic": format_optional_float(summary.constrained_aic),
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
                "branch_length": format_optional_float(row.branch_length),
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
            confidence=stable_float(estimate.confidence),
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
    observed_hosts_by_taxon: dict[str, str],
) -> list[HostSwitchBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[HostSwitchBranchRow] = []

    def visit(node) -> None:
        parent_estimate = estimate_by_node[node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[node_signature(child)]
            parent_host = resolve_host_state(
                descendant_taxa=parent_estimate.descendant_taxa,
                candidate_hosts=parent_estimate.state_set,
                observed_hosts_by_taxon=observed_hosts_by_taxon,
                fallback_host=parent_estimate.most_likely_state,
            )
            child_host = resolve_host_state(
                descendant_taxa=child_estimate.descendant_taxa,
                candidate_hosts=child_estimate.state_set,
                observed_hosts_by_taxon=observed_hosts_by_taxon,
                fallback_host=child_estimate.most_likely_state,
            )
            overlapping_hosts = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = parent_host != child_host
            transition = f"{parent_host}->{child_host}" if changed else ""
            branch_rows.append(
                HostSwitchBranchRow(
                    branch_id=child_estimate.node,
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=list(child_estimate.descendant_taxa),
                    branch_length=child.branch_length,
                    parent_most_likely_host=parent_host,
                    child_most_likely_host=child_host,
                    parent_host_set=list(parent_estimate.state_set),
                    child_host_set=list(child_estimate.state_set),
                    overlapping_hosts=overlapping_hosts,
                    changed=changed,
                    transition=transition,
                    certainty_class=transition_certainty_class(
                        changed=changed,
                        overlapping_hosts=overlapping_hosts,
                        parent_host_set=parent_estimate.state_set,
                        child_host_set=child_estimate.state_set,
                    ),
                    parent_confidence=stable_float(parent_estimate.confidence),
                    child_confidence=stable_float(child_estimate.confidence),
                    transition_allowed=(
                        not changed
                        or (parent_host, child_host) in allowed_transition_pairs
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
        log_likelihood=stable_float(report.log_likelihood),
        parameter_count=report.parameter_count,
        aic=stable_float(report.aic),
        root_host=root_estimate.most_likely_state,
        root_confidence=stable_float(root_estimate.confidence),
    )


def _build_unsupported_claim_rows(
    *,
    unconstrained_report: DiscreteAncestralReport,
    constrained_report: DiscreteAncestralReport | None,
    allowed_transition_pairs: set[tuple[str, str]],
    observed_hosts_by_taxon: dict[str, str],
) -> list[UnsupportedHostSwitchClaimRow]:
    if constrained_report is None:
        return []
    unconstrained_rows = {
        row.branch_id: row
        for row in _build_branch_rows(
            unconstrained_report,
            allowed_transition_pairs=allowed_transition_pairs,
            observed_hosts_by_taxon=observed_hosts_by_taxon,
        )
        if row.changed and not row.transition_allowed
    }
    constrained_rows = {
        row.branch_id: row
        for row in _build_branch_rows(
            constrained_report,
            allowed_transition_pairs=allowed_transition_pairs,
            observed_hosts_by_taxon=observed_hosts_by_taxon,
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


def _build_exclusion_rows(audit) -> list[HostSwitchExclusionRow]:
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
            stable_float(constrained_report.log_likelihood)
            if constrained_report is not None
            and constrained_report.log_likelihood is not None
            else None
        ),
        unconstrained_log_likelihood=stable_float(
            unconstrained_report.log_likelihood or 0.0
        ),
        constrained_aic=(
            stable_float(constrained_report.aic)
            if constrained_report is not None and constrained_report.aic is not None
            else None
        ),
        unconstrained_aic=stable_float(unconstrained_report.aic or 0.0),
        preferred_constraint=preferred_constraint,
        unsupported_switch_claim_count=len(unsupported_claim_rows),
        root_host=root_estimate.most_likely_state,
        root_confidence=stable_float(root_estimate.confidence),
        warning_count=len(warnings),
    )

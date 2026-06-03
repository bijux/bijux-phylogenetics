from __future__ import annotations

from bijux_phylogenetics.ancestral.discrete import DiscreteAncestralReport
from bijux_phylogenetics.io.newick import loads_newick
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .contracts import (
    HostStateNodeRow,
    HostSwitchBranchRow,
    HostSwitchCountRow,
    HostSwitchExclusionRow,
    HostSwitchFitRow,
    HostSwitchSummary,
    UnsupportedHostSwitchClaimRow,
)
from .shared import (
    node_signature,
    resolve_host_state,
    stable_float,
    transition_certainty_class,
)


def build_node_rows(report: DiscreteAncestralReport) -> list[HostStateNodeRow]:
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


def build_branch_rows(
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


def build_count_rows(
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


def build_fit_rows(
    *,
    unconstrained_report: DiscreteAncestralReport,
    constrained_report: DiscreteAncestralReport | None,
) -> list[HostSwitchFitRow]:
    rows = [build_fit_row("unconstrained", unconstrained_report)]
    if constrained_report is not None:
        rows.append(build_fit_row("constrained", constrained_report))
    return rows


def build_fit_row(
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


def build_unsupported_claim_rows(
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
        for row in build_branch_rows(
            unconstrained_report,
            allowed_transition_pairs=allowed_transition_pairs,
            observed_hosts_by_taxon=observed_hosts_by_taxon,
        )
        if row.changed and not row.transition_allowed
    }
    constrained_rows = {
        row.branch_id: row
        for row in build_branch_rows(
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


def build_exclusion_rows(audit) -> list[HostSwitchExclusionRow]:
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


def build_summary(
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
            row.ambiguous for row in build_node_rows(active_report)
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

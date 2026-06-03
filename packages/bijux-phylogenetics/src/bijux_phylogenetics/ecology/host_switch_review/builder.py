from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    audit_discrete_state_coding,
)

from .analysis import (
    build_branch_rows,
    build_count_rows,
    build_exclusion_rows,
    build_fit_rows,
    build_node_rows,
    build_summary,
    build_unsupported_claim_rows,
)
from .contracts import HostSwitchingReport
from .shared import load_allowed_host_transitions


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
    exclusion_rows = build_exclusion_rows(coding_audit)
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
    node_rows = build_node_rows(active_report)
    branch_rows = build_branch_rows(
        active_report,
        allowed_transition_pairs=active_allowed_transition_pairs,
        observed_hosts_by_taxon=observed_hosts_by_taxon,
    )
    count_rows = build_count_rows(branch_rows)
    fit_rows = build_fit_rows(
        unconstrained_report=unconstrained_report,
        constrained_report=constrained_report,
    )
    unsupported_claim_rows = build_unsupported_claim_rows(
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
    summary = build_summary(
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

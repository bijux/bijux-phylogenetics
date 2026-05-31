from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.discrete import reconstruct_discrete_ancestral_states

from .aggregation import _summarize_transition_rows
from .branch_analysis import _build_transition_branch_rows
from .contracts import (
    AncestralTransitionExclusion,
    AncestralTransitionReport,
    AncestralTransitionSummary,
)


def summarize_ancestral_transitions(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "fitch",
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
) -> AncestralTransitionReport:
    """Count inferred discrete ancestral transitions for one rooted tree."""
    reconstruction = reconstruct_discrete_ancestral_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    branch_rows = _build_transition_branch_rows(reconstruction)
    transition_rows = _summarize_transition_rows(branch_rows)
    exclusions = [
        AncestralTransitionExclusion(
            taxon=taxon,
            reason="missing_discrete_trait_state",
        )
        for taxon in reconstruction.dropped_missing_taxa
    ]
    warnings = list(reconstruction.warnings)
    if any(row.uncertain_change_count > 0 for row in transition_rows):
        warnings.append(
            "one or more inferred ancestral transitions remain uncertain because parent and child state sets overlap"
        )
    return AncestralTransitionReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_column=reconstruction.taxon_column,
        model=reconstruction.model,
        state_ordering=reconstruction.state_ordering,
        ordered_states=reconstruction.ordered_states,
        taxon_count=reconstruction.taxon_count,
        branch_rows=branch_rows,
        transition_rows=transition_rows,
        exclusions=exclusions,
        warnings=warnings,
    )


def summarize_ancestral_transition_report(
    report: AncestralTransitionReport,
) -> AncestralTransitionSummary:
    """Summarize the main review facts for one ancestral transition report."""
    changed_rows = [row for row in report.branch_rows if row.changed]
    certain_change_count = sum(
        row.certainty_class == "certain_change" for row in changed_rows
    )
    uncertain_change_count = sum(
        row.certainty_class == "uncertain_change" for row in changed_rows
    )
    return AncestralTransitionSummary(
        trait=report.trait,
        taxon_column=report.taxon_column,
        model=report.model,
        state_ordering=report.state_ordering,
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(report.exclusions),
        total_branch_count=len(report.branch_rows),
        changed_branch_count=len(changed_rows),
        certain_change_count=certain_change_count,
        uncertain_change_count=uncertain_change_count,
        unique_transition_count=len(report.transition_rows),
        warning_count=len(report.warnings),
    )

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteAncestralReport,
    reconstruct_discrete_ancestral_states,
)
from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.io.newick import loads_newick


@dataclass(frozen=True, slots=True)
class AncestralTransitionBranchRow:
    """One branchwise ancestral transition review row."""

    parent_node: str
    child_node: str
    child_descendant_taxa: list[str]
    branch_length: float | None
    parent_most_likely_state: str
    child_most_likely_state: str
    parent_state_set: list[str]
    child_state_set: list[str]
    overlapping_states: list[str]
    changed: bool
    transition: str
    certainty_class: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionCountRow:
    """One transition-pair count summary for a discrete ancestral reconstruction."""

    transition: str
    source_state: str
    target_state: str
    certain_change_count: int
    uncertain_change_count: int
    total_change_count: int


@dataclass(frozen=True, slots=True)
class AncestralTransitionExclusion:
    """One excluded tip from ancestral transition counting."""

    taxon: str
    reason: str


@dataclass(frozen=True, slots=True)
class AncestralTransitionSummary:
    """Reviewer-facing summary for one ancestral transition count report."""

    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    analyzed_taxon_count: int
    excluded_taxon_count: int
    total_branch_count: int
    changed_branch_count: int
    certain_change_count: int
    uncertain_change_count: int
    unique_transition_count: int
    warning_count: int


@dataclass(slots=True)
class AncestralTransitionReport:
    """Discrete ancestral transition counts for one analyzed rooted tree."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_column: str
    model: str
    state_ordering: str
    ordered_states: list[str]
    taxon_count: int
    branch_rows: list[AncestralTransitionBranchRow]
    transition_rows: list[AncestralTransitionCountRow]
    exclusions: list[AncestralTransitionExclusion]
    warnings: list[str]


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


def write_ancestral_transition_summary_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one summary ledger for ancestral transition counting."""
    summary = summarize_ancestral_transition_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "total_branch_count",
            "changed_branch_count",
            "certain_change_count",
            "uncertain_change_count",
            "unique_transition_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "total_branch_count": str(summary.total_branch_count),
                "changed_branch_count": str(summary.changed_branch_count),
                "certain_change_count": str(summary.certain_change_count),
                "uncertain_change_count": str(summary.uncertain_change_count),
                "unique_transition_count": str(summary.unique_transition_count),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_ancestral_transition_branch_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one branchwise ancestral transition ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_most_likely_state",
            "child_most_likely_state",
            "parent_state_set",
            "child_state_set",
            "overlapping_states",
            "changed",
            "transition",
            "certainty_class",
        ],
        rows=[
            {
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": _format_optional_float(row.branch_length),
                "parent_most_likely_state": row.parent_most_likely_state,
                "child_most_likely_state": row.child_most_likely_state,
                "parent_state_set": ",".join(row.parent_state_set),
                "child_state_set": ",".join(row.child_state_set),
                "overlapping_states": ",".join(row.overlapping_states),
                "changed": str(row.changed).lower(),
                "transition": row.transition,
                "certainty_class": row.certainty_class,
            }
            for row in report.branch_rows
        ],
    )


def write_ancestral_transition_count_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one transition-pair count ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "transition",
            "source_state",
            "target_state",
            "certain_change_count",
            "uncertain_change_count",
            "total_change_count",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_state": row.source_state,
                "target_state": row.target_state,
                "certain_change_count": str(row.certain_change_count),
                "uncertain_change_count": str(row.uncertain_change_count),
                "total_change_count": str(row.total_change_count),
            }
            for row in report.transition_rows
        ],
    )


def write_ancestral_transition_exclusion_table(
    path: Path,
    report: AncestralTransitionReport,
) -> Path:
    """Write one excluded-tip ledger for ancestral transition counting."""
    return write_ancestral_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {
                "taxon": row.taxon,
                "reason": row.reason,
            }
            for row in report.exclusions
        ],
    )


def _build_transition_branch_rows(
    report: DiscreteAncestralReport,
) -> list[AncestralTransitionBranchRow]:
    tree = loads_newick(report.analysis_tree_newick)
    estimate_by_node = {estimate.node: estimate for estimate in report.estimates}
    branch_rows: list[AncestralTransitionBranchRow] = []

    def visit(node: TreeNode) -> None:
        parent_estimate = estimate_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimate_by_node[_node_signature(child)]
            overlapping_states = sorted(
                set(parent_estimate.state_set) & set(child_estimate.state_set)
            )
            changed = (
                parent_estimate.most_likely_state
                != child_estimate.most_likely_state
            )
            branch_rows.append(
                AncestralTransitionBranchRow(
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=child_estimate.descendant_taxa,
                    branch_length=child.branch_length,
                    parent_most_likely_state=parent_estimate.most_likely_state,
                    child_most_likely_state=child_estimate.most_likely_state,
                    parent_state_set=parent_estimate.state_set,
                    child_state_set=child_estimate.state_set,
                    overlapping_states=overlapping_states,
                    changed=changed,
                    transition=(
                        f"{parent_estimate.most_likely_state}->{child_estimate.most_likely_state}"
                        if changed
                        else ""
                    ),
                    certainty_class=_transition_certainty_class(
                        changed=changed,
                        overlapping_states=overlapping_states,
                        parent_state_set=parent_estimate.state_set,
                        child_state_set=child_estimate.state_set,
                    ),
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def _summarize_transition_rows(
    branch_rows: list[AncestralTransitionBranchRow],
) -> list[AncestralTransitionCountRow]:
    grouped: dict[str, list[AncestralTransitionBranchRow]] = {}
    for row in branch_rows:
        if row.changed:
            grouped.setdefault(row.transition, []).append(row)
    transition_rows: list[AncestralTransitionCountRow] = []
    for transition in sorted(grouped):
        transition_branch_rows = grouped[transition]
        source_state, target_state = transition.split("->", maxsplit=1)
        certain_change_count = sum(
            row.certainty_class == "certain_change" for row in transition_branch_rows
        )
        uncertain_change_count = sum(
            row.certainty_class == "uncertain_change"
            for row in transition_branch_rows
        )
        transition_rows.append(
            AncestralTransitionCountRow(
                transition=transition,
                source_state=source_state,
                target_state=target_state,
                certain_change_count=certain_change_count,
                uncertain_change_count=uncertain_change_count,
                total_change_count=len(transition_branch_rows),
            )
        )
    return transition_rows


def _transition_certainty_class(
    *,
    changed: bool,
    overlapping_states: list[str],
    parent_state_set: list[str],
    child_state_set: list[str],
) -> str:
    if changed:
        if not overlapping_states:
            return "certain_change"
        return "uncertain_change"
    if len(parent_state_set) == 1 and len(child_state_set) == 1:
        return "certain_no_change"
    return "uncertain_no_change"


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    descendant_taxa: list[str] = []
    for child in node.children:
        descendant_taxa.extend(_node_signature_taxa(child))
    return "|".join(sorted(descendant_taxa))


def _node_signature_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    descendant_taxa: list[str] = []
    for child in node.children:
        descendant_taxa.extend(_node_signature_taxa(child))
    return descendant_taxa


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return format(value, ".15g")

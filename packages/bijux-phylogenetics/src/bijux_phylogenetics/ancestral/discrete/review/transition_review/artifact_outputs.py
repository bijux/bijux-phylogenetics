from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.ancestral.common import write_ancestral_rows

from .aggregation import _format_optional_float
from .contracts import AncestralTransitionReport, AncestralTransitionTreeSetReport
from .single_tree_review import summarize_ancestral_transition_report
from .tree_set_review import summarize_ancestral_transition_tree_set_report


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
    report: AncestralTransitionReport | AncestralTransitionTreeSetReport,
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


def write_ancestral_transition_tree_set_summary_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one summary ledger for ancestral transition counting across a tree set."""
    summary = summarize_ancestral_transition_tree_set_report(report)
    return write_ancestral_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "model",
            "state_ordering",
            "total_tree_count",
            "burnin_tree_count",
            "kept_tree_count",
            "shared_tree_taxon_count",
            "analysis_taxon_count",
            "rooted_topology_count",
            "unrooted_topology_count",
            "transition_pair_count",
            "topology_sensitive_transition_pair_count",
            "uncertainty_sensitive_transition_pair_count",
            "warning_count",
        ],
        rows=[
            {
                "trait": summary.trait,
                "taxon_column": summary.taxon_column,
                "model": summary.model,
                "state_ordering": summary.state_ordering,
                "total_tree_count": str(summary.total_tree_count),
                "burnin_tree_count": str(summary.burnin_tree_count),
                "kept_tree_count": str(summary.kept_tree_count),
                "shared_tree_taxon_count": str(summary.shared_tree_taxon_count),
                "analysis_taxon_count": str(summary.analysis_taxon_count),
                "rooted_topology_count": str(summary.rooted_topology_count),
                "unrooted_topology_count": str(summary.unrooted_topology_count),
                "transition_pair_count": str(summary.transition_pair_count),
                "topology_sensitive_transition_pair_count": str(
                    summary.topology_sensitive_transition_pair_count
                ),
                "uncertainty_sensitive_transition_pair_count": str(
                    summary.uncertainty_sensitive_transition_pair_count
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_ancestral_transition_tree_set_tree_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one retained-tree transition summary ledger."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
            "total_branch_count",
            "changed_branch_count",
            "certain_change_count",
            "uncertain_change_count",
        ],
        rows=[
            {
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
                "total_branch_count": str(row.total_branch_count),
                "changed_branch_count": str(row.changed_branch_count),
                "certain_change_count": str(row.certain_change_count),
                "uncertain_change_count": str(row.uncertain_change_count),
            }
            for row in report.tree_rows
        ],
    )


def write_ancestral_transition_tree_set_branch_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one per-branch ancestral transition ledger across retained trees."""
    return write_ancestral_rows(
        path,
        columns=[
            "source_tree_index",
            "post_burnin_index",
            "rooted_topology_id",
            "unrooted_topology_id",
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
                "source_tree_index": str(row.source_tree_index),
                "post_burnin_index": str(row.post_burnin_index),
                "rooted_topology_id": row.rooted_topology_id,
                "unrooted_topology_id": row.unrooted_topology_id,
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


def write_ancestral_transition_tree_set_count_table(
    path: Path,
    report: AncestralTransitionTreeSetReport,
) -> Path:
    """Write one transition-pair summary ledger across retained trees."""
    return write_ancestral_rows(
        path,
        columns=[
            "transition",
            "source_state",
            "target_state",
            "tree_presence_count",
            "tree_presence_fraction",
            "mean_certain_change_count",
            "mean_uncertain_change_count",
            "mean_total_change_count",
            "minimum_total_change_count",
            "maximum_total_change_count",
            "lower_95_total_change_count",
            "upper_95_total_change_count",
            "stability_class",
        ],
        rows=[
            {
                "transition": row.transition,
                "source_state": row.source_state,
                "target_state": row.target_state,
                "tree_presence_count": str(row.tree_presence_count),
                "tree_presence_fraction": str(row.tree_presence_fraction),
                "mean_certain_change_count": str(row.mean_certain_change_count),
                "mean_uncertain_change_count": str(row.mean_uncertain_change_count),
                "mean_total_change_count": str(row.mean_total_change_count),
                "minimum_total_change_count": str(row.minimum_total_change_count),
                "maximum_total_change_count": str(row.maximum_total_change_count),
                "lower_95_total_change_count": str(row.lower_95_total_change_count),
                "upper_95_total_change_count": str(row.upper_95_total_change_count),
                "stability_class": row.stability_class,
            }
            for row in report.transition_rows
        ],
    )

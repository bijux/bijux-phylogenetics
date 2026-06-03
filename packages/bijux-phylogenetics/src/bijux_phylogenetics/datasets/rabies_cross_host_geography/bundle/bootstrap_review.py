from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.compare.presentation import (
    ComparisonReportBuildResult,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.trees import BootstrapTreeSetSummaryReport

from ..shared import _format_number


def _write_bootstrap_tree_comparison_summary(
    path: Path,
    comparison_report: ComparisonReportBuildResult,
) -> Path:
    high_support_conflict_count = len(
        [
            row
            for row in comparison_report.support.conflicting_clades
            if row.conflict_classification == "high_support_conflict"
        ]
    )
    row = {
        "left_tree": comparison_report.topology.left_path.name,
        "right_tree": comparison_report.topology.right_path.name,
        "shared_taxon_count": str(len(comparison_report.topology.shared_taxa)),
        "rooted_rf_distance": str(
            comparison_report.topology.rooted_robinson_foulds_distance
        ),
        "rooted_normalized_rf": _format_number(
            comparison_report.topology.rooted_normalized_robinson_foulds
        ),
        "topology_equal": "true"
        if comparison_report.topology.topology_equal
        else "false",
        "same_unrooted_topology": (
            "true" if comparison_report.topology.same_unrooted_topology else "false"
        ),
        "same_taxa_different_rooting": (
            "true"
            if comparison_report.topology.same_taxa_different_rooting
            else "false"
        ),
        "same_topology_different_branch_lengths": (
            "true"
            if comparison_report.topology.same_topology_different_branch_lengths
            else "false"
        ),
        "support_conflict_count": str(
            len(comparison_report.support.conflicting_clades)
        ),
        "high_support_conflict_count": str(high_support_conflict_count),
        "branch_score_distance": _format_number(
            comparison_report.branch_lengths.branch_score.branch_score_distance
        ),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_stable_bootstrap_summary_table(
    path: Path,
    report: BootstrapTreeSetSummaryReport,
) -> Path:
    row = {
        "tree_count": str(report.tree_count),
        "shared_taxon_count": str(len(report.shared_taxa)),
        "rooted_topology_count": str(report.diversity.rooted_topology_count),
        "dominant_topology_frequency": _format_number(
            report.diversity.dominant_topology_frequency
        ),
        "effective_topology_count": _format_number(
            report.diversity.effective_topology_count
        ),
        "mean_robinson_foulds_distance": _format_number(
            report.diversity.mean_robinson_foulds_distance
        ),
        "mean_normalized_robinson_foulds_distance": _format_number(
            report.diversity.mean_normalized_robinson_foulds_distance
        ),
        "consensus_threshold": _format_number(report.consensus_threshold),
        "robust_support_threshold": _format_number(report.robust_support_threshold),
        "unstable_branch_count": str(report.unstable_branch_count),
        "warning_count": str(len(report.warnings)),
        "consensus_newick": report.consensus.consensus_newick,
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _stabilize_bundle_report_paths(path: Path, *, output_root: Path) -> Path:
    text = path.read_text(encoding="utf-8")
    normalized_root = output_root.as_posix().rstrip("/") + "/"
    path.write_text(text.replace(normalized_root, ""), encoding="utf-8")
    return path

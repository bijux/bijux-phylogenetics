from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.presentation import (
    ComparisonReportBuildResult,
    build_tree_comparison_report,
)
from bijux_phylogenetics.compare.topology import write_tree_comparison_table
from bijux_phylogenetics.datasets.rabies_cross_host_geography.models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
)
from bijux_phylogenetics.trees import (
    BootstrapTreeSetArtifactReport,
    write_bootstrap_tree_set_artifacts,
)

from .bootstrap_review import (
    _stabilize_bundle_report_paths,
    _write_bootstrap_tree_comparison_summary,
    _write_stable_bootstrap_summary_table,
)


@dataclass(frozen=True)
class BootstrapReviewArtifacts:
    output_root: Path
    artifact_report: BootstrapTreeSetArtifactReport
    summary_path: Path
    comparison_report: ComparisonReportBuildResult
    comparison_table_path: Path
    comparison_summary_path: Path


def _write_bootstrap_review_artifacts(
    output_root: Path,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    *,
    rooted_tree_path: Path,
) -> BootstrapReviewArtifacts:
    bootstrap_output_root = output_root / "bootstrap-review"
    bootstrap_artifacts = write_bootstrap_tree_set_artifacts(
        report.fasta_to_tree.bootstrap_workflow.output_paths["bootstrap_trees"],
        out_dir=bootstrap_output_root,
        prefix="bootstrap-review",
        consensus_threshold=report.config.bootstrap_consensus_threshold,
        robust_support_threshold=report.config.bootstrap_robust_support_threshold,
        max_tree_count=report.config.max_bootstrap_tree_count,
        memory_warning_threshold_bytes=report.config.memory_warning_threshold_bytes,
    )
    bootstrap_summary_path = _write_stable_bootstrap_summary_table(
        bootstrap_output_root / "bootstrap-review.summary.tsv",
        bootstrap_artifacts.summary_report,
    )
    bootstrap_tree_comparison_report = build_tree_comparison_report(
        rooted_tree_path,
        bootstrap_artifacts.output_paths["consensus_tree"],
        out_path=bootstrap_output_root
        / "rooted-tree-vs-bootstrap-consensus.report.html",
    )
    _stabilize_bundle_report_paths(
        bootstrap_tree_comparison_report.output_path,
        output_root=output_root,
    )
    bootstrap_tree_comparison_table_path = write_tree_comparison_table(
        bootstrap_output_root / "rooted-tree-vs-bootstrap-consensus.comparison.tsv",
        rooted_tree_path,
        bootstrap_artifacts.output_paths["consensus_tree"],
    )
    bootstrap_tree_comparison_summary_path = _write_bootstrap_tree_comparison_summary(
        bootstrap_output_root / "rooted-tree-vs-bootstrap-consensus.summary.tsv",
        bootstrap_tree_comparison_report,
    )
    return BootstrapReviewArtifacts(
        output_root=bootstrap_output_root,
        artifact_report=bootstrap_artifacts,
        summary_path=bootstrap_summary_path,
        comparison_report=bootstrap_tree_comparison_report,
        comparison_table_path=bootstrap_tree_comparison_table_path,
        comparison_summary_path=bootstrap_tree_comparison_summary_path,
    )

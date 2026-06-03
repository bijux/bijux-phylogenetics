from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
)
from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.trees import BootstrapTreeSetArtifactReport

from ..models import RabiesCrossHostGeographyPanelWorkflowReport
from ..shared import _checksum, _format_number


def _write_workflow_summary_table(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    clade_row_count: int,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    scientific_finding_count: int,
) -> Path:
    support = report.fasta_to_tree.support_summary
    host_summary = report.host_switching.summary
    geography_summary = report.biogeography_report.state_report.summary
    migration_summary = report.biogeography_report.event_report.summary
    bootstrap_summary = bootstrap_artifacts.summary_report
    row = {
        "dataset_id": report.dataset.dataset_id,
        "sequence_count": str(report.dataset.sequence_count),
        "sequence_type": report.dataset.sequence_type,
        "inferred_sequence_type": report.fasta_to_tree.sequence_type,
        "selected_model": report.fasta_to_tree.selected_model,
        "input_repair_applied": (
            "true" if report.fasta_to_tree.input_repair is not None else "false"
        ),
        "aligned_quality_score": _format_number(report.aligned_quality.quality_score),
        "trimmed_quality_score": _format_number(report.trimmed_quality.quality_score),
        "minimum_support": _format_number(support.minimum_support),
        "maximum_support": _format_number(support.maximum_support),
        "median_support": _format_number(support.median_support),
        "weakly_supported_clade_count": str(support.weakly_supported_clade_count),
        "clade_row_count": str(clade_row_count),
        "bootstrap_tree_count": str(bootstrap_summary.tree_count),
        "bootstrap_topology_count": str(
            bootstrap_summary.diversity.rooted_topology_count
        ),
        "bootstrap_unstable_branch_count": str(bootstrap_summary.unstable_branch_count),
        "bootstrap_consensus_rooted_rf_distance": str(
            bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
        ),
        "bootstrap_consensus_same_unrooted_topology": (
            "true"
            if bootstrap_tree_comparison_report.topology.same_unrooted_topology
            else "false"
        ),
        "bootstrap_consensus_high_support_conflict_count": str(
            len(
                [
                    row
                    for row in bootstrap_tree_comparison_report.support.conflicting_clades
                    if row.conflict_classification == "high_support_conflict"
                ]
            )
        ),
        "outgroup_taxa": ",".join(report.dataset.outgroup_taxa),
        "root_host": host_summary.root_host,
        "root_host_confidence": _format_number(host_summary.root_confidence),
        "host_switch_count": str(host_summary.host_switch_count),
        "certain_host_switch_count": str(host_summary.certain_host_switch_count),
        "uncertain_host_switch_count": str(host_summary.uncertain_host_switch_count),
        "root_region": geography_summary.root_region,
        "root_region_probability": _format_number(
            geography_summary.root_region_probability
        ),
        "changed_region_branch_count": str(geography_summary.changed_branch_count),
        "migration_event_count": str(migration_summary.event_count),
        "strongly_supported_migration_event_count": str(
            migration_summary.strongly_supported_event_count
        ),
        "comparative_response": comparative_summary_row.response,
        "comparative_formula": comparative_summary_row.formula,
        "comparative_selected_model": comparative_summary_row.selected_model,
        "comparative_pgls_lambda": _format_number(comparative_summary_row.pgls_lambda),
        "comparative_pgls_r_squared": _format_number(
            comparative_summary_row.pgls_r_squared
        ),
        "comparative_branch_repair_count": str(len(report.comparative_branch_repairs)),
        "conclusion_stable_count": str(
            report.conclusion_stability_report.summary.stable_count
        ),
        "conclusion_weak_count": str(
            report.conclusion_stability_report.summary.weak_count
        ),
        "conclusion_unstable_count": str(
            report.conclusion_stability_report.summary.unstable_count
        ),
        "timeout_seconds": _format_number(report.config.timeout_seconds),
        "max_bootstrap_tree_count": (
            ""
            if report.config.max_bootstrap_tree_count is None
            else str(report.config.max_bootstrap_tree_count)
        ),
        "max_report_table_rows": (
            ""
            if report.config.max_report_table_rows is None
            else str(report.config.max_report_table_rows)
        ),
        "memory_warning_threshold_bytes": (
            ""
            if report.config.memory_warning_threshold_bytes is None
            else str(report.config.memory_warning_threshold_bytes)
        ),
        "config_check_count": str(len(report.config_audit_rows)),
        "scientific_finding_count": str(scientific_finding_count),
    }
    return write_taxon_rows(path, columns=list(row.keys()), rows=[row])


def _write_resource_observation_table(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
) -> Path:
    rows = [
        {
            "resource_surface": "workflow",
            "max_bootstrap_tree_count": (
                ""
                if report.config.max_bootstrap_tree_count is None
                else str(report.config.max_bootstrap_tree_count)
            ),
            "memory_warning_threshold_bytes": (
                ""
                if report.config.memory_warning_threshold_bytes is None
                else str(report.config.memory_warning_threshold_bytes)
            ),
            "tree_count_observed": str(bootstrap_artifacts.summary_report.tree_count),
            "warning_messages": " | ".join(
                bootstrap_artifacts.budget_report.warning_messages
            ),
            "error_messages": " | ".join(
                bootstrap_artifacts.budget_report.error_messages
            ),
        }
    ]
    return write_taxon_rows(path, columns=list(rows[0].keys()), rows=rows)


def _write_manifest(
    path: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    clade_row_count: int,
    scientific_finding_count: int,
    bundle_paths: dict[str, Path],
) -> Path:
    manifest = {
        "report_kind": "rabies_cross_host_geography_workflow_bundle",
        "dataset_id": report.dataset.dataset_id,
        "input_checksums": {
            "workflow-config.json": _checksum(report.dataset.workflow_config_path),
            "sequences.fasta": _checksum(report.dataset.sequences_path),
            "metadata.csv": _checksum(report.dataset.metadata_path),
            "region-centroids.csv": _checksum(report.dataset.centroids_path),
        },
        "output_checksums": {
            key: _checksum(value) for key, value in bundle_paths.items()
        },
        "metrics": {
            "sequence_count": report.dataset.sequence_count,
            "selected_model": report.fasta_to_tree.selected_model,
            "minimum_support": report.fasta_to_tree.support_summary.minimum_support,
            "maximum_support": report.fasta_to_tree.support_summary.maximum_support,
            "host_switch_count": report.host_switching.summary.host_switch_count,
            "migration_event_count": report.biogeography_report.event_report.summary.event_count,
            "root_host": report.host_switching.summary.root_host,
            "root_region": report.biogeography_report.state_report.summary.root_region,
            "clade_row_count": clade_row_count,
            "bootstrap_tree_count": bootstrap_artifacts.summary_report.tree_count,
            "bootstrap_unstable_branch_count": (
                bootstrap_artifacts.summary_report.unstable_branch_count
            ),
            "bootstrap_consensus_rooted_rf_distance": (
                bootstrap_tree_comparison_report.topology.rooted_robinson_foulds_distance
            ),
            "bootstrap_consensus_same_unrooted_topology": (
                bootstrap_tree_comparison_report.topology.same_unrooted_topology
            ),
            "bootstrap_consensus_high_support_conflict_count": len(
                [
                    row
                    for row in bootstrap_tree_comparison_report.support.conflicting_clades
                    if row.conflict_classification == "high_support_conflict"
                ]
            ),
            "comparative_selected_model": comparative_summary_row.selected_model,
            "comparative_pgls_lambda": comparative_summary_row.pgls_lambda,
            "comparative_pgls_r_squared": comparative_summary_row.pgls_r_squared,
            "conclusion_stable_count": (
                report.conclusion_stability_report.summary.stable_count
            ),
            "conclusion_weak_count": report.conclusion_stability_report.summary.weak_count,
            "conclusion_unstable_count": (
                report.conclusion_stability_report.summary.unstable_count
            ),
            "config_check_count": len(report.config_audit_rows),
            "scientific_finding_count": scientific_finding_count,
        },
    }
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return path

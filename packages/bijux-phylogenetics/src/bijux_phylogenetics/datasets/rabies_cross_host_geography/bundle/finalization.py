from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.comparative.reporting.analysis_package import (
    ComparativeAnalysisSummaryRow,
    ComparativeInterpretationRow,
)
from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.datasets.rabies_cross_host_geography.models import (
    RabiesCrossHostGeographyPanelWorkflowReport,
    RabiesScientificFindingRow,
)
from bijux_phylogenetics.trees import BootstrapTreeSetArtifactReport

from .findings import _build_scientific_finding_rows, _write_scientific_findings_table
from .integrated_report import _write_integrated_report
from .package_ledger import (
    _write_manifest,
    _write_resource_observation_table,
    _write_workflow_summary_table,
)


@dataclass(frozen=True)
class FinalBundleArtifacts:
    scientific_finding_rows: list[RabiesScientificFindingRow]
    scientific_findings_path: Path
    workflow_summary_path: Path
    resource_observations_path: Path
    final_report_path: Path
    final_manifest_path: Path


def _write_final_bundle_artifacts(
    output_root: Path,
    *,
    report: RabiesCrossHostGeographyPanelWorkflowReport,
    clade_row_count: int,
    bootstrap_artifacts: BootstrapTreeSetArtifactReport,
    bootstrap_tree_comparison_report: ComparisonReportBuildResult,
    comparative_summary_row: ComparativeAnalysisSummaryRow,
    comparative_interpretation_rows: list[ComparativeInterpretationRow],
    comparative_branch_repair_count: int,
    bundle_paths: dict[str, Path],
) -> FinalBundleArtifacts:
    scientific_finding_rows = _build_scientific_finding_rows(
        report=report,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        comparative_summary_row=comparative_summary_row,
        comparative_interpretation_rows=comparative_interpretation_rows,
    )
    scientific_findings_path = _write_scientific_findings_table(
        output_root / "scientific-findings.tsv",
        scientific_finding_rows,
    )
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report=report,
        clade_row_count=clade_row_count,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        comparative_summary_row=comparative_summary_row,
        scientific_finding_count=len(scientific_finding_rows),
    )
    resource_observations_path = _write_resource_observation_table(
        output_root / "resource-observations.tsv",
        report=report,
        bootstrap_artifacts=bootstrap_artifacts,
    )
    final_report_path = _write_integrated_report(
        output_root / "rabies-cross-host-geography-report.html",
        report=report,
        workflow_summary_path=workflow_summary_path,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        clade_row_count=clade_row_count,
        comparative_summary_row=comparative_summary_row,
        comparative_interpretation_rows=comparative_interpretation_rows,
        comparative_branch_repair_count=comparative_branch_repair_count,
        scientific_finding_rows=scientific_finding_rows,
        max_report_table_rows=report.config.max_report_table_rows,
    )
    final_manifest_path = _write_manifest(
        output_root / "rabies-cross-host-geography.manifest.json",
        report=report,
        comparative_summary_row=comparative_summary_row,
        bootstrap_artifacts=bootstrap_artifacts,
        bootstrap_tree_comparison_report=bootstrap_tree_comparison_report,
        clade_row_count=clade_row_count,
        scientific_finding_count=len(scientific_finding_rows),
        bundle_paths={
            **bundle_paths,
            "workflow_summary": workflow_summary_path,
            "scientific_findings": scientific_findings_path,
            "final_report": final_report_path,
        },
    )
    return FinalBundleArtifacts(
        scientific_finding_rows=scientific_finding_rows,
        scientific_findings_path=scientific_findings_path,
        workflow_summary_path=workflow_summary_path,
        resource_observations_path=resource_observations_path,
        final_report_path=final_report_path,
        final_manifest_path=final_manifest_path,
    )

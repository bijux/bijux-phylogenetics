from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..models import RabiesMethodSensitivityPanelWorkflowReport
from .package_ledger import (
    _write_clade_table,
    _write_conclusion_summary_table,
    _write_parallel_execution_summary_table,
    _write_preprocessing_comparison_table,
    _write_variant_summary_table,
    _write_workflow_summary_table,
)
from .package_manifest import _write_resolved_config
from .variant_artifacts import _copy_output, _copy_task_logs, _write_variant_outputs


@dataclass(slots=True)
class RabiesMethodSensitivityWorkflowBundleArtifacts:
    """Owned workflow artifact set written into the reviewer-facing bundle."""

    workflow_summary_path: Path
    variant_summary_path: Path
    parallel_summary_path: Path
    preprocessing_comparison_path: Path
    stable_clades_path: Path
    changed_clades_path: Path
    conclusion_summary_path: Path
    config_path: Path
    execution_record_path: Path
    task_logs_root: Path
    variants_root: Path


def _write_workflow_bundle_artifacts(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
) -> RabiesMethodSensitivityWorkflowBundleArtifacts:
    workflow_summary_path = _write_workflow_summary_table(
        output_root / "workflow-summary.tsv",
        report,
    )
    variant_summary_path = _write_variant_summary_table(
        output_root / "variant-summary.tsv",
        report,
    )
    parallel_summary_path = _write_parallel_execution_summary_table(
        output_root / "parallel-execution-summary.tsv",
        report,
    )
    preprocessing_comparison_path = _write_preprocessing_comparison_table(
        output_root / "preprocessing-rooted-comparisons.tsv",
        report.preprocessing_comparison_rows,
    )
    stable_clades_path = _write_clade_table(
        output_root / "stable-clades.tsv",
        report.stable_clade_rows,
    )
    changed_clades_path = _write_clade_table(
        output_root / "changed-clades.tsv",
        report.changed_clade_rows,
    )
    conclusion_summary_path = _write_conclusion_summary_table(
        output_root / "method-conclusion-summary.tsv",
        report.conclusion_rows,
    )
    config_path = _write_resolved_config(
        output_root / "workflow-config.resolved.json",
        report,
    )
    execution_record_path = _copy_output(
        report.execution_record_path,
        output_root / report.execution_record_path.name,
    )
    task_logs_root = _copy_task_logs(output_root / "parallel-logs", report.task_records)
    variants_root = _write_variant_outputs(
        output_root / "variants",
        report.variant_runs,
    )
    return RabiesMethodSensitivityWorkflowBundleArtifacts(
        workflow_summary_path=workflow_summary_path,
        variant_summary_path=variant_summary_path,
        parallel_summary_path=parallel_summary_path,
        preprocessing_comparison_path=preprocessing_comparison_path,
        stable_clades_path=stable_clades_path,
        changed_clades_path=changed_clades_path,
        conclusion_summary_path=conclusion_summary_path,
        config_path=config_path,
        execution_record_path=execution_record_path,
        task_logs_root=task_logs_root,
        variants_root=variants_root,
    )

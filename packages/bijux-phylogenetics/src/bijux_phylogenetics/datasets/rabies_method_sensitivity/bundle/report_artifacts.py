from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..models import RabiesMethodSensitivityPanelWorkflowReport
from ..reporting import _write_report
from .manifest_artifacts import RabiesMethodSensitivityBundleManifestArtifacts
from .reproducibility_artifacts import (
    RabiesMethodSensitivityReproducibilityArtifacts,
)
from .slurm_artifacts import RabiesMethodSensitivitySlurmBundleArtifacts
from .workflow_artifacts import RabiesMethodSensitivityWorkflowBundleArtifacts


@dataclass(slots=True)
class RabiesMethodSensitivityBundleReportArtifacts:
    """Owned HTML report artifacts for the reviewer-facing workflow bundle."""

    report_path: Path
    linked_files: tuple[Path, ...]
    linked_artifact_count: int
    html_size_bytes: int
    linked_artifact_bytes: int
    total_output_bytes: int


def _report_bundle_paths(
    workflow_artifacts: RabiesMethodSensitivityWorkflowBundleArtifacts,
    manifest_artifacts: RabiesMethodSensitivityBundleManifestArtifacts,
    slurm_artifacts: RabiesMethodSensitivitySlurmBundleArtifacts,
    reproducibility_artifacts: RabiesMethodSensitivityReproducibilityArtifacts,
) -> dict[str, Path]:
    return {
        "workflow_summary": workflow_artifacts.workflow_summary_path,
        "variant_summary": workflow_artifacts.variant_summary_path,
        "parallel_summary": workflow_artifacts.parallel_summary_path,
        "preprocessing_comparison": (workflow_artifacts.preprocessing_comparison_path),
        "stable_clades": workflow_artifacts.stable_clades_path,
        "changed_clades": workflow_artifacts.changed_clades_path,
        "conclusion_summary": workflow_artifacts.conclusion_summary_path,
        "config": workflow_artifacts.config_path,
        "execution_record": workflow_artifacts.execution_record_path,
        "workflow_manifest": manifest_artifacts.manifest_path,
        "slurm_job_plan": slurm_artifacts.job_plan_path,
        "slurm_assumptions": slurm_artifacts.assumptions_path,
        "slurm_summary": slurm_artifacts.summary_path,
        "slurm_array_partitions": slurm_artifacts.array_partitions_path,
        "slurm_array_members": slurm_artifacts.array_members_path,
        "slurm_array_strategy": slurm_artifacts.array_strategy_path,
        "slurm_job_evidence_index": slurm_artifacts.job_evidence_report.index_path,
        "slurm_job_evidence_summary": slurm_artifacts.job_evidence_report.summary_path,
        "slurm_storage_categories": slurm_artifacts.storage_categories_path,
        "slurm_storage_variants": slurm_artifacts.storage_variants_path,
        "slurm_storage_summary": slurm_artifacts.storage_summary_path,
        "slurm_storage_report": slurm_artifacts.storage_report_path,
        "slurm_output_explosion_checks": slurm_artifacts.output_explosion_checks_path,
        "slurm_output_explosion_variants": (
            slurm_artifacts.output_explosion_variants_path
        ),
        "slurm_output_explosion_summary": (
            slurm_artifacts.output_explosion_summary_path
        ),
        "slurm_output_explosion_report": slurm_artifacts.output_explosion_report_path,
        "slurm_tree_retention_checks": slurm_artifacts.tree_retention_checks_path,
        "slurm_tree_retention_files": slurm_artifacts.tree_retention_files_path,
        "slurm_tree_retention_summary": slurm_artifacts.tree_retention_summary_path,
        "slurm_tree_retention_report": slurm_artifacts.tree_retention_report_path,
        "slurm_merge_checks": slurm_artifacts.merge_checks_path,
        "slurm_merge_variants": slurm_artifacts.merge_variants_path,
        "slurm_merge_summary": slurm_artifacts.merge_summary_path,
        "slurm_merge_report": slurm_artifacts.merge_report_path,
        "slurm_output_freshness": slurm_artifacts.output_freshness_path,
        "slurm_output_freshness_checks": slurm_artifacts.output_freshness_checks_path,
        "slurm_output_freshness_summary": (
            slurm_artifacts.output_freshness_summary_path
        ),
        "slurm_job_status": slurm_artifacts.job_status_path,
        "slurm_partition_status": slurm_artifacts.partition_status_path,
        "slurm_workflow_status": slurm_artifacts.workflow_status_path,
        "slurm_failure_recovery_jobs": slurm_artifacts.failure_recovery_jobs_path,
        "slurm_failure_recovery_partitions": (
            slurm_artifacts.failure_recovery_partitions_path
        ),
        "slurm_failure_recovery_summary": slurm_artifacts.failure_recovery_summary_path,
        "slurm_failure_recovery_report": slurm_artifacts.failure_recovery_report_path,
        "reproducibility_checks": reproducibility_artifacts.checks_path,
        "reproducibility_variant_audit": reproducibility_artifacts.variant_audit_path,
        "reproducibility_audit": reproducibility_artifacts.audit_path,
    }


def _report_linked_files(
    workflow_artifacts: RabiesMethodSensitivityWorkflowBundleArtifacts,
    manifest_artifacts: RabiesMethodSensitivityBundleManifestArtifacts,
    slurm_artifacts: RabiesMethodSensitivitySlurmBundleArtifacts,
    reproducibility_artifacts: RabiesMethodSensitivityReproducibilityArtifacts,
) -> tuple[Path, ...]:
    return (
        workflow_artifacts.workflow_summary_path,
        workflow_artifacts.variant_summary_path,
        workflow_artifacts.parallel_summary_path,
        workflow_artifacts.preprocessing_comparison_path,
        workflow_artifacts.stable_clades_path,
        workflow_artifacts.changed_clades_path,
        workflow_artifacts.conclusion_summary_path,
        workflow_artifacts.config_path,
        workflow_artifacts.execution_record_path,
        manifest_artifacts.manifest_path,
        manifest_artifacts.report_manifest_path,
        slurm_artifacts.job_plan_path,
        slurm_artifacts.assumptions_path,
        slurm_artifacts.summary_path,
        slurm_artifacts.array_partitions_path,
        slurm_artifacts.array_members_path,
        slurm_artifacts.array_strategy_path,
        slurm_artifacts.job_evidence_report.index_path,
        slurm_artifacts.job_evidence_report.summary_path,
        slurm_artifacts.storage_categories_path,
        slurm_artifacts.storage_variants_path,
        slurm_artifacts.storage_summary_path,
        slurm_artifacts.storage_report_path,
        slurm_artifacts.output_explosion_checks_path,
        slurm_artifacts.output_explosion_variants_path,
        slurm_artifacts.output_explosion_summary_path,
        slurm_artifacts.output_explosion_report_path,
        slurm_artifacts.tree_retention_checks_path,
        slurm_artifacts.tree_retention_files_path,
        slurm_artifacts.tree_retention_summary_path,
        slurm_artifacts.tree_retention_report_path,
        slurm_artifacts.merge_checks_path,
        slurm_artifacts.merge_variants_path,
        slurm_artifacts.merge_summary_path,
        slurm_artifacts.merge_report_path,
        slurm_artifacts.output_freshness_path,
        slurm_artifacts.output_freshness_checks_path,
        slurm_artifacts.output_freshness_summary_path,
        slurm_artifacts.job_status_path,
        slurm_artifacts.partition_status_path,
        slurm_artifacts.workflow_status_path,
        slurm_artifacts.failure_recovery_jobs_path,
        slurm_artifacts.failure_recovery_partitions_path,
        slurm_artifacts.failure_recovery_summary_path,
        slurm_artifacts.failure_recovery_report_path,
        reproducibility_artifacts.checks_path,
        reproducibility_artifacts.variant_audit_path,
        reproducibility_artifacts.audit_path,
    )


def _write_bundle_report_artifacts(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    workflow_artifacts: RabiesMethodSensitivityWorkflowBundleArtifacts,
    manifest_artifacts: RabiesMethodSensitivityBundleManifestArtifacts,
    slurm_artifacts: RabiesMethodSensitivitySlurmBundleArtifacts,
    reproducibility_artifacts: RabiesMethodSensitivityReproducibilityArtifacts,
) -> RabiesMethodSensitivityBundleReportArtifacts:
    linked_files = _report_linked_files(
        workflow_artifacts,
        manifest_artifacts,
        slurm_artifacts,
        reproducibility_artifacts,
    )
    report_path = _write_report(
        output_root / "rabies-method-sensitivity-report.html",
        report=report,
        bundle_paths=_report_bundle_paths(
            workflow_artifacts,
            manifest_artifacts,
            slurm_artifacts,
            reproducibility_artifacts,
        ),
        report_manifest_path=manifest_artifacts.report_manifest_path,
        reproducibility_report=reproducibility_artifacts.report,
        slurm_planning_report=slurm_artifacts.planning_report,
        slurm_array_strategy_report=slurm_artifacts.array_strategy_report,
        slurm_job_evidence_report=slurm_artifacts.job_evidence_report,
        slurm_storage_report=slurm_artifacts.storage_report,
        slurm_output_explosion_report=slurm_artifacts.output_explosion_report,
        slurm_tree_retention_report=slurm_artifacts.tree_retention_report,
        slurm_merge_report=slurm_artifacts.merge_report,
        slurm_output_freshness_report=slurm_artifacts.output_freshness_report,
        slurm_status_report=slurm_artifacts.status_report,
        slurm_failure_recovery_report=slurm_artifacts.failure_recovery_report,
    )
    html_size_bytes = report_path.stat().st_size
    linked_artifact_bytes = sum(path.stat().st_size for path in linked_files)
    total_output_bytes = html_size_bytes + linked_artifact_bytes
    return RabiesMethodSensitivityBundleReportArtifacts(
        report_path=report_path,
        linked_files=linked_files,
        linked_artifact_count=len(linked_files),
        html_size_bytes=html_size_bytes,
        linked_artifact_bytes=linked_artifact_bytes,
        total_output_bytes=total_output_bytes,
    )

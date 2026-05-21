from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..models import RabiesMethodSensitivityPanelWorkflowReport
from ..reporting import _write_report_manifest
from .package_manifest import _sha256, _write_manifest
from .slurm_artifacts import RabiesMethodSensitivitySlurmBundleArtifacts
from .workflow_artifacts import RabiesMethodSensitivityWorkflowBundleArtifacts


@dataclass(slots=True)
class RabiesMethodSensitivityBundleManifestArtifacts:
    """Owned manifest artifacts for the reviewer-facing workflow bundle."""

    manifest_path: Path
    report_manifest_path: Path


def _workflow_manifest_paths(
    workflow_artifacts: RabiesMethodSensitivityWorkflowBundleArtifacts,
    slurm_artifacts: RabiesMethodSensitivitySlurmBundleArtifacts,
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
        "task_logs_root": workflow_artifacts.task_logs_root,
        "variants_root": workflow_artifacts.variants_root,
        "slurm_job_plan": slurm_artifacts.job_plan_path,
        "slurm_assumptions": slurm_artifacts.assumptions_path,
        "slurm_summary": slurm_artifacts.summary_path,
        "slurm_array_partitions": slurm_artifacts.array_partitions_path,
        "slurm_array_members": slurm_artifacts.array_members_path,
        "slurm_array_strategy": slurm_artifacts.array_strategy_path,
        "slurm_array_scripts_root": slurm_artifacts.array_scripts_root,
        "slurm_job_evidence_root": slurm_artifacts.job_evidence_report.evidence_root,
        "slurm_job_evidence_index": slurm_artifacts.job_evidence_report.index_path,
        "slurm_job_evidence_summary": slurm_artifacts.job_evidence_report.summary_path,
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
    }


def _report_manifest_paths(
    workflow_artifacts: RabiesMethodSensitivityWorkflowBundleArtifacts,
    slurm_artifacts: RabiesMethodSensitivitySlurmBundleArtifacts,
    *,
    manifest_path: Path,
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
        "workflow_manifest": manifest_path,
        "slurm_job_plan": slurm_artifacts.job_plan_path,
        "slurm_assumptions": slurm_artifacts.assumptions_path,
        "slurm_summary": slurm_artifacts.summary_path,
        "slurm_array_partitions": slurm_artifacts.array_partitions_path,
        "slurm_array_members": slurm_artifacts.array_members_path,
        "slurm_array_strategy": slurm_artifacts.array_strategy_path,
        "slurm_job_evidence_index": slurm_artifacts.job_evidence_report.index_path,
        "slurm_job_evidence_summary": slurm_artifacts.job_evidence_report.summary_path,
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
    }


def _write_bundle_manifest_artifacts(
    output_root: Path,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    workflow_artifacts: RabiesMethodSensitivityWorkflowBundleArtifacts,
    slurm_artifacts: RabiesMethodSensitivitySlurmBundleArtifacts,
) -> RabiesMethodSensitivityBundleManifestArtifacts:
    manifest_path = _write_manifest(
        output_root / "rabies-method-sensitivity.manifest.json",
        report=report,
        bundle_paths=_workflow_manifest_paths(workflow_artifacts, slurm_artifacts),
    )
    report_manifest_path = _write_report_manifest(
        output_root
        / "report-artifacts"
        / "rabies-method-sensitivity-report.manifest.json",
        report=report,
        bundle_paths=_report_manifest_paths(
            workflow_artifacts,
            slurm_artifacts,
            manifest_path=manifest_path,
        ),
        sha256=_sha256,
    )
    return RabiesMethodSensitivityBundleManifestArtifacts(
        manifest_path=manifest_path,
        report_manifest_path=report_manifest_path,
    )

from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from ..audit import RabiesMethodSensitivityReproducibilityAuditReport
from ..models import RabiesMethodSensitivityPanelWorkflowReport
from ..slurm import (
    RabiesMethodSensitivitySlurmArrayStrategyReport,
    RabiesMethodSensitivitySlurmFailureRecoveryReport,
    RabiesMethodSensitivitySlurmJobEvidenceReport,
    RabiesMethodSensitivitySlurmMergeReport,
    RabiesMethodSensitivitySlurmOutputExplosionReport,
    RabiesMethodSensitivitySlurmOutputFreshnessReport,
    RabiesMethodSensitivitySlurmPlanningReport,
    RabiesMethodSensitivitySlurmStatusReport,
    RabiesMethodSensitivitySlurmStorageReport,
    RabiesMethodSensitivitySlurmTreeRetentionReport,
)
from .sections import (
    _build_report_artifact_links,
    _build_report_embedded_json,
    _build_report_sections,
    _build_report_summary_metrics,
)
from .shared import _relative_bundle_path


def _write_report(
    path: Path,
    *,
    report: RabiesMethodSensitivityPanelWorkflowReport,
    bundle_paths: dict[str, Path],
    report_manifest_path: Path,
    reproducibility_report: RabiesMethodSensitivityReproducibilityAuditReport,
    slurm_planning_report: RabiesMethodSensitivitySlurmPlanningReport,
    slurm_array_strategy_report: RabiesMethodSensitivitySlurmArrayStrategyReport,
    slurm_job_evidence_report: RabiesMethodSensitivitySlurmJobEvidenceReport,
    slurm_storage_report: RabiesMethodSensitivitySlurmStorageReport,
    slurm_output_explosion_report: RabiesMethodSensitivitySlurmOutputExplosionReport,
    slurm_tree_retention_report: RabiesMethodSensitivitySlurmTreeRetentionReport,
    slurm_merge_report: RabiesMethodSensitivitySlurmMergeReport,
    slurm_output_freshness_report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
    slurm_status_report: RabiesMethodSensitivitySlurmStatusReport,
    slurm_failure_recovery_report: RabiesMethodSensitivitySlurmFailureRecoveryReport,
) -> Path:
    sections = _build_report_sections(
        report=report,
        bundle_paths=bundle_paths,
        reproducibility_report=reproducibility_report,
        slurm_planning_report=slurm_planning_report,
        slurm_array_strategy_report=slurm_array_strategy_report,
        slurm_job_evidence_report=slurm_job_evidence_report,
        slurm_storage_report=slurm_storage_report,
        slurm_output_explosion_report=slurm_output_explosion_report,
        slurm_tree_retention_report=slurm_tree_retention_report,
        slurm_merge_report=slurm_merge_report,
        slurm_output_freshness_report=slurm_output_freshness_report,
        slurm_status_report=slurm_status_report,
        slurm_failure_recovery_report=slurm_failure_recovery_report,
    )
    artifact_links = _build_report_artifact_links(path, bundle_paths)
    return write_html_report(
        title="Rabies Method-Sensitivity Report",
        sections=sections,
        out_path=path,
        embedded_json=_build_report_embedded_json(
            path=path,
            report=report,
            report_manifest_path=report_manifest_path,
            reproducibility_report=reproducibility_report,
            slurm_planning_report=slurm_planning_report,
            slurm_array_strategy_report=slurm_array_strategy_report,
            slurm_job_evidence_report=slurm_job_evidence_report,
            slurm_storage_report=slurm_storage_report,
            slurm_output_explosion_report=slurm_output_explosion_report,
            slurm_tree_retention_report=slurm_tree_retention_report,
            slurm_merge_report=slurm_merge_report,
            slurm_output_freshness_report=slurm_output_freshness_report,
            slurm_status_report=slurm_status_report,
            slurm_failure_recovery_report=slurm_failure_recovery_report,
        ),
        summary_metrics=_build_report_summary_metrics(
            report=report,
            report_manifest_path=report_manifest_path,
            reproducibility_report=reproducibility_report,
            slurm_planning_report=slurm_planning_report,
            slurm_array_strategy_report=slurm_array_strategy_report,
            slurm_job_evidence_report=slurm_job_evidence_report,
            slurm_storage_report=slurm_storage_report,
            slurm_output_explosion_report=slurm_output_explosion_report,
            slurm_tree_retention_report=slurm_tree_retention_report,
            slurm_merge_report=slurm_merge_report,
            slurm_output_freshness_report=slurm_output_freshness_report,
            slurm_status_report=slurm_status_report,
            slurm_failure_recovery_report=slurm_failure_recovery_report,
        ),
        artifact_links=[
            *artifact_links,
            (
                "report-manifest",
                _relative_bundle_path(path, report_manifest_path),
                f"{report_manifest_path.stat().st_size} bytes",
            ),
        ],
    )

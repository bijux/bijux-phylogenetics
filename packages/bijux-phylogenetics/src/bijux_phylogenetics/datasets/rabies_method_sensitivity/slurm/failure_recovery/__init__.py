from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    RabiesMethodSensitivitySlurmFailureRecoveryJobRow,
    RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow,
    RabiesMethodSensitivitySlurmFailureRecoveryReport,
)
from .report_builder import build_rabies_method_sensitivity_slurm_failure_recovery_report
from .shared import _write_tsv

__all__ = [
    "RabiesMethodSensitivitySlurmFailureRecoveryJobRow",
    "RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow",
    "RabiesMethodSensitivitySlurmFailureRecoveryReport",
    "build_rabies_method_sensitivity_slurm_failure_recovery_report",
    "write_rabies_method_sensitivity_slurm_failure_recovery_html_report",
    "write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table",
    "write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table",
    "write_rabies_method_sensitivity_slurm_failure_recovery_summary_json",
]


def write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmFailureRecoveryReport,
) -> Path:
    """Write the per-job failure-recovery ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "partition_id",
            "array_index",
            "variant_id",
            "current_status",
            "evidence_class",
            "rerunnable",
            "recovery_action",
            "recovery_scope",
            "likely_cause_code",
            "likely_cause_detail",
            "prerequisite",
            "task_status",
            "error_code",
            "error_message",
            "task_log_path",
            "output_root",
            "detail",
        ),
        rows=[
            {
                "partition_id": row.partition_id,
                "array_index": row.array_index,
                "variant_id": row.variant_id,
                "current_status": row.current_status,
                "evidence_class": row.evidence_class,
                "rerunnable": str(row.rerunnable).lower(),
                "recovery_action": row.recovery_action,
                "recovery_scope": row.recovery_scope,
                "likely_cause_code": row.likely_cause_code,
                "likely_cause_detail": row.likely_cause_detail,
                "prerequisite": row.prerequisite,
                "task_status": row.task_status,
                "error_code": row.error_code,
                "error_message": row.error_message,
                "task_log_path": row.task_log_path,
                "output_root": row.output_root,
                "detail": row.detail,
            }
            for row in report.jobs
        ],
    )


def write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmFailureRecoveryReport,
) -> Path:
    """Write the partition-level failure-recovery ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "partition_id",
            "script_path",
            "overall_status",
            "recovery_action",
            "job_count",
            "rerunnable_job_count",
            "blocked_job_count",
            "failed_job_count",
            "stale_job_count",
            "pending_job_count",
            "variant_ids",
            "likely_cause_codes",
        ),
        rows=[
            {
                "partition_id": row.partition_id,
                "script_path": row.script_path,
                "overall_status": row.overall_status,
                "recovery_action": row.recovery_action,
                "job_count": row.job_count,
                "rerunnable_job_count": row.rerunnable_job_count,
                "blocked_job_count": row.blocked_job_count,
                "failed_job_count": row.failed_job_count,
                "stale_job_count": row.stale_job_count,
                "pending_job_count": row.pending_job_count,
                "variant_ids": ",".join(row.variant_ids),
                "likely_cause_codes": ",".join(row.likely_cause_codes),
            }
            for row in report.partitions
        ],
    )


def write_rabies_method_sensitivity_slurm_failure_recovery_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmFailureRecoveryReport,
) -> Path:
    """Write the machine-readable failure-recovery summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_failure_recovery_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmFailureRecoveryReport,
) -> Path:
    """Write the reviewer-facing failure-recovery report."""
    return write_html_report(
        title="Rabies Slurm Failure Recovery Report",
        sections=[
            (
                "recovery-summary",
                "\n".join(
                    [
                        f"overall_recovery_status: {report.overall_recovery_status}",
                        f"workflow_status: {report.workflow_status}",
                        f"active_run_state: {report.active_run_state}",
                        f"rerunnable_job_count: {report.rerunnable_job_count}",
                        f"blocked_job_count: {report.blocked_job_count}",
                        f"recovery_partition_count: {report.recovery_partition_count}",
                    ]
                ),
            ),
            (
                "rerunnable-jobs",
                "none"
                if report.rerunnable_job_count == 0
                else "\n".join(
                    (
                        f"{row.variant_id}: {row.recovery_action}; "
                        f"{row.likely_cause_code}; {row.likely_cause_detail}"
                    )
                    for row in report.jobs
                    if row.rerunnable
                ),
            ),
            (
                "blocked-jobs",
                "none"
                if report.blocked_job_count == 0
                else "\n".join(
                    f"{row.variant_id}: {row.likely_cause_detail}"
                    for row in report.jobs
                    if row.recovery_action == "wait_for_live_workflow"
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "overall_recovery_status": report.overall_recovery_status,
            "workflow_status": report.workflow_status,
            "active_run_state": report.active_run_state,
            "job_count": report.job_count,
            "rerunnable_job_count": report.rerunnable_job_count,
            "blocked_job_count": report.blocked_job_count,
            "failed_job_count": report.failed_job_count,
            "stale_job_count": report.stale_job_count,
            "pending_job_count": report.pending_job_count,
            "recovery_partition_count": report.recovery_partition_count,
        },
        summary_metrics=[
            ("overall status", report.overall_recovery_status),
            ("rerunnable jobs", report.rerunnable_job_count),
            ("blocked jobs", report.blocked_job_count),
            ("failed jobs", report.failed_job_count),
            ("stale jobs", report.stale_job_count),
            ("recovery partitions", report.recovery_partition_count),
        ],
        artifact_links=[
            ("failure recovery jobs", "slurm-failure-recovery-jobs.tsv", None),
            (
                "failure recovery partitions",
                "slurm-failure-recovery-partitions.tsv",
                None,
            ),
            ("failure recovery summary", "slurm-failure-recovery-report.json", None),
        ],
    )

from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import RabiesMethodSensitivitySlurmFailureRecoveryReport


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

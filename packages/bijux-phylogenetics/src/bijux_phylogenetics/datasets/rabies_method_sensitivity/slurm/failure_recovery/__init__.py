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
from .inputs import load_failure_recovery_inputs
from .recovery_policy import (
    build_failure_recovery_job_row,
    build_failure_recovery_partition_row,
    derive_overall_recovery_status,
)
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

def build_rabies_method_sensitivity_slurm_failure_recovery_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmFailureRecoveryReport:
    """Classify rerunnable jobs and likely failure causes from governed batch evidence."""
    loaded_inputs = load_failure_recovery_inputs(bundle_root)
    bundle_root = loaded_inputs.bundle_root
    config = loaded_inputs.config
    checks = list(loaded_inputs.checks)

    def add_check(
        check_id: str,
        *,
        surface: str,
        condition: bool,
        expected: object,
        observed: object,
        detail: str,
    ) -> None:
        checks.append((check_id, surface, condition, expected, observed, detail))

    configured_variant_ids = loaded_inputs.configured_variant_ids

    job_rows = tuple(
        build_failure_recovery_job_row(bundle_root=bundle_root, job_status_row=row)
        for row in loaded_inputs.job_status_rows
    )
    jobs_by_partition: dict[str, list[RabiesMethodSensitivitySlurmFailureRecoveryJobRow]] = {}
    for row in job_rows:
        jobs_by_partition.setdefault(row.partition_id, []).append(row)

    partition_rows = tuple(
        build_failure_recovery_partition_row(
            partition_status_row=partition_row,
            job_rows=jobs_by_partition.get(str(partition_row["partition_id"]), []),
        )
        for partition_row in loaded_inputs.partition_status_rows
    )
    add_check(
        "partition-status:coverage",
        surface="partition-status",
        condition=len(partition_rows) == len(loaded_inputs.partition_status_rows),
        expected=len(loaded_inputs.partition_status_rows),
        observed=len(partition_rows),
        detail="failure-recovery partition rows cover the written partition-status surface",
    )
    rerunnable_job_count = sum(1 for row in job_rows if row.rerunnable)
    blocked_job_count = sum(
        1 for row in job_rows if row.recovery_action == "wait_for_live_workflow"
    )
    no_action_job_count = sum(1 for row in job_rows if row.recovery_action == "no_action")
    failed_job_count = sum(1 for row in job_rows if row.current_status == "failed")
    stale_job_count = sum(1 for row in job_rows if row.current_status == "stale")
    pending_job_count = sum(1 for row in job_rows if row.current_status == "pending")
    recovery_partition_count = sum(
        1 for row in partition_rows if row.recovery_action != "no_action"
    )

    overall_recovery_status = derive_overall_recovery_status(
        rerunnable_job_count=rerunnable_job_count,
        blocked_job_count=blocked_job_count,
    )
    add_check(
        "recovery-summary:job-count",
        surface="recovery-summary",
        condition=len(job_rows) == len(configured_variant_ids),
        expected=len(configured_variant_ids),
        observed=len(job_rows),
        detail="recovery decisions cover every configured variant exactly once",
    )
    add_check(
        "recovery-summary:overall-status",
        surface="recovery-summary",
        condition=overall_recovery_status
        == derive_overall_recovery_status(
            rerunnable_job_count=rerunnable_job_count,
            blocked_job_count=blocked_job_count,
        ),
        expected=derive_overall_recovery_status(
            rerunnable_job_count=rerunnable_job_count,
            blocked_job_count=blocked_job_count,
        ),
        observed=overall_recovery_status,
        detail="overall recovery status matches the job-level rerun decisions",
    )

    return RabiesMethodSensitivitySlurmFailureRecoveryReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        workflow_status=str(loaded_inputs.workflow_status["workflow_status"]),
        active_run_state=str(loaded_inputs.workflow_status["active_run_state"]),
        overall_recovery_status=overall_recovery_status,
        partition_count=len(partition_rows),
        job_count=len(job_rows),
        rerunnable_job_count=rerunnable_job_count,
        blocked_job_count=blocked_job_count,
        no_action_job_count=no_action_job_count,
        failed_job_count=failed_job_count,
        stale_job_count=stale_job_count,
        pending_job_count=pending_job_count,
        recovery_partition_count=recovery_partition_count,
        check_count=len(checks),
        failed_check_count=sum(1 for row in checks if not row[2]),
        partitions=partition_rows,
        jobs=job_rows,
    )


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

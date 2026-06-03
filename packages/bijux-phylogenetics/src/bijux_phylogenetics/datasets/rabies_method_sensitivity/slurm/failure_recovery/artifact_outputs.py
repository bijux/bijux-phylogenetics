from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmFailureRecoveryReport
from .shared import _write_tsv


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

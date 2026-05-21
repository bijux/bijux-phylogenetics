from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmStatusReport
from .shared import relative_bundle_path, write_tsv


def write_rabies_method_sensitivity_slurm_job_status_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmStatusReport,
) -> Path:
    """Write one per-job resume-state ledger."""
    return write_tsv(
        path,
        fieldnames=(
            "partition_id",
            "array_index",
            "variant_id",
            "status",
            "evidence_class",
            "script_path",
            "task_log_path",
            "output_root",
            "task_status",
            "execution_record_status",
            "output_file_count",
            "output_byte_count",
            "missing_required_file_count",
            "missing_required_files",
            "detail",
            "output_freshness_status",
            "output_freshness_reason_codes",
        ),
        rows=[
            {
                "partition_id": row.partition_id,
                "array_index": row.array_index,
                "variant_id": row.variant_id,
                "status": row.status,
                "evidence_class": row.evidence_class,
                "script_path": row.script_path,
                "task_log_path": row.task_log_path,
                "output_root": row.output_root,
                "task_status": row.task_status,
                "execution_record_status": row.execution_record_status,
                "output_file_count": row.output_file_count,
                "output_byte_count": row.output_byte_count,
                "missing_required_file_count": row.missing_required_file_count,
                "missing_required_files": ",".join(row.missing_required_files),
                "detail": row.detail,
                "output_freshness_status": row.output_freshness_status,
                "output_freshness_reason_codes": ",".join(
                    row.output_freshness_reason_codes
                ),
            }
            for row in report.jobs
        ],
    )


def write_rabies_method_sensitivity_slurm_partition_status_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmStatusReport,
) -> Path:
    """Write one partition-level resume-state ledger."""
    return write_tsv(
        path,
        fieldnames=(
            "partition_id",
            "script_path",
            "job_count",
            "completed_job_count",
            "failed_job_count",
            "pending_job_count",
            "stale_job_count",
            "overall_status",
            "variant_ids",
        ),
        rows=[
            {
                "partition_id": row.partition_id,
                "script_path": row.script_path,
                "job_count": row.job_count,
                "completed_job_count": row.completed_job_count,
                "failed_job_count": row.failed_job_count,
                "pending_job_count": row.pending_job_count,
                "stale_job_count": row.stale_job_count,
                "overall_status": row.overall_status,
                "variant_ids": ",".join(row.variant_ids),
            }
            for row in report.partitions
        ],
    )


def write_rabies_method_sensitivity_slurm_status_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmStatusReport,
) -> Path:
    """Write the structured workflow resume-state summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    payload["execution_record_path"] = relative_bundle_path(
        report.bundle_root,
        report.execution_record_path,
    )
    payload["active_marker_path"] = relative_bundle_path(
        report.bundle_root,
        report.active_marker_path,
    )
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path

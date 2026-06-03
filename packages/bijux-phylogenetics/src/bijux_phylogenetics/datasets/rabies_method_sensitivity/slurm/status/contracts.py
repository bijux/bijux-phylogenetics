from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmJobStatusRow:
    """One real job-state classification for one planned rabies workflow job."""

    partition_id: str
    array_index: int
    variant_id: str
    status: str
    evidence_class: str
    script_path: str
    task_log_path: str
    output_root: str
    task_status: str | None
    execution_record_status: str | None
    output_file_count: int
    output_byte_count: int
    missing_required_file_count: int
    missing_required_files: tuple[str, ...]
    detail: str
    output_freshness_status: str
    output_freshness_reason_codes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmPartitionStatusRow:
    """One aggregate job-state summary for one Slurm array partition."""

    partition_id: str
    script_path: str
    job_count: int
    completed_job_count: int
    failed_job_count: int
    pending_job_count: int
    stale_job_count: int
    overall_status: str
    variant_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmStatusReport:
    """One reviewer-facing status summary for one resumable rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    execution_record_path: Path
    active_marker_path: Path
    active_run_state: str
    workflow_status: str
    partition_count: int
    job_count: int
    completed_job_count: int
    failed_job_count: int
    pending_job_count: int
    stale_job_count: int
    output_freshness_check_count: int
    failed_output_freshness_check_count: int
    fresh_output_job_count: int
    stale_output_job_count: int
    partitions: tuple[RabiesMethodSensitivitySlurmPartitionStatusRow, ...]
    jobs: tuple[RabiesMethodSensitivitySlurmJobStatusRow, ...]

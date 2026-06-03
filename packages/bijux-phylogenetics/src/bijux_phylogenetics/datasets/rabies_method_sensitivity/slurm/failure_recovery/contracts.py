from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmFailureRecoveryJobRow:
    """One recovery decision for one planned rabies Slurm job."""

    partition_id: str
    array_index: int
    variant_id: str
    current_status: str
    evidence_class: str
    rerunnable: bool
    recovery_action: str
    recovery_scope: str
    likely_cause_code: str
    likely_cause_detail: str
    prerequisite: str
    task_status: str | None
    error_code: str | None
    error_message: str | None
    task_log_path: str
    output_root: str
    detail: str


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow:
    """One partition-level rerun recommendation over planned rabies jobs."""

    partition_id: str
    script_path: str
    overall_status: str
    recovery_action: str
    job_count: int
    rerunnable_job_count: int
    blocked_job_count: int
    failed_job_count: int
    stale_job_count: int
    pending_job_count: int
    variant_ids: tuple[str, ...]
    likely_cause_codes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmFailureRecoveryReport:
    """One reviewer-facing failure-recovery summary over the governed rabies batch workflow."""

    dataset_id: str
    workflow_prefix: str
    bundle_root: Path
    workflow_status: str
    active_run_state: str
    overall_recovery_status: str
    partition_count: int
    job_count: int
    rerunnable_job_count: int
    blocked_job_count: int
    no_action_job_count: int
    failed_job_count: int
    stale_job_count: int
    pending_job_count: int
    recovery_partition_count: int
    check_count: int
    failed_check_count: int
    partitions: tuple[RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow, ...]
    jobs: tuple[RabiesMethodSensitivitySlurmFailureRecoveryJobRow, ...]

from __future__ import annotations

from pathlib import Path

from .contracts import (
    RabiesMethodSensitivitySlurmFailureRecoveryJobRow,
    RabiesMethodSensitivitySlurmFailureRecoveryReport,
)
from .inputs import load_failure_recovery_inputs
from .recovery_policy import (
    build_failure_recovery_job_row,
    build_failure_recovery_partition_row,
    derive_overall_recovery_status,
)


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
    jobs_by_partition: dict[
        str, list[RabiesMethodSensitivitySlurmFailureRecoveryJobRow]
    ] = {}
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
    no_action_job_count = sum(
        1 for row in job_rows if row.recovery_action == "no_action"
    )
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

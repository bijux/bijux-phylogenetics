from __future__ import annotations

from pathlib import Path

from .contracts import (
    RabiesMethodSensitivitySlurmFailureRecoveryJobRow,
    RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow,
)
from .shared import (
    _TERMINAL_FAILURE_CODES,
    _normalize_optional,
    _parse_task_log,
)


def build_failure_recovery_job_row(
    *,
    bundle_root: Path,
    job_status_row: dict[str, str],
) -> RabiesMethodSensitivitySlurmFailureRecoveryJobRow:
    task_log_path = bundle_root / str(job_status_row["task_log_path"])
    task_log = _parse_task_log(task_log_path) if task_log_path.is_file() else {}
    current_status = str(job_status_row["status"])
    evidence_class = str(job_status_row["evidence_class"])
    error_code = _normalize_optional(task_log.get("error_code"))
    error_message = _normalize_optional(task_log.get("error_message"))

    rerunnable = False
    recovery_action = "no_action"
    recovery_scope = "job"
    likely_cause_code = "none"
    likely_cause_detail = "job completed cleanly and does not need recovery"
    prerequisite = "none"

    if current_status == "completed":
        pass
    elif current_status == "pending":
        rerunnable = False
        recovery_action = "wait_for_live_workflow"
        recovery_scope = "workflow"
        likely_cause_code = "workflow_still_running"
        likely_cause_detail = (
            "the workflow is still active and this job should not be rerun yet"
        )
        prerequisite = "wait_for_terminal_state"
    elif current_status == "failed":
        rerunnable = True
        recovery_action = "rerun_variant"
        likely_cause_code, likely_cause_detail = classify_terminal_failure(
            error_code=error_code,
            error_message=error_message,
        )
        prerequisite = (
            "restore_inputs_then_rerun"
            if likely_cause_code == "missing_input"
            else "inspect_task_log_then_rerun"
        )
    elif current_status == "stale":
        rerunnable = True
        likely_cause_code, likely_cause_detail, recovery_action, prerequisite = (
            classify_stale_recovery(job_status_row)
        )
    else:
        rerunnable = True
        recovery_action = "rerun_variant"
        likely_cause_code = "unknown_state"
        likely_cause_detail = (
            "job state could not be classified cleanly from the governed evidence"
        )
        prerequisite = "inspect_status_ledgers_then_rerun"

    return RabiesMethodSensitivitySlurmFailureRecoveryJobRow(
        partition_id=str(job_status_row["partition_id"]),
        array_index=int(job_status_row["array_index"]),
        variant_id=str(job_status_row["variant_id"]),
        current_status=current_status,
        evidence_class=evidence_class,
        rerunnable=rerunnable,
        recovery_action=recovery_action,
        recovery_scope=recovery_scope,
        likely_cause_code=likely_cause_code,
        likely_cause_detail=likely_cause_detail,
        prerequisite=prerequisite,
        task_status=_normalize_optional(job_status_row.get("task_status")),
        error_code=error_code,
        error_message=error_message,
        task_log_path=str(job_status_row["task_log_path"]),
        output_root=str(job_status_row["output_root"]),
        detail=str(job_status_row["detail"]),
    )


def build_failure_recovery_partition_row(
    *,
    partition_status_row: dict[str, str],
    job_rows: list[RabiesMethodSensitivitySlurmFailureRecoveryJobRow],
) -> RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow:
    rerunnable_job_count = sum(1 for row in job_rows if row.rerunnable)
    blocked_job_count = sum(
        1 for row in job_rows if row.recovery_action == "wait_for_live_workflow"
    )
    failed_job_count = sum(1 for row in job_rows if row.current_status == "failed")
    stale_job_count = sum(1 for row in job_rows if row.current_status == "stale")
    pending_job_count = sum(1 for row in job_rows if row.current_status == "pending")
    recovery_action = "no_action"
    if blocked_job_count > 0 and rerunnable_job_count == 0:
        recovery_action = "wait_for_live_workflow"
    elif rerunnable_job_count == len(job_rows) and len(job_rows) > 0:
        recovery_action = "rerun_partition"
    elif rerunnable_job_count > 0:
        recovery_action = "rerun_selected_jobs"
    likely_cause_codes = tuple(
        sorted(
            {
                row.likely_cause_code
                for row in job_rows
                if row.likely_cause_code != "none"
            }
        )
    )
    return RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow(
        partition_id=str(partition_status_row["partition_id"]),
        script_path=str(partition_status_row["script_path"]),
        overall_status=str(partition_status_row["overall_status"]),
        recovery_action=recovery_action,
        job_count=len(job_rows),
        rerunnable_job_count=rerunnable_job_count,
        blocked_job_count=blocked_job_count,
        failed_job_count=failed_job_count,
        stale_job_count=stale_job_count,
        pending_job_count=pending_job_count,
        variant_ids=tuple(row.variant_id for row in job_rows),
        likely_cause_codes=likely_cause_codes,
    )


def derive_overall_recovery_status(
    *,
    rerunnable_job_count: int,
    blocked_job_count: int,
) -> str:
    if rerunnable_job_count > 0:
        return "recovery_needed"
    if blocked_job_count > 0:
        return "workflow_active"
    return "clean"


def classify_terminal_failure(
    *,
    error_code: str | None,
    error_message: str | None,
) -> tuple[str, str]:
    if error_code in _TERMINAL_FAILURE_CODES:
        code = _TERMINAL_FAILURE_CODES[error_code]
        if code == "task_timeout":
            return (code, "engine execution timed out before the variant finished")
        if code == "incomplete_outputs":
            return (code, "engine reported success but durable outputs were missing")
        if code == "stale_running_marker":
            return (code, "a stale running marker blocked a clean variant completion")
        return (
            code,
            "variant execution failed and the task log recorded a terminal error",
        )
    if error_message is not None and "No such file" in error_message:
        return ("missing_input", "required input or intermediate files were missing")
    if error_message is not None and "timeout" in error_message.lower():
        return (
            "task_timeout",
            "engine execution timed out before the variant finished",
        )
    if error_message is not None and "killed" in error_message.lower():
        return ("task_terminated", "engine execution was terminated before completion")
    return ("task_failure", "task log recorded a terminal variant failure")


def classify_stale_recovery(
    job_status_row: dict[str, str],
) -> tuple[str, str, str, str]:
    evidence_class = str(job_status_row["evidence_class"])
    freshness_status = str(job_status_row["output_freshness_status"])
    freshness_codes = {
        code
        for code in str(job_status_row["output_freshness_reason_codes"]).split(",")
        if code
    }
    if freshness_status == "stale":
        return (
            "output_drift",
            "completed outputs no longer match current packaged inputs or workflow settings",
            "refresh_inputs_then_rerun",
            "confirm_current_inputs_then_rerun",
        )
    if evidence_class == "stale-running-marker":
        return (
            "stale_running_marker",
            "a dead running marker remained after the workflow stopped",
            "clear_stale_marker_then_rerun",
            "clear_stale_marker",
        )
    if evidence_class == "incomplete-success-output":
        return (
            "incomplete_outputs",
            "success was recorded but one or more required durable outputs are missing",
            "rerun_variant",
            "inspect_missing_outputs_then_rerun",
        )
    if evidence_class == "failed-workflow-gap":
        return (
            "workflow_gap",
            "the workflow failed before this planned job reached a terminal record",
            "rerun_variant",
            "inspect_parallel_failure_then_rerun",
        )
    if evidence_class == "abandoned-partial-output":
        detail = "the job wrote partial evidence without a matching terminal record"
        if "missing_required_outputs" in freshness_codes:
            detail = "partial outputs are missing required files and should be rerun from scratch"
        return (
            "abandoned_partial_output",
            detail,
            "rerun_variant",
            "clean_partial_outputs_then_rerun",
        )
    return (
        "stale_job_state",
        "job evidence is stale and needs a fresh rerun to recover a coherent state",
        "rerun_variant",
        "inspect_status_ledgers_then_rerun",
    )

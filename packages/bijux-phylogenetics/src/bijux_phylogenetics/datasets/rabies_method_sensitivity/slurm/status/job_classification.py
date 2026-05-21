from __future__ import annotations

from pathlib import Path

from ..freshness import RabiesMethodSensitivitySlurmOutputFreshnessRow
from .contracts import RabiesMethodSensitivitySlurmJobStatusRow
from .shared import (
    collect_output_observations,
    missing_required_variant_outputs,
    parse_task_log,
)


def build_job_status_row(
    *,
    bundle_root: Path,
    member_row: dict[str, str],
    execution_task_row: dict[str, object] | None,
    workflow_status: str | None,
    active_run_state: str,
    freshness_row: RabiesMethodSensitivitySlurmOutputFreshnessRow | None,
) -> RabiesMethodSensitivitySlurmJobStatusRow:
    variant_id = str(member_row["variant_id"])
    task_log_path = bundle_root / str(member_row["task_log_path"])
    output_root = bundle_root / str(member_row["bundle_output_directory"])
    task_log = parse_task_log(task_log_path) if task_log_path.is_file() else {}
    task_status = resolve_task_status(
        execution_task_row=execution_task_row,
        task_log=task_log,
    )
    output_file_count, output_byte_count = collect_output_observations(output_root)
    missing_required_files = missing_required_variant_outputs(output_root, variant_id)
    has_output_evidence = output_root.exists()
    has_terminal_failure = task_status == "failed"
    has_terminal_success = task_status == "succeeded"
    has_stale_marker = active_run_state == "stale"

    if has_terminal_failure:
        status = "failed"
        evidence_class = "terminal-failure"
        detail = "task log or execution record marks the job as failed"
    elif has_terminal_success and not missing_required_files:
        status = "completed"
        evidence_class = "terminal-success"
        detail = "task log and durable outputs confirm a complete job result"
    elif has_terminal_success and missing_required_files:
        status = "stale"
        evidence_class = "incomplete-success-output"
        detail = (
            "success was recorded but one or more required durable outputs are missing"
        )
    elif active_run_state == "live":
        status = "pending"
        evidence_class = (
            "live-partial-output" if has_output_evidence else "planned-not-started"
        )
        detail = (
            "workflow is still running and this job has partial execution evidence"
            if has_output_evidence
            else "workflow is still running and this planned job has not written durable evidence yet"
        )
    elif has_stale_marker or has_output_evidence or task_log_path.is_file():
        status = "stale"
        evidence_class = (
            "stale-running-marker" if has_stale_marker else "abandoned-partial-output"
        )
        detail = (
            "a dead running marker remains without a live workflow owner"
            if has_stale_marker
            else "the job wrote partial evidence without a matching terminal success or failure record"
        )
    else:
        status = "pending"
        evidence_class = "planned-not-started"
        detail = "the job is planned but has not started writing durable evidence"

    if task_status is None and workflow_status == "failed" and status == "pending":
        status = "stale"
        evidence_class = "failed-workflow-gap"
        detail = (
            "the workflow ended in failure before this job reached a terminal record"
        )
    if (
        freshness_row is not None
        and freshness_row.freshness_status == "stale"
        and status == "completed"
    ):
        status = "stale"
        evidence_class = "stale-output-drift"
        detail = freshness_row.stale_reason_detail

    return RabiesMethodSensitivitySlurmJobStatusRow(
        partition_id=str(member_row["partition_id"]),
        array_index=int(member_row["array_index"]),
        variant_id=variant_id,
        status=status,
        evidence_class=evidence_class,
        script_path=str(member_row["script_path"]),
        task_log_path=str(member_row["task_log_path"]),
        output_root=str(member_row["bundle_output_directory"]),
        task_status=task_status,
        execution_record_status=(
            None
            if execution_task_row is None or execution_task_row.get("status") is None
            else str(execution_task_row["status"])
        ),
        output_file_count=output_file_count,
        output_byte_count=output_byte_count,
        missing_required_file_count=len(missing_required_files),
        missing_required_files=missing_required_files,
        detail=detail,
        output_freshness_status=(
            "unknown" if freshness_row is None else freshness_row.freshness_status
        ),
        output_freshness_reason_codes=(
            () if freshness_row is None else freshness_row.stale_reason_codes
        ),
    )


def resolve_task_status(
    *,
    execution_task_row: dict[str, object] | None,
    task_log: dict[str, str],
) -> str | None:
    task_status = task_log.get("status")
    if task_status is not None:
        return task_status
    if execution_task_row is None or execution_task_row.get("status") is None:
        return None
    return str(execution_task_row["status"])

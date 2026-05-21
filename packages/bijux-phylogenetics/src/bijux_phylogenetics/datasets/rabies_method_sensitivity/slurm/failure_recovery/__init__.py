from __future__ import annotations

from dataclasses import asdict
import csv
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    RabiesMethodSensitivitySlurmFailureRecoveryJobRow,
    RabiesMethodSensitivitySlurmFailureRecoveryPartitionRow,
    RabiesMethodSensitivitySlurmFailureRecoveryReport,
)

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

_CONFIG_FILENAME = "workflow-config.resolved.json"
_SLURM_JOB_STATUS_FILENAME = "slurm-job-status.tsv"
_SLURM_PARTITION_STATUS_FILENAME = "slurm-partition-status.tsv"
_SLURM_WORKFLOW_STATUS_FILENAME = "slurm-workflow-status.json"
_TERMINAL_FAILURE_CODES = {
    "parallel_variant_failed": "task_failure",
    "engine_process_timeout": "task_timeout",
    "engine_output_missing": "incomplete_outputs",
    "engine_workflow_running_marker_invalid": "stale_running_marker",
}


def build_rabies_method_sensitivity_slurm_failure_recovery_report(
    bundle_root: Path,
) -> RabiesMethodSensitivitySlurmFailureRecoveryReport:
    """Classify rerunnable jobs and likely failure causes from governed batch evidence."""
    bundle_root = bundle_root.resolve()
    config = _load_json(bundle_root / _CONFIG_FILENAME)
    job_status_rows = _read_tsv_rows(bundle_root / _SLURM_JOB_STATUS_FILENAME)
    partition_status_rows = _read_tsv_rows(bundle_root / _SLURM_PARTITION_STATUS_FILENAME)
    workflow_status = _load_json(bundle_root / _SLURM_WORKFLOW_STATUS_FILENAME)

    checks: list[tuple[str, str, bool, object, object, str]] = []

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

    configured_variant_ids = sorted(
        str(row["variant_id"]) for row in list(config.get("variants", []))
    )
    observed_variant_ids = sorted(str(row["variant_id"]) for row in job_status_rows)
    add_check(
        "job-status:variant-coverage",
        surface="job-status",
        condition=observed_variant_ids == configured_variant_ids,
        expected=configured_variant_ids,
        observed=observed_variant_ids,
        detail="job-status rows cover the configured variant ids",
    )

    job_rows = tuple(
        _build_job_row(bundle_root=bundle_root, job_status_row=row)
        for row in job_status_rows
    )
    jobs_by_partition: dict[str, list[RabiesMethodSensitivitySlurmFailureRecoveryJobRow]] = {}
    for row in job_rows:
        jobs_by_partition.setdefault(row.partition_id, []).append(row)

    partition_rows = tuple(
        _build_partition_row(
            partition_status_row=partition_row,
            job_rows=jobs_by_partition.get(str(partition_row["partition_id"]), []),
        )
        for partition_row in partition_status_rows
    )
    add_check(
        "partition-status:coverage",
        surface="partition-status",
        condition=len(partition_rows) == len(partition_status_rows),
        expected=len(partition_status_rows),
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

    overall_recovery_status = "clean"
    if blocked_job_count > 0 and rerunnable_job_count == 0:
        overall_recovery_status = "workflow_active"
    elif rerunnable_job_count > 0:
        overall_recovery_status = "recovery_needed"
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
        == _derive_overall_recovery_status(
            rerunnable_job_count=rerunnable_job_count,
            blocked_job_count=blocked_job_count,
        ),
        expected=_derive_overall_recovery_status(
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
        workflow_status=str(workflow_status["workflow_status"]),
        active_run_state=str(workflow_status["active_run_state"]),
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


def _build_job_row(
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
        likely_cause_detail = "the workflow is still active and this job should not be rerun yet"
        prerequisite = "wait_for_terminal_state"
    elif current_status == "failed":
        rerunnable = True
        recovery_action = "rerun_variant"
        likely_cause_code, likely_cause_detail = _classify_terminal_failure(
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
            _classify_stale_recovery(job_status_row)
        )
    else:
        rerunnable = True
        recovery_action = "rerun_variant"
        likely_cause_code = "unknown_state"
        likely_cause_detail = "job state could not be classified cleanly from the governed evidence"
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


def _build_partition_row(
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
        sorted({row.likely_cause_code for row in job_rows if row.likely_cause_code != "none"})
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


def _derive_overall_recovery_status(
    *,
    rerunnable_job_count: int,
    blocked_job_count: int,
) -> str:
    if rerunnable_job_count > 0:
        return "recovery_needed"
    if blocked_job_count > 0:
        return "workflow_active"
    return "clean"


def _classify_terminal_failure(
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
        return (code, "variant execution failed and the task log recorded a terminal error")
    if error_message is not None and "No such file" in error_message:
        return ("missing_input", "required input or intermediate files were missing")
    if error_message is not None and "timeout" in error_message.lower():
        return ("task_timeout", "engine execution timed out before the variant finished")
    if error_message is not None and "killed" in error_message.lower():
        return ("task_terminated", "engine execution was terminated before completion")
    return ("task_failure", "task log recorded a terminal variant failure")


def _classify_stale_recovery(
    job_status_row: dict[str, str],
) -> tuple[str, str, str, str]:
    evidence_class = str(job_status_row["evidence_class"])
    freshness_status = str(job_status_row["output_freshness_status"])
    freshness_codes = {
        code for code in str(job_status_row["output_freshness_reason_codes"]).split(",") if code
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


def _parse_task_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        if ": " not in raw_line:
            continue
        key, value = raw_line.split(": ", 1)
        values[key.strip()] = value.strip()
    return values


def _normalize_optional(value: object | None) -> str | None:
    if value is None:
        return None
    text = str(value)
    return None if text == "" else text


def _load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return path

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.engines.common import (
    active_engine_run_is_live,
    engine_active_marker_path,
    load_active_engine_run,
)

from ..freshness import (
    RabiesMethodSensitivitySlurmOutputFreshnessRow,
    build_rabies_method_sensitivity_slurm_output_freshness_report,
)
from .contracts import (
    RabiesMethodSensitivitySlurmJobStatusRow,
    RabiesMethodSensitivitySlurmPartitionStatusRow,
    RabiesMethodSensitivitySlurmStatusReport,
)
from .shared import (
    _CONFIG_FILENAME,
    _EXECUTION_RECORD_FILENAME,
    _SLURM_ARRAY_MEMBERS_FILENAME,
    _SLURM_ARRAY_PARTITIONS_FILENAME,
    collect_output_observations,
    load_json,
    missing_required_variant_outputs,
    parse_task_log,
    read_tsv_rows,
    relative_bundle_path,
    write_tsv,
)

__all__ = [
    "RabiesMethodSensitivitySlurmJobStatusRow",
    "RabiesMethodSensitivitySlurmPartitionStatusRow",
    "RabiesMethodSensitivitySlurmStatusReport",
    "build_rabies_method_sensitivity_slurm_status_report",
    "write_rabies_method_sensitivity_slurm_job_status_table",
    "write_rabies_method_sensitivity_slurm_partition_status_table",
    "write_rabies_method_sensitivity_slurm_status_json",
]

def build_rabies_method_sensitivity_slurm_status_report(
    bundle_root: Path,
    *,
    dataset: object | None = None,
) -> RabiesMethodSensitivitySlurmStatusReport:
    """Classify planned rabies workflow jobs by real completion state."""
    bundle_root = bundle_root.resolve()
    config = load_json(bundle_root / _CONFIG_FILENAME)
    execution_record_path = bundle_root / _EXECUTION_RECORD_FILENAME
    execution_record = (
        load_json(execution_record_path) if execution_record_path.is_file() else {}
    )
    active_marker_path = engine_active_marker_path(execution_record_path)
    active_record = load_active_engine_run(execution_record_path)
    active_run_state = "absent"
    if active_record is not None:
        active_run_state = "live" if active_engine_run_is_live(active_record) else "stale"

    partition_rows = read_tsv_rows(bundle_root / _SLURM_ARRAY_PARTITIONS_FILENAME)
    member_rows = read_tsv_rows(bundle_root / _SLURM_ARRAY_MEMBERS_FILENAME)
    freshness_report = build_rabies_method_sensitivity_slurm_output_freshness_report(
        bundle_root,
        dataset=dataset,
    )
    freshness_rows_by_variant = {
        row.variant_id: row for row in freshness_report.jobs
    }
    execution_task_rows = {
        str(row["variant_id"]): row
        for row in list(execution_record.get("task_records", []))
    }
    job_rows = tuple(
        _build_job_status_row(
            bundle_root=bundle_root,
            member_row=member_row,
            execution_task_row=execution_task_rows.get(str(member_row["variant_id"])),
            workflow_status=(
                str(execution_record.get("status"))
                if execution_record.get("status") is not None
                else None
            ),
            active_run_state=active_run_state,
            freshness_row=freshness_rows_by_variant.get(str(member_row["variant_id"])),
        )
        for member_row in member_rows
    )
    jobs_by_partition: dict[str, list[RabiesMethodSensitivitySlurmJobStatusRow]] = {}
    for row in job_rows:
        jobs_by_partition.setdefault(row.partition_id, []).append(row)
    partition_status_rows = tuple(
        _build_partition_status_row(
            partition_row=partition_row,
            job_rows=jobs_by_partition.get(str(partition_row["partition_id"]), []),
        )
        for partition_row in partition_rows
    )
    return RabiesMethodSensitivitySlurmStatusReport(
        dataset_id=str(config["dataset_id"]),
        workflow_prefix=str(config["workflow_prefix"]),
        bundle_root=bundle_root,
        execution_record_path=execution_record_path,
        active_marker_path=active_marker_path,
        active_run_state=active_run_state,
        workflow_status=str(execution_record.get("status", "unknown")),
        partition_count=len(partition_status_rows),
        job_count=len(job_rows),
        completed_job_count=sum(1 for row in job_rows if row.status == "completed"),
        failed_job_count=sum(1 for row in job_rows if row.status == "failed"),
        pending_job_count=sum(1 for row in job_rows if row.status == "pending"),
        stale_job_count=sum(1 for row in job_rows if row.status == "stale"),
        output_freshness_check_count=freshness_report.check_count,
        failed_output_freshness_check_count=freshness_report.failed_check_count,
        fresh_output_job_count=freshness_report.fresh_job_count,
        stale_output_job_count=freshness_report.stale_job_count,
        partitions=partition_status_rows,
        jobs=job_rows,
    )


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


def _build_job_status_row(
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
    task_status = _resolve_task_status(
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
        detail = "success was recorded but one or more required durable outputs are missing"
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
            "stale-running-marker"
            if has_stale_marker
            else "abandoned-partial-output"
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
        detail = "the workflow ended in failure before this job reached a terminal record"
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


def _build_partition_status_row(
    *,
    partition_row: dict[str, str],
    job_rows: list[RabiesMethodSensitivitySlurmJobStatusRow],
) -> RabiesMethodSensitivitySlurmPartitionStatusRow:
    completed_job_count = sum(1 for row in job_rows if row.status == "completed")
    failed_job_count = sum(1 for row in job_rows if row.status == "failed")
    pending_job_count = sum(1 for row in job_rows if row.status == "pending")
    stale_job_count = sum(1 for row in job_rows if row.status == "stale")
    overall_status = _partition_overall_status(
        job_count=len(job_rows),
        completed_job_count=completed_job_count,
        failed_job_count=failed_job_count,
        pending_job_count=pending_job_count,
        stale_job_count=stale_job_count,
    )
    return RabiesMethodSensitivitySlurmPartitionStatusRow(
        partition_id=str(partition_row["partition_id"]),
        script_path=str(partition_row["script_path"]),
        job_count=len(job_rows),
        completed_job_count=completed_job_count,
        failed_job_count=failed_job_count,
        pending_job_count=pending_job_count,
        stale_job_count=stale_job_count,
        overall_status=overall_status,
        variant_ids=tuple(row.variant_id for row in sorted(job_rows, key=_job_sort_key)),
    )


def _partition_overall_status(
    *,
    job_count: int,
    completed_job_count: int,
    failed_job_count: int,
    pending_job_count: int,
    stale_job_count: int,
) -> str:
    if completed_job_count == job_count:
        return "completed"
    non_zero_states = [
        name
        for name, count in (
            ("failed", failed_job_count),
            ("pending", pending_job_count),
            ("stale", stale_job_count),
        )
        if count > 0
    ]
    if len(non_zero_states) == 1 and completed_job_count == 0:
        return non_zero_states[0]
    return "mixed"


def _job_sort_key(row: RabiesMethodSensitivitySlurmJobStatusRow) -> tuple[int, str]:
    return row.array_index, row.variant_id


def _resolve_task_status(
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


from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.common import (
    active_engine_run_is_live,
    engine_active_marker_path,
    load_active_engine_run,
)

from ..freshness import build_rabies_method_sensitivity_slurm_output_freshness_report
from .contracts import (
    RabiesMethodSensitivitySlurmJobStatusRow,
    RabiesMethodSensitivitySlurmPartitionStatusRow,
    RabiesMethodSensitivitySlurmStatusReport,
)
from .job_classification import build_job_status_row
from .shared import (
    _CONFIG_FILENAME,
    _EXECUTION_RECORD_FILENAME,
    _SLURM_ARRAY_MEMBERS_FILENAME,
    _SLURM_ARRAY_PARTITIONS_FILENAME,
    load_json,
    read_tsv_rows,
)


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
        active_run_state = (
            "live" if active_engine_run_is_live(active_record) else "stale"
        )

    partition_rows = read_tsv_rows(bundle_root / _SLURM_ARRAY_PARTITIONS_FILENAME)
    member_rows = read_tsv_rows(bundle_root / _SLURM_ARRAY_MEMBERS_FILENAME)
    freshness_report = build_rabies_method_sensitivity_slurm_output_freshness_report(
        bundle_root,
        dataset=dataset,
    )
    freshness_rows_by_variant = {row.variant_id: row for row in freshness_report.jobs}
    execution_task_rows = {
        str(row["variant_id"]): row
        for row in list(execution_record.get("task_records", []))
    }
    job_rows = tuple(
        build_job_status_row(
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
        build_partition_status_row(
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


def build_partition_status_row(
    *,
    partition_row: dict[str, str],
    job_rows: list[RabiesMethodSensitivitySlurmJobStatusRow],
) -> RabiesMethodSensitivitySlurmPartitionStatusRow:
    completed_job_count = sum(1 for row in job_rows if row.status == "completed")
    failed_job_count = sum(1 for row in job_rows if row.status == "failed")
    pending_job_count = sum(1 for row in job_rows if row.status == "pending")
    stale_job_count = sum(1 for row in job_rows if row.status == "stale")
    overall_status = partition_overall_status(
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
        variant_ids=tuple(row.variant_id for row in sorted(job_rows, key=job_sort_key)),
    )


def partition_overall_status(
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


def job_sort_key(row: RabiesMethodSensitivitySlurmJobStatusRow) -> tuple[int, str]:
    return row.array_index, row.variant_id

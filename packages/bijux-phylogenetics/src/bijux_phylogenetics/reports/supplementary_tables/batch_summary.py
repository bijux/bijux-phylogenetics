from __future__ import annotations

import csv
import json
from pathlib import Path

from .columns import batch_summary_table_columns
from .models import (
    SupplementaryBatchSummaryRow,
    SupplementaryBatchSummaryTableResult,
)
from .shared import stringify_list, table_delimiter, write_dict_rows


def _split_batch_values(value: str) -> list[str]:
    stripped = value.strip()
    if not stripped:
        return []
    for separator in ("|", "; ", ";"):
        if separator in stripped:
            return [item for item in stripped.split(separator) if item]
    return [stripped]


def _read_bundle_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=table_delimiter(path)))


def _read_bundle_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def _optional_int(value: str | None) -> int | None:
    if value is None:
        return None
    stripped = value.strip()
    if not stripped:
        return None
    return int(stripped)


def _serialize_batch_summary_row(
    row: SupplementaryBatchSummaryRow,
) -> dict[str, str]:
    return {
        "row_scope": row.row_scope,
        "dataset_id": row.dataset_id,
        "dataset_label": row.dataset_label,
        "workflow_status": row.workflow_status,
        "variant_id": "" if row.variant_id is None else row.variant_id,
        "label": "" if row.label is None else row.label,
        "execution_mode": "" if row.execution_mode is None else row.execution_mode,
        "task_status": "" if row.task_status is None else row.task_status,
        "job_status": "" if row.job_status is None else row.job_status,
        "output_freshness_status": (
            "" if row.output_freshness_status is None else row.output_freshness_status
        ),
        "recovery_action": "" if row.recovery_action is None else row.recovery_action,
        "merge_status": "" if row.merge_status is None else row.merge_status,
        "evidence_status": "" if row.evidence_status is None else row.evidence_status,
        "reproducibility_status": (
            "" if row.reproducibility_status is None else row.reproducibility_status
        ),
        "selected_model": "" if row.selected_model is None else row.selected_model,
        "output_root": "" if row.output_root is None else row.output_root,
        "task_log_path": "" if row.task_log_path is None else row.task_log_path,
        "evidence_json_path": (
            "" if row.evidence_json_path is None else row.evidence_json_path
        ),
        "evidence_html_path": (
            "" if row.evidence_html_path is None else row.evidence_html_path
        ),
        "variant_count": "" if row.variant_count is None else str(row.variant_count),
        "successful_variant_count": (
            ""
            if row.successful_variant_count is None
            else str(row.successful_variant_count)
        ),
        "failed_variant_count": (
            "" if row.failed_variant_count is None else str(row.failed_variant_count)
        ),
        "output_file_count": (
            "" if row.output_file_count is None else str(row.output_file_count)
        ),
        "output_byte_count": (
            "" if row.output_byte_count is None else str(row.output_byte_count)
        ),
        "artifact_file_count": (
            "" if row.artifact_file_count is None else str(row.artifact_file_count)
        ),
        "linked_artifact_count": (
            "" if row.linked_artifact_count is None else str(row.linked_artifact_count)
        ),
        "linked_artifact_bytes": (
            "" if row.linked_artifact_bytes is None else str(row.linked_artifact_bytes)
        ),
        "issue_count": str(row.issue_count),
        "issues": stringify_list(row.issues),
        "error_code": "" if row.error_code is None else row.error_code,
        "error_message": "" if row.error_message is None else row.error_message,
        "job_evidence_warning_count": (
            ""
            if row.job_evidence_warning_count is None
            else str(row.job_evidence_warning_count)
        ),
        "warning_count": str(row.warning_count),
        "warnings": stringify_list(row.warnings),
    }


def _write_batch_summary_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryBatchSummaryRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_batch_summary_row(row) for row in rows],
    )


def _bundle_output_totals(bundle_root: Path) -> tuple[int, int]:
    file_paths = [path for path in bundle_root.rglob("*") if path.is_file()]
    return len(file_paths), sum(path.stat().st_size for path in file_paths)


def _maybe_path(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return stripped if stripped else None


def _build_batch_variant_row(
    *,
    dataset_id: str,
    dataset_label: str,
    workflow_status: str,
    task_row: dict[str, str],
    task_record: dict[str, object],
    variant_summary_row: dict[str, str],
    job_status_row: dict[str, str],
    freshness_row: dict[str, str],
    recovery_row: dict[str, str],
    merge_row: dict[str, str],
    evidence_row: dict[str, str],
    reproducibility_row: dict[str, str],
) -> SupplementaryBatchSummaryRow:
    issues = (
        _split_batch_values(merge_row.get("issues", ""))
        + _split_batch_values(reproducibility_row.get("issues", ""))
        + _split_batch_values(freshness_row.get("stale_reason_codes", ""))
        + _split_batch_values(job_status_row.get("output_freshness_reason_codes", ""))
    )
    warnings: list[str] = []
    task_status = _maybe_path(task_row.get("status"))
    job_status = _maybe_path(job_status_row.get("status"))
    freshness_status = _maybe_path(freshness_row.get("freshness_status"))
    recovery_action = _maybe_path(recovery_row.get("recovery_action"))
    merge_status = _maybe_path(merge_row.get("merge_status"))
    evidence_status = _maybe_path(merge_row.get("evidence_status"))
    reproducibility_status = _maybe_path(reproducibility_row.get("status"))
    error_code = _maybe_path(task_row.get("error_code")) or _maybe_path(
        task_record.get("error_code") if isinstance(task_record, dict) else None
    )
    error_message = _maybe_path(
        task_record.get("error_message") if isinstance(task_record, dict) else None
    )
    if task_status is not None and task_status != "succeeded":
        warnings.append(f"task-status:{task_status}")
    if error_code is not None:
        warnings.append(f"task-error:{error_code}")
    if job_status is not None and job_status != "completed":
        warnings.append(f"job-status:{job_status}")
    if freshness_status is not None and freshness_status != "fresh":
        warnings.append(f"freshness-status:{freshness_status}")
    if recovery_action is not None and recovery_action != "no_action":
        warnings.append(f"recovery-action:{recovery_action}")
    if merge_status is not None and merge_status != "merged":
        warnings.append(f"merge-status:{merge_status}")
    if evidence_status is not None and evidence_status != "present":
        warnings.append(f"evidence-status:{evidence_status}")
    if reproducibility_status is not None and reproducibility_status != "passed":
        warnings.append(f"reproducibility-status:{reproducibility_status}")
    job_evidence_warning_count = _optional_int(evidence_row.get("warning_count"))
    if job_evidence_warning_count is not None and job_evidence_warning_count > 0:
        warnings.append(f"job-evidence-warning-count:{job_evidence_warning_count}")
    warnings.extend(issue for issue in issues if issue not in warnings)
    return SupplementaryBatchSummaryRow(
        row_scope="variant",
        dataset_id=dataset_id,
        dataset_label=dataset_label,
        workflow_status=workflow_status,
        variant_id=str(task_row["variant_id"]),
        label=_maybe_path(task_row.get("label")),
        execution_mode=_maybe_path(task_row.get("execution_mode")),
        task_status=task_status,
        job_status=job_status,
        output_freshness_status=freshness_status,
        recovery_action=recovery_action,
        merge_status=merge_status,
        evidence_status=evidence_status,
        reproducibility_status=reproducibility_status,
        selected_model=_maybe_path(variant_summary_row.get("selected_model")),
        output_root=_maybe_path(job_status_row.get("output_root"))
        or _maybe_path(
            task_record.get("output_root") if isinstance(task_record, dict) else None
        ),
        task_log_path=_maybe_path(task_row.get("log_path"))
        or _maybe_path(job_status_row.get("task_log_path")),
        evidence_json_path=_maybe_path(merge_row.get("evidence_json_path"))
        or _maybe_path(evidence_row.get("evidence_json_path")),
        evidence_html_path=_maybe_path(merge_row.get("evidence_html_path"))
        or _maybe_path(evidence_row.get("evidence_html_path")),
        variant_count=None,
        successful_variant_count=None,
        failed_variant_count=None,
        output_file_count=_optional_int(job_status_row.get("output_file_count")),
        output_byte_count=_optional_int(job_status_row.get("output_byte_count")),
        artifact_file_count=_optional_int(evidence_row.get("artifact_file_count")),
        linked_artifact_count=None,
        linked_artifact_bytes=None,
        issue_count=len(issues),
        issues=issues,
        error_code=error_code,
        error_message=error_message,
        job_evidence_warning_count=job_evidence_warning_count,
        warning_count=len(warnings),
        warnings=warnings,
    )


def _build_batch_dataset_row(
    *,
    bundle_root: Path,
    dataset_id: str,
    dataset_label: str,
    workflow_status: str,
    execution_mode: str,
    workflow_run: dict[str, object],
    workflow_status_summary: dict[str, object],
    failure_recovery_summary: dict[str, object],
    merge_summary: dict[str, object],
    reproducibility_summary: dict[str, object],
    report_manifest: dict[str, object],
) -> SupplementaryBatchSummaryRow:
    total_file_count, total_byte_count = _bundle_output_totals(bundle_root)
    warnings: list[str] = []
    issues: list[str] = []
    active_run_state = str(workflow_status_summary.get("active_run_state", "unknown"))
    overall_recovery_status = str(
        failure_recovery_summary.get("overall_recovery_status", "unknown")
    )
    merge_status = str(merge_summary.get("merge_status", "unknown"))
    reproducibility_status = (
        "passed" if bool(reproducibility_summary.get("all_passed", False)) else "failed"
    )
    failed_variant_count = len(list(workflow_run.get("failed_variants", [])))
    successful_variant_count = len(list(workflow_run.get("successful_variants", [])))
    if workflow_status != "succeeded":
        warnings.append(f"workflow-status:{workflow_status}")
    if active_run_state != "absent":
        warnings.append(f"active-run-state:{active_run_state}")
    if overall_recovery_status != "clean":
        warnings.append(f"recovery-status:{overall_recovery_status}")
    if merge_status != "merge-ready":
        warnings.append(f"merge-status:{merge_status}")
    if reproducibility_status != "passed":
        warnings.append("reproducibility-status:failed")
    failed_job_count = int(workflow_status_summary.get("failed_job_count", 0))
    stale_job_count = int(workflow_status_summary.get("stale_job_count", 0))
    stale_output_job_count = int(
        workflow_status_summary.get("stale_output_job_count", 0)
    )
    if failed_job_count > 0:
        issues.append(f"failed-job-count:{failed_job_count}")
    if stale_job_count > 0:
        issues.append(f"stale-job-count:{stale_job_count}")
    if stale_output_job_count > 0:
        issues.append(f"stale-output-job-count:{stale_output_job_count}")
    if int(reproducibility_summary.get("failed_variant_count", 0)) > 0:
        issues.append(
            "reproducibility-failed-variant-count:"
            + str(reproducibility_summary["failed_variant_count"])
        )
    warnings.extend(issue for issue in issues if issue not in warnings)
    linked_artifacts = report_manifest.get("linked_artifacts", {})
    linked_artifact_bytes = 0
    if isinstance(linked_artifacts, dict):
        linked_artifact_bytes = sum(
            int(payload.get("byte_count", 0))
            for payload in linked_artifacts.values()
            if isinstance(payload, dict)
        )
    return SupplementaryBatchSummaryRow(
        row_scope="dataset",
        dataset_id=dataset_id,
        dataset_label=dataset_label,
        workflow_status=workflow_status,
        variant_id=None,
        label=dataset_label,
        execution_mode=execution_mode,
        task_status=None,
        job_status=None,
        output_freshness_status=None,
        recovery_action=overall_recovery_status,
        merge_status=merge_status,
        evidence_status=None,
        reproducibility_status=reproducibility_status,
        selected_model=None,
        output_root=".",
        task_log_path=None,
        evidence_json_path=None,
        evidence_html_path=None,
        variant_count=int(workflow_run.get("variant_count", 0)),
        successful_variant_count=successful_variant_count,
        failed_variant_count=failed_variant_count,
        output_file_count=total_file_count,
        output_byte_count=total_byte_count,
        artifact_file_count=None,
        linked_artifact_count=int(report_manifest.get("linked_artifact_count", 0)),
        linked_artifact_bytes=linked_artifact_bytes,
        issue_count=len(issues),
        issues=issues,
        error_code=None,
        error_message=None,
        job_evidence_warning_count=None,
        warning_count=len(warnings),
        warnings=warnings,
    )


def write_supplementary_batch_summary_table(
    path: Path,
    *,
    workflow_bundle_root: Path,
) -> SupplementaryBatchSummaryTableResult:
    """Write one supplementary batch summary table from a written workflow bundle."""
    bundle_root = workflow_bundle_root.resolve()
    workflow_run = _read_bundle_json(
        bundle_root / "rabies-method-sensitivity-panel.run.json"
    )
    workflow_manifest = _read_bundle_json(
        bundle_root / "rabies-method-sensitivity.manifest.json"
    )
    report_manifest = _read_bundle_json(
        bundle_root
        / "report-artifacts"
        / "rabies-method-sensitivity-report.manifest.json"
    )
    workflow_summary_rows = _read_bundle_rows(bundle_root / "workflow-summary.tsv")
    task_rows = _read_bundle_rows(bundle_root / "parallel-execution-summary.tsv")
    variant_summary_rows = _read_bundle_rows(bundle_root / "variant-summary.tsv")
    job_status_rows = _read_bundle_rows(bundle_root / "slurm-job-status.tsv")
    freshness_rows = _read_bundle_rows(bundle_root / "slurm-output-freshness.tsv")
    recovery_rows = _read_bundle_rows(bundle_root / "slurm-failure-recovery-jobs.tsv")
    merge_rows = _read_bundle_rows(bundle_root / "slurm-merge-variants.tsv")
    evidence_rows = _read_bundle_rows(bundle_root / "slurm-job-evidence.tsv")
    reproducibility_rows = _read_bundle_rows(
        bundle_root / "reproducibility-variants.tsv"
    )
    workflow_status_summary = _read_bundle_json(
        bundle_root / "slurm-workflow-status.json"
    )
    failure_recovery_summary = _read_bundle_json(
        bundle_root / "slurm-failure-recovery-report.json"
    )
    merge_summary = _read_bundle_json(bundle_root / "slurm-merge-report.json")
    reproducibility_summary = _read_bundle_json(
        bundle_root / "reproducibility-audit.json"
    )

    workflow_summary = workflow_summary_rows[0]
    dataset_id = str(workflow_summary["dataset_id"])
    dataset_label = str(workflow_manifest.get("label", dataset_id))
    task_records = workflow_run.get("task_records", [])
    task_record_by_variant = {
        str(record["variant_id"]): record
        for record in task_records
        if isinstance(record, dict)
    }
    task_row_by_variant = {str(row["variant_id"]): row for row in task_rows}
    variant_summary_by_variant = {
        str(row["variant_id"]): row for row in variant_summary_rows
    }
    job_status_by_variant = {str(row["variant_id"]): row for row in job_status_rows}
    freshness_by_variant = {str(row["variant_id"]): row for row in freshness_rows}
    recovery_by_variant = {str(row["variant_id"]): row for row in recovery_rows}
    merge_by_variant = {str(row["variant_id"]): row for row in merge_rows}
    evidence_by_variant = {str(row["variant_id"]): row for row in evidence_rows}
    reproducibility_by_variant = {
        str(row["variant_id"]): row for row in reproducibility_rows
    }
    workflow_status = str(workflow_run.get("status", "unknown"))
    execution_mode = str(workflow_run.get("execution_mode", "unknown"))
    variant_ids = [
        str(record["variant_id"])
        for record in task_records
        if isinstance(record, dict) and "variant_id" in record
    ]
    rows = [
        _build_batch_dataset_row(
            bundle_root=bundle_root,
            dataset_id=dataset_id,
            dataset_label=dataset_label,
            workflow_status=workflow_status,
            execution_mode=execution_mode,
            workflow_run=workflow_run,
            workflow_status_summary=workflow_status_summary,
            failure_recovery_summary=failure_recovery_summary,
            merge_summary=merge_summary,
            reproducibility_summary=reproducibility_summary,
            report_manifest=report_manifest,
        )
    ]
    rows.extend(
        _build_batch_variant_row(
            dataset_id=dataset_id,
            dataset_label=dataset_label,
            workflow_status=workflow_status,
            task_row=task_row_by_variant[variant_id],
            task_record=task_record_by_variant[variant_id],
            variant_summary_row=variant_summary_by_variant[variant_id],
            job_status_row=job_status_by_variant[variant_id],
            freshness_row=freshness_by_variant[variant_id],
            recovery_row=recovery_by_variant[variant_id],
            merge_row=merge_by_variant[variant_id],
            evidence_row=evidence_by_variant[variant_id],
            reproducibility_row=reproducibility_by_variant[variant_id],
        )
        for variant_id in variant_ids
    )
    columns = batch_summary_table_columns()
    _write_batch_summary_rows(path, columns=columns, rows=rows)
    return SupplementaryBatchSummaryTableResult(
        output_path=path,
        row_count=len(rows),
        dataset_row_count=sum(1 for row in rows if row.row_scope == "dataset"),
        variant_row_count=sum(1 for row in rows if row.row_scope == "variant"),
        workflow_status=workflow_status,
        warning_count=sum(row.warning_count for row in rows),
        columns=columns,
        rows=rows,
    )

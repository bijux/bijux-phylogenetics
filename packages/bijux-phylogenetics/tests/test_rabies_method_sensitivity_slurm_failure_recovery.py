from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.datasets.rabies_method_sensitivity.slurm.failure_recovery import (
    build_rabies_method_sensitivity_slurm_failure_recovery_report,
    write_rabies_method_sensitivity_slurm_failure_recovery_html_report,
    write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table,
    write_rabies_method_sensitivity_slurm_failure_recovery_summary_json,
)


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _write_tsv(
    path: Path,
    *,
    fieldnames: tuple[str, ...],
    rows: list[dict[str, object]],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _write_task_log(
    path: Path,
    *,
    variant_id: str,
    status: str,
    error_code: str | None = None,
    error_message: str | None = None,
) -> None:
    lines = [
        f"variant_id: {variant_id}",
        "execution_mode: parallel",
        f"status: {status}",
        f"output_root: variants/{variant_id}",
    ]
    if error_code is not None:
        lines.append(f"error_code: {error_code}")
    if error_message is not None:
        lines.append(f"error_message: {error_message}")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_synthetic_bundle(bundle_root: Path) -> None:
    _write_json(
        bundle_root / "workflow-config.resolved.json",
        {
            "dataset_id": "rabies_method_sensitivity_panel",
            "workflow_prefix": "rabies-method-sensitivity-panel",
            "variants": [
                {"variant_id": "auto-gap-threshold"},
                {"variant_id": "auto-gappyout"},
                {"variant_id": "ginsi-gappyout"},
            ],
        },
    )
    _write_tsv(
        bundle_root / "slurm-job-status.tsv",
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
                "partition_id": "compact-mafft-auto-standard",
                "array_index": 0,
                "variant_id": "auto-gap-threshold",
                "status": "failed",
                "evidence_class": "terminal-failure",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "task_log_path": "parallel-logs/auto-gap-threshold.log",
                "output_root": "variants/auto-gap-threshold",
                "task_status": "failed",
                "execution_record_status": "failed",
                "output_file_count": 3,
                "output_byte_count": 1024,
                "missing_required_file_count": 1,
                "missing_required_files": "rooted-fasttree.nwk",
                "detail": "task log or execution record marks the job as failed",
                "output_freshness_status": "unknown",
                "output_freshness_reason_codes": "",
            },
            {
                "partition_id": "compact-mafft-auto-standard",
                "array_index": 1,
                "variant_id": "auto-gappyout",
                "status": "stale",
                "evidence_class": "stale-output-drift",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "task_log_path": "parallel-logs/auto-gappyout.log",
                "output_root": "variants/auto-gappyout",
                "task_status": "succeeded",
                "execution_record_status": "succeeded",
                "output_file_count": 12,
                "output_byte_count": 8192,
                "missing_required_file_count": 0,
                "missing_required_files": "",
                "detail": "outputs are stale against current inputs",
                "output_freshness_status": "stale",
                "output_freshness_reason_codes": "input_checksum_drift",
            },
            {
                "partition_id": "compact-mafft-ginsi-elevated",
                "array_index": 0,
                "variant_id": "ginsi-gappyout",
                "status": "pending",
                "evidence_class": "planned-not-started",
                "script_path": "slurm-arrays/compact-mafft-ginsi-elevated.sbatch",
                "task_log_path": "parallel-logs/ginsi-gappyout.log",
                "output_root": "variants/ginsi-gappyout",
                "task_status": "",
                "execution_record_status": "",
                "output_file_count": 0,
                "output_byte_count": 0,
                "missing_required_file_count": 12,
                "missing_required_files": "fasttree.nwk",
                "detail": "workflow is still running and this planned job has not written durable evidence yet",
                "output_freshness_status": "unknown",
                "output_freshness_reason_codes": "",
            },
        ],
    )
    _write_tsv(
        bundle_root / "slurm-partition-status.tsv",
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
                "partition_id": "compact-mafft-auto-standard",
                "script_path": "slurm-arrays/compact-mafft-auto-standard.sbatch",
                "job_count": 2,
                "completed_job_count": 0,
                "failed_job_count": 1,
                "pending_job_count": 0,
                "stale_job_count": 1,
                "overall_status": "failed",
                "variant_ids": "auto-gap-threshold,auto-gappyout",
            },
            {
                "partition_id": "compact-mafft-ginsi-elevated",
                "script_path": "slurm-arrays/compact-mafft-ginsi-elevated.sbatch",
                "job_count": 1,
                "completed_job_count": 0,
                "failed_job_count": 0,
                "pending_job_count": 1,
                "stale_job_count": 0,
                "overall_status": "pending",
                "variant_ids": "ginsi-gappyout",
            },
        ],
    )
    _write_json(
        bundle_root / "slurm-workflow-status.json",
        {
            "workflow_status": "failed",
            "active_run_state": "live",
        },
    )
    _write_task_log(
        bundle_root / "parallel-logs" / "auto-gap-threshold.log",
        variant_id="auto-gap-threshold",
        status="failed",
        error_code="engine_process_timeout",
        error_message="engine timeout after 3600 seconds",
    )
    _write_task_log(
        bundle_root / "parallel-logs" / "auto-gappyout.log",
        variant_id="auto-gappyout",
        status="succeeded",
    )


def test_build_rabies_method_sensitivity_slurm_failure_recovery_report_identifies_reruns(
    tmp_path: Path,
) -> None:
    _write_synthetic_bundle(tmp_path)

    report = build_rabies_method_sensitivity_slurm_failure_recovery_report(tmp_path)

    assert report.dataset_id == "rabies_method_sensitivity_panel"
    assert report.overall_recovery_status == "recovery_needed"
    assert report.rerunnable_job_count == 2
    assert report.blocked_job_count == 1
    assert report.recovery_partition_count == 2
    failed_row = next(
        row for row in report.jobs if row.variant_id == "auto-gap-threshold"
    )
    assert failed_row.rerunnable is True
    assert failed_row.recovery_action == "rerun_variant"
    assert failed_row.likely_cause_code == "task_timeout"
    stale_row = next(row for row in report.jobs if row.variant_id == "auto-gappyout")
    assert stale_row.recovery_action == "refresh_inputs_then_rerun"
    pending_row = next(row for row in report.jobs if row.variant_id == "ginsi-gappyout")
    assert pending_row.rerunnable is False
    assert pending_row.recovery_action == "wait_for_live_workflow"


def test_build_rabies_method_sensitivity_slurm_failure_recovery_report_is_clean_for_packaged_outputs() -> (
    None
):
    from bijux_phylogenetics.datasets import (
        load_rabies_method_sensitivity_panel_dataset,
    )

    dataset = load_rabies_method_sensitivity_panel_dataset()
    report = build_rabies_method_sensitivity_slurm_failure_recovery_report(
        dataset.reference_output_root
    )

    assert report.job_count == 4
    assert report.overall_recovery_status == "clean"
    assert report.rerunnable_job_count == 0
    assert report.blocked_job_count == 0


def test_write_rabies_method_sensitivity_slurm_failure_recovery_artifacts(
    tmp_path: Path,
) -> None:
    _write_synthetic_bundle(tmp_path)
    report = build_rabies_method_sensitivity_slurm_failure_recovery_report(tmp_path)

    jobs_path = write_rabies_method_sensitivity_slurm_failure_recovery_jobs_table(
        tmp_path / "slurm-failure-recovery-jobs.tsv",
        report,
    )
    partitions_path = (
        write_rabies_method_sensitivity_slurm_failure_recovery_partitions_table(
            tmp_path / "slurm-failure-recovery-partitions.tsv",
            report,
        )
    )
    summary_path = write_rabies_method_sensitivity_slurm_failure_recovery_summary_json(
        tmp_path / "slurm-failure-recovery-report.json",
        report,
    )
    html_path = write_rabies_method_sensitivity_slurm_failure_recovery_html_report(
        tmp_path / "slurm-failure-recovery-report.html",
        report,
    )

    assert "recovery_action" in jobs_path.read_text(encoding="utf-8")
    assert "rerunnable_job_count" in partitions_path.read_text(encoding="utf-8")
    payload = json.loads(summary_path.read_text(encoding="utf-8"))
    assert payload["overall_recovery_status"] == "recovery_needed"
    html = html_path.read_text(encoding="utf-8")
    assert "Rabies Slurm Failure Recovery Report" in html
    assert "slurm-failure-recovery-report.json" in html

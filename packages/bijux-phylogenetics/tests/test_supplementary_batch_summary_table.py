from __future__ import annotations

import csv
import json
from pathlib import Path
import shutil

from bijux_phylogenetics.reports import write_supplementary_batch_summary_table

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")
RABIES_EXPECTED_BUNDLE = (
    Path(__file__).resolve().parent.parent
    / "src"
    / "bijux_phylogenetics"
    / "resources"
    / "datasets"
    / "pathogens"
    / "rabies_method_sensitivity_panel"
    / "expected"
)


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def _read_tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_tsv_rows(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def _copy_expected_bundle(tmp_path: Path) -> Path:
    destination = tmp_path / "workflow-bundle"
    shutil.copytree(RABIES_EXPECTED_BUNDLE, destination)
    return destination


def test_write_supplementary_batch_summary_table_reads_governed_bundle(
    tmp_path: Path,
) -> None:
    output_path = tmp_path / "supplementary-batch-summary.tsv"

    result = write_supplementary_batch_summary_table(
        output_path,
        workflow_bundle_root=RABIES_EXPECTED_BUNDLE,
    )

    assert result.output_path == output_path
    assert result.row_count == 5
    assert result.dataset_row_count == 1
    assert result.variant_row_count == 4
    assert result.workflow_status == "succeeded"
    dataset_row = next(row for row in result.rows if row.row_scope == "dataset")
    assert dataset_row.dataset_id == "rabies_method_sensitivity_panel"
    assert dataset_row.workflow_status == "succeeded"
    assert dataset_row.variant_count == 4
    assert dataset_row.successful_variant_count == 4
    assert dataset_row.failed_variant_count == 0
    assert dataset_row.linked_artifact_count == 32
    variant_rows = [row for row in result.rows if row.row_scope == "variant"]
    assert all(row.output_file_count == 14 for row in variant_rows)
    assert all(row.output_byte_count == 42065 for row in variant_rows)
    assert all(row.selected_model == "TPM2u+F+G4" for row in variant_rows)
    assert all(row.job_status == "completed" for row in variant_rows)
    assert all(row.output_freshness_status == "fresh" for row in variant_rows)
    assert all("job-evidence-warning-count:4" in row.warnings for row in variant_rows)
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 5
    assert rows[0]["row_scope"] == "dataset"
    assert rows[1]["row_scope"] == "variant"


def test_write_supplementary_batch_summary_table_surfaces_failed_variant_bundle_state(
    tmp_path: Path,
) -> None:
    bundle_root = _copy_expected_bundle(tmp_path)
    output_path = tmp_path / "supplementary-batch-summary-failed.tsv"

    run_path = bundle_root / "rabies-method-sensitivity-panel.run.json"
    run_payload = json.loads(run_path.read_text(encoding="utf-8"))
    run_payload["status"] = "failed"
    run_payload["successful_variants"] = [
        variant
        for variant in run_payload["successful_variants"]
        if variant != "auto-gap-threshold"
    ]
    run_payload["failed_variants"] = ["auto-gap-threshold"]
    run_payload["task_records"][0]["status"] = "failed"
    run_payload["task_records"][0]["error_code"] = "parallel_variant_failed"
    run_payload["task_records"][0]["error_message"] = "simulated batch failure"
    run_path.write_text(json.dumps(run_payload, indent=2) + "\n", encoding="utf-8")

    workflow_status_path = bundle_root / "slurm-workflow-status.json"
    workflow_status_payload = json.loads(
        workflow_status_path.read_text(encoding="utf-8")
    )
    workflow_status_payload["workflow_status"] = "failed"
    workflow_status_payload["failed_job_count"] = 1
    workflow_status_payload["fresh_output_job_count"] = 3
    workflow_status_payload["stale_output_job_count"] = 1
    workflow_status_path.write_text(
        json.dumps(workflow_status_payload, indent=2) + "\n",
        encoding="utf-8",
    )

    recovery_summary_path = bundle_root / "slurm-failure-recovery-report.json"
    recovery_summary_payload = json.loads(
        recovery_summary_path.read_text(encoding="utf-8")
    )
    recovery_summary_payload["workflow_status"] = "failed"
    recovery_summary_payload["overall_recovery_status"] = "recovery_needed"
    recovery_summary_payload["rerunnable_job_count"] = 1
    recovery_summary_payload["failed_job_count"] = 1
    recovery_summary_path.write_text(
        json.dumps(recovery_summary_payload, indent=2) + "\n",
        encoding="utf-8",
    )

    merge_summary_path = bundle_root / "slurm-merge-report.json"
    merge_summary_payload = json.loads(merge_summary_path.read_text(encoding="utf-8"))
    merge_summary_payload["merge_status"] = "merge-blocked"
    merge_summary_payload["merge_ready"] = False
    merge_summary_payload["failed_variant_count"] = 1
    merge_summary_payload["mergeable_variant_count"] = 3
    merge_summary_path.write_text(
        json.dumps(merge_summary_payload, indent=2) + "\n",
        encoding="utf-8",
    )

    reproducibility_summary_path = bundle_root / "reproducibility-audit.json"
    reproducibility_summary_payload = json.loads(
        reproducibility_summary_path.read_text(encoding="utf-8")
    )
    reproducibility_summary_payload["all_passed"] = False
    reproducibility_summary_payload["failed_variant_count"] = 1
    reproducibility_summary_path.write_text(
        json.dumps(reproducibility_summary_payload, indent=2) + "\n",
        encoding="utf-8",
    )

    parallel_rows = _read_tsv_rows(bundle_root / "parallel-execution-summary.tsv")
    parallel_rows[0]["status"] = "failed"
    parallel_rows[0]["error_code"] = "parallel_variant_failed"
    _write_tsv_rows(bundle_root / "parallel-execution-summary.tsv", parallel_rows)

    job_status_rows = _read_tsv_rows(bundle_root / "slurm-job-status.tsv")
    job_status_rows[0]["status"] = "failed"
    job_status_rows[0]["task_status"] = "failed"
    job_status_rows[0]["execution_record_status"] = "failed"
    job_status_rows[0]["output_freshness_status"] = "stale"
    job_status_rows[0]["output_freshness_reason_codes"] = "input_digest_changed"
    _write_tsv_rows(bundle_root / "slurm-job-status.tsv", job_status_rows)

    freshness_rows = _read_tsv_rows(bundle_root / "slurm-output-freshness.tsv")
    freshness_rows[0]["freshness_status"] = "stale"
    freshness_rows[0]["inputs_match"] = "false"
    freshness_rows[0]["stale_reason_count"] = "1"
    freshness_rows[0]["stale_reason_codes"] = "input_digest_changed"
    freshness_rows[0]["stale_reason_detail"] = "current inputs no longer match"
    _write_tsv_rows(bundle_root / "slurm-output-freshness.tsv", freshness_rows)

    recovery_rows = _read_tsv_rows(bundle_root / "slurm-failure-recovery-jobs.tsv")
    recovery_rows[0]["current_status"] = "failed"
    recovery_rows[0]["rerunnable"] = "true"
    recovery_rows[0]["recovery_action"] = "rerun_variant"
    recovery_rows[0]["likely_cause_code"] = "task_failure"
    recovery_rows[0]["task_status"] = "failed"
    recovery_rows[0]["error_code"] = "parallel_variant_failed"
    recovery_rows[0]["error_message"] = "simulated batch failure"
    _write_tsv_rows(bundle_root / "slurm-failure-recovery-jobs.tsv", recovery_rows)

    merge_rows = _read_tsv_rows(bundle_root / "slurm-merge-variants.tsv")
    merge_rows[0]["merge_status"] = "blocked"
    merge_rows[0]["included_in_merge"] = "false"
    merge_rows[0]["issue_count"] = "2"
    merge_rows[0]["issues"] = "output-not-fresh|task-failed"
    _write_tsv_rows(bundle_root / "slurm-merge-variants.tsv", merge_rows)

    reproducibility_rows = _read_tsv_rows(bundle_root / "reproducibility-variants.tsv")
    reproducibility_rows[0]["status"] = "failed"
    reproducibility_rows[0]["issues"] = "output_digest_changed"
    _write_tsv_rows(bundle_root / "reproducibility-variants.tsv", reproducibility_rows)

    result = write_supplementary_batch_summary_table(
        output_path,
        workflow_bundle_root=bundle_root,
    )

    dataset_row = next(row for row in result.rows if row.row_scope == "dataset")
    failed_variant_row = next(
        row for row in result.rows if row.variant_id == "auto-gap-threshold"
    )
    assert result.workflow_status == "failed"
    assert dataset_row.workflow_status == "failed"
    assert "workflow-status:failed" in dataset_row.warnings
    assert "recovery-status:recovery_needed" in dataset_row.warnings
    assert "merge-status:merge-blocked" in dataset_row.warnings
    assert "reproducibility-status:failed" in dataset_row.warnings
    assert failed_variant_row.task_status == "failed"
    assert failed_variant_row.job_status == "failed"
    assert failed_variant_row.output_freshness_status == "stale"
    assert failed_variant_row.recovery_action == "rerun_variant"
    assert failed_variant_row.merge_status == "blocked"
    assert failed_variant_row.reproducibility_status == "failed"
    assert failed_variant_row.error_code == "parallel_variant_failed"
    assert failed_variant_row.error_message == "simulated batch failure"
    assert "task-status:failed" in failed_variant_row.warnings
    assert "task-error:parallel_variant_failed" in failed_variant_row.warnings
    assert "freshness-status:stale" in failed_variant_row.warnings
    assert "recovery-action:rerun_variant" in failed_variant_row.warnings
    assert "merge-status:blocked" in failed_variant_row.warnings
    assert "reproducibility-status:failed" in failed_variant_row.warnings
    assert "output-not-fresh" in failed_variant_row.issues
    assert "output_digest_changed" in failed_variant_row.issues

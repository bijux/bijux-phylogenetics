from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import (
    RabiesMethodSensitivityOutputFreshnessCheckRow,
    RabiesMethodSensitivitySlurmOutputFreshnessReport,
    RabiesMethodSensitivitySlurmOutputFreshnessRow,
)
from .interfaces import DatasetLike
from .policy import build_freshness_row, evaluate_freshness_checks
from .shared import (
    CONFIG_FILENAME,
    SLURM_ARRAY_MEMBERS_FILENAME,
    load_json,
    read_tsv_rows,
    sha256,
    write_tsv,
)

__all__ = [
    "RabiesMethodSensitivityOutputFreshnessCheckRow",
    "RabiesMethodSensitivitySlurmOutputFreshnessRow",
    "RabiesMethodSensitivitySlurmOutputFreshnessReport",
    "build_rabies_method_sensitivity_slurm_output_freshness_report",
    "write_rabies_method_sensitivity_slurm_output_freshness_checks_table",
    "write_rabies_method_sensitivity_slurm_output_freshness_json",
    "write_rabies_method_sensitivity_slurm_output_freshness_table",
]

def build_rabies_method_sensitivity_slurm_output_freshness_report(
    bundle_root: Path,
    *,
    dataset: DatasetLike | None = None,
) -> RabiesMethodSensitivitySlurmOutputFreshnessReport:
    """Detect whether bundle outputs still match the current packaged workflow state."""
    bundle_root = bundle_root.resolve()
    resolved_config = load_json(bundle_root / CONFIG_FILENAME)
    if dataset is None:
        from ...config import (
            load_rabies_method_sensitivity_panel_dataset,
        )

        dataset = load_rabies_method_sensitivity_panel_dataset()
    member_rows = read_tsv_rows(bundle_root / SLURM_ARRAY_MEMBERS_FILENAME)
    selected_variant_ids = tuple(
        str(value) for value in list(resolved_config.get("selected_variant_ids", []))
    )
    if not selected_variant_ids:
        selected_variant_ids = tuple(str(row["variant_id"]) for row in member_rows)
    input_checksums = {
        "sequences.fasta": sha256(dataset.sequences_path),
        "metadata.csv": sha256(dataset.metadata_path),
    }
    (
        checks,
        bundle_input_ok,
        bundle_workflow_settings_ok,
        variant_check_rows_by_variant,
    ) = evaluate_freshness_checks(
        resolved_config=resolved_config,
        dataset=dataset,
        selected_variant_ids=selected_variant_ids,
        input_checksums=input_checksums,
    )

    global_failed_reason_codes = tuple(
        row.check_id
        for row in checks
        if row.scope == "workflow" and row.status == "failed"
    )
    global_failed_reason_detail = "; ".join(
        row.detail
        for row in checks
        if row.scope == "workflow" and row.status == "failed"
    )
    freshness_rows = tuple(
        build_freshness_row(
            member_row=member_row,
            bundle_input_ok=bundle_input_ok,
            bundle_workflow_settings_ok=bundle_workflow_settings_ok,
            global_failed_reason_codes=global_failed_reason_codes,
            global_failed_reason_detail=global_failed_reason_detail,
            variant_check_rows=variant_check_rows_by_variant.get(
                str(member_row["variant_id"]),
                [],
            ),
        )
        for member_row in member_rows
    )
    failed_check_count = sum(1 for row in checks if row.status == "failed")
    fresh_job_count = sum(1 for row in freshness_rows if row.freshness_status == "fresh")
    stale_job_count = sum(1 for row in freshness_rows if row.freshness_status == "stale")

    return RabiesMethodSensitivitySlurmOutputFreshnessReport(
        dataset_id=dataset.dataset_id,
        workflow_prefix=dataset.workflow_prefix,
        bundle_root=bundle_root,
        all_outputs_fresh=failed_check_count == 0 and stale_job_count == 0,
        selected_variant_ids=selected_variant_ids,
        check_count=len(checks),
        failed_check_count=failed_check_count,
        fresh_job_count=fresh_job_count,
        stale_job_count=stale_job_count,
        checks=tuple(checks),
        jobs=freshness_rows,
    )


def write_rabies_method_sensitivity_slurm_output_freshness_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
) -> Path:
    """Write one per-job output-freshness ledger."""
    return write_tsv(
        path,
        fieldnames=(
            "partition_id",
            "array_index",
            "variant_id",
            "freshness_status",
            "inputs_match",
            "workflow_settings_match",
            "variant_settings_match",
            "stale_reason_count",
            "stale_reason_codes",
            "stale_reason_detail",
        ),
        rows=[
            {
                "partition_id": row.partition_id,
                "array_index": row.array_index,
                "variant_id": row.variant_id,
                "freshness_status": row.freshness_status,
                "inputs_match": str(row.inputs_match).lower(),
                "workflow_settings_match": str(row.workflow_settings_match).lower(),
                "variant_settings_match": str(row.variant_settings_match).lower(),
                "stale_reason_count": row.stale_reason_count,
                "stale_reason_codes": ",".join(row.stale_reason_codes),
                "stale_reason_detail": row.stale_reason_detail,
            }
            for row in report.jobs
        ],
    )


def write_rabies_method_sensitivity_slurm_output_freshness_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
) -> Path:
    """Write one check-level ledger for output freshness."""
    return write_tsv(
        path,
        fieldnames=(
            "check_id",
            "surface",
            "scope",
            "status",
            "expected",
            "observed",
            "detail",
        ),
        rows=[
            {
                "check_id": row.check_id,
                "surface": row.surface,
                "scope": row.scope,
                "status": row.status,
                "expected": row.expected,
                "observed": row.observed,
                "detail": row.detail,
            }
            for row in report.checks
        ],
    )


def write_rabies_method_sensitivity_slurm_output_freshness_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputFreshnessReport,
) -> Path:
    """Write the structured output-freshness summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path



from __future__ import annotations

from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmOutputFreshnessReport
from .interfaces import DatasetLike
from .policy import build_freshness_row, evaluate_freshness_checks
from .shared import (
    CONFIG_FILENAME,
    SLURM_ARRAY_MEMBERS_FILENAME,
    load_json,
    read_tsv_rows,
    sha256,
)


def build_rabies_method_sensitivity_slurm_output_freshness_report(
    bundle_root: Path,
    *,
    dataset: DatasetLike | None = None,
) -> RabiesMethodSensitivitySlurmOutputFreshnessReport:
    """Detect whether bundle outputs still match the current packaged workflow state."""
    bundle_root = bundle_root.resolve()
    resolved_config = load_json(bundle_root / CONFIG_FILENAME)
    if dataset is None:
        from ...config import load_rabies_method_sensitivity_panel_dataset

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
    fresh_job_count = sum(
        1 for row in freshness_rows if row.freshness_status == "fresh"
    )
    stale_job_count = sum(
        1 for row in freshness_rows if row.freshness_status == "stale"
    )

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

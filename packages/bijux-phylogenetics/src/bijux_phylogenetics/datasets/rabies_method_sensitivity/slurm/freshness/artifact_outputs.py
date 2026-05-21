from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmOutputFreshnessReport
from .shared import write_tsv


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

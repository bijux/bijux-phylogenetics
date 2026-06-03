from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmStorageReport
from .shared import _write_tsv


def write_rabies_method_sensitivity_slurm_storage_categories_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write the category-level retained-storage estimate."""
    return _write_tsv(
        path,
        fieldnames=(
            "category_id",
            "category_label",
            "variant_file_count",
            "workflow_file_count",
            "total_file_count",
            "variant_byte_count",
            "workflow_byte_count",
            "total_byte_count",
            "estimated_storage_mib",
            "detail",
        ),
        rows=[asdict(row) for row in report.categories],
    )


def write_rabies_method_sensitivity_slurm_storage_variants_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write one per-variant retained-storage ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "output_file_count",
            "log_file_count",
            "tree_file_count",
            "posterior_sample_file_count",
            "report_file_count",
            "total_file_count",
            "output_byte_count",
            "log_byte_count",
            "tree_byte_count",
            "posterior_sample_byte_count",
            "report_byte_count",
            "total_byte_count",
            "estimated_storage_mib",
        ),
        rows=[asdict(row) for row in report.variants],
    )


def write_rabies_method_sensitivity_slurm_storage_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write the machine-readable retained-storage summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path

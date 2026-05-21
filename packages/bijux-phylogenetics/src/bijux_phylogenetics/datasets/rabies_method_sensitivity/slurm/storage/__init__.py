from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    RabiesMethodSensitivitySlurmStorageAssumptionRow,
    RabiesMethodSensitivitySlurmStorageCategoryRow,
    RabiesMethodSensitivitySlurmStorageReport,
    RabiesMethodSensitivitySlurmStorageVariantRow,
)
from .report_builder import build_rabies_method_sensitivity_slurm_storage_report
from .shared import _write_tsv

__all__ = [
    "RabiesMethodSensitivitySlurmStorageAssumptionRow",
    "RabiesMethodSensitivitySlurmStorageCategoryRow",
    "RabiesMethodSensitivitySlurmStorageReport",
    "RabiesMethodSensitivitySlurmStorageVariantRow",
    "build_rabies_method_sensitivity_slurm_storage_report",
    "write_rabies_method_sensitivity_slurm_storage_categories_table",
    "write_rabies_method_sensitivity_slurm_storage_html_report",
    "write_rabies_method_sensitivity_slurm_storage_summary_json",
    "write_rabies_method_sensitivity_slurm_storage_variants_table",
]


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


def write_rabies_method_sensitivity_slurm_storage_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmStorageReport,
) -> Path:
    """Write the reviewer-facing retained-storage estimate."""
    category_lines = [
        (
            f"{row.category_id}: {row.total_file_count} files, "
            f"{row.total_byte_count} bytes, {row.estimated_storage_mib} MiB"
        )
        for row in report.categories
    ]
    largest_variant = next(
        row for row in report.variants if row.variant_id == report.largest_variant_id
    )
    return write_html_report(
        title="Rabies Slurm Storage Report",
        sections=[
            (
                "storage-summary",
                "\n".join(
                    [
                        f"variant_count: {report.variant_count}",
                        f"total_file_count: {report.total_file_count}",
                        f"total_byte_count: {report.total_byte_count}",
                        f"total_estimated_storage_mib: {report.total_estimated_storage_mib}",
                        (
                            "variant_scoped_byte_count: "
                            f"{report.variant_scoped_byte_count}"
                        ),
                        (
                            "workflow_shared_byte_count: "
                            f"{report.workflow_shared_byte_count}"
                        ),
                    ]
                ),
            ),
            ("storage-categories", "\n".join(category_lines)),
            (
                "largest-variant",
                "\n".join(
                    [
                        f"variant_id: {largest_variant.variant_id}",
                        f"total_file_count: {largest_variant.total_file_count}",
                        f"total_byte_count: {largest_variant.total_byte_count}",
                        (
                            "estimated_storage_mib: "
                            f"{largest_variant.estimated_storage_mib}"
                        ),
                    ]
                ),
            ),
            (
                "posterior-samples",
                (
                    "this governed workflow writes no posterior chains or posterior "
                    "tree sets, so posterior_samples remains an explicit zero-valued "
                    "storage category"
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "variant_count": report.variant_count,
            "total_byte_count": report.total_byte_count,
            "total_estimated_storage_mib": report.total_estimated_storage_mib,
            "largest_variant_id": report.largest_variant_id,
            "largest_variant_total_byte_count": report.largest_variant_total_byte_count,
            "output_byte_count": report.output_byte_count,
            "log_byte_count": report.log_byte_count,
            "tree_byte_count": report.tree_byte_count,
            "posterior_sample_byte_count": report.posterior_sample_byte_count,
            "report_byte_count": report.report_byte_count,
        },
        summary_metrics=[
            ("estimated storage MiB", report.total_estimated_storage_mib),
            ("total files", report.total_file_count),
            ("workflow outputs bytes", report.output_byte_count),
            ("log bytes", report.log_byte_count),
            ("tree bytes", report.tree_byte_count),
            ("posterior bytes", report.posterior_sample_byte_count),
            ("report bytes", report.report_byte_count),
        ],
        artifact_links=[
            ("storage categories", "slurm-storage-categories.tsv", None),
            ("storage variants", "slurm-storage-variants.tsv", None),
            ("storage summary", "slurm-storage-report.json", None),
        ],
    )

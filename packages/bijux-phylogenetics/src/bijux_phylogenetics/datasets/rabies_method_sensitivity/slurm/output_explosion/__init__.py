from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from .contracts import (
    RabiesMethodSensitivitySlurmOutputExplosionCheckRow,
    RabiesMethodSensitivitySlurmOutputExplosionReport,
    RabiesMethodSensitivitySlurmOutputExplosionVariantRow,
)
from .report_builder import (
    build_rabies_method_sensitivity_slurm_output_explosion_report,
)
from .shared import (
    _format_float,
    _write_tsv,
)

__all__ = [
    "RabiesMethodSensitivitySlurmOutputExplosionCheckRow",
    "RabiesMethodSensitivitySlurmOutputExplosionReport",
    "RabiesMethodSensitivitySlurmOutputExplosionVariantRow",
    "build_rabies_method_sensitivity_slurm_output_explosion_report",
    "write_rabies_method_sensitivity_slurm_output_explosion_checks_table",
    "write_rabies_method_sensitivity_slurm_output_explosion_html_report",
    "write_rabies_method_sensitivity_slurm_output_explosion_summary_json",
    "write_rabies_method_sensitivity_slurm_output_explosion_variants_table",
]


def write_rabies_method_sensitivity_slurm_output_explosion_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the check-level output explosion ledger."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_slurm_output_explosion_variants_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the per-variant output explosion assessment."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "risk_status",
            "estimated_output_mib",
            "estimated_storage_mib",
            "tree_file_count",
            "tree_byte_count",
            "posterior_sample_file_count",
            "posterior_sample_byte_count",
            "report_byte_count",
            "output_share",
            "issue_count",
            "issues",
        ),
        rows=[
            {
                "variant_id": row.variant_id,
                "risk_status": row.risk_status,
                "estimated_output_mib": row.estimated_output_mib,
                "estimated_storage_mib": row.estimated_storage_mib,
                "tree_file_count": row.tree_file_count,
                "tree_byte_count": row.tree_byte_count,
                "posterior_sample_file_count": row.posterior_sample_file_count,
                "posterior_sample_byte_count": row.posterior_sample_byte_count,
                "report_byte_count": row.report_byte_count,
                "output_share": _format_float(row.output_share),
                "issue_count": row.issue_count,
                "issues": " | ".join(row.issues),
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_slurm_output_explosion_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the machine-readable output explosion report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_output_explosion_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmOutputExplosionReport,
) -> Path:
    """Write the reviewer-facing output explosion report."""
    return write_html_report(
        title="Rabies Slurm Output Explosion Report",
        sections=[
            (
                "risk-summary",
                "\n".join(
                    [
                        f"overall_risk_status: {report.overall_risk_status}",
                        f"bootstrap_replicates: {report.bootstrap_replicates}",
                        f"variant_count: {report.variant_count}",
                        (
                            "total_estimated_output_mib: "
                            f"{report.total_estimated_output_mib}"
                        ),
                        (
                            "total_estimated_storage_mib: "
                            f"{report.total_estimated_storage_mib}"
                        ),
                    ]
                ),
            ),
            (
                "global-issues",
                "none"
                if report.global_issue_count == 0
                else "\n".join(report.global_issues),
            ),
            (
                "high-risk-variants",
                "none"
                if report.high_risk_variant_count == 0
                else "\n".join(
                    f"{row.variant_id}: {'; '.join(row.issues)}"
                    for row in report.variants
                    if row.risk_status == "high"
                ),
            ),
            (
                "warning-variants",
                "none"
                if report.warning_variant_count == 0
                else "\n".join(
                    f"{row.variant_id}: {'; '.join(row.issues)}"
                    for row in report.variants
                    if row.risk_status == "warning"
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "overall_risk_status": report.overall_risk_status,
            "bootstrap_replicates": report.bootstrap_replicates,
            "variant_count": report.variant_count,
            "warning_variant_count": report.warning_variant_count,
            "high_risk_variant_count": report.high_risk_variant_count,
            "total_estimated_output_mib": report.total_estimated_output_mib,
            "total_estimated_storage_mib": report.total_estimated_storage_mib,
            "total_tree_byte_count": report.total_tree_byte_count,
            "total_posterior_sample_byte_count": report.total_posterior_sample_byte_count,
            "total_report_byte_count": report.total_report_byte_count,
            "largest_variant_id": report.largest_variant_id,
            "largest_variant_output_share": report.largest_variant_output_share,
        },
        summary_metrics=[
            ("overall risk", report.overall_risk_status),
            ("warning variants", report.warning_variant_count),
            ("high-risk variants", report.high_risk_variant_count),
            ("global issues", report.global_issue_count),
            ("estimated output MiB", report.total_estimated_output_mib),
            ("estimated storage MiB", report.total_estimated_storage_mib),
            ("posterior bytes", report.total_posterior_sample_byte_count),
        ],
        artifact_links=[
            ("output explosion checks", "slurm-output-explosion-checks.tsv", None),
            ("output explosion variants", "slurm-output-explosion-variants.tsv", None),
            ("output explosion summary", "slurm-output-explosion-report.json", None),
        ],
    )

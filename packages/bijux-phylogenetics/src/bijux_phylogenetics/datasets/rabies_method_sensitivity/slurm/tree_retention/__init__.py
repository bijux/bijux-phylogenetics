from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report
from .contracts import (
    RabiesMethodSensitivitySlurmTreeRetentionCheckRow,
    RabiesMethodSensitivitySlurmTreeRetentionFileRow,
    RabiesMethodSensitivitySlurmTreeRetentionReport,
)
from .report_builder import build_rabies_method_sensitivity_slurm_tree_retention_report
from .shared import _write_tsv

__all__ = [
    "RabiesMethodSensitivitySlurmTreeRetentionCheckRow",
    "RabiesMethodSensitivitySlurmTreeRetentionFileRow",
    "RabiesMethodSensitivitySlurmTreeRetentionReport",
    "build_rabies_method_sensitivity_slurm_tree_retention_report",
    "write_rabies_method_sensitivity_slurm_tree_retention_checks_table",
    "write_rabies_method_sensitivity_slurm_tree_retention_files_table",
    "write_rabies_method_sensitivity_slurm_tree_retention_html_report",
    "write_rabies_method_sensitivity_slurm_tree_retention_summary_json",
]


def write_rabies_method_sensitivity_slurm_tree_retention_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the check-level tree-retention ledger."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_slurm_tree_retention_files_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the per-file tree-retention recommendations."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "relative_path",
            "artifact_scope",
            "tree_count",
            "byte_count",
            "thinning_policy",
            "thinning_interval",
            "retained_tree_count",
            "compression_policy",
            "recommended_suffix",
            "issue_count",
            "issues",
        ),
        rows=[
            {
                "variant_id": row.variant_id,
                "relative_path": row.relative_path,
                "artifact_scope": row.artifact_scope,
                "tree_count": row.tree_count,
                "byte_count": row.byte_count,
                "thinning_policy": row.thinning_policy,
                "thinning_interval": row.thinning_interval,
                "retained_tree_count": row.retained_tree_count,
                "compression_policy": row.compression_policy,
                "recommended_suffix": row.recommended_suffix,
                "issue_count": row.issue_count,
                "issues": " | ".join(row.issues),
            }
            for row in report.files
        ],
    )


def write_rabies_method_sensitivity_slurm_tree_retention_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the machine-readable tree-retention policy report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_tree_retention_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmTreeRetentionReport,
) -> Path:
    """Write the reviewer-facing tree-retention policy report."""
    return write_html_report(
        title="Rabies Slurm Tree Retention Report",
        sections=[
            (
                "policy-summary",
                "\n".join(
                    [
                        f"overall_policy_status: {report.overall_policy_status}",
                        f"variant_count: {report.variant_count}",
                        f"file_count: {report.file_count}",
                        f"tree_set_file_count: {report.tree_set_file_count}",
                        (
                            "thinning_required_file_count: "
                            f"{report.thinning_required_file_count}"
                        ),
                        (
                            "compression_required_file_count: "
                            f"{report.compression_required_file_count}"
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
                "required-actions",
                "none"
                if report.overall_policy_status != "required"
                else "\n".join(
                    f"{row.relative_path}: {'; '.join(row.issues)}"
                    for row in report.files
                    if row.thinning_policy == "thin_required"
                    or row.compression_policy == "compress_required"
                ),
            ),
            (
                "recommended-actions",
                "none"
                if report.overall_policy_status == "no_action"
                else "\n".join(
                    f"{row.relative_path}: {'; '.join(row.issues)}"
                    for row in report.files
                    if row.issue_count > 0
                ),
            ),
        ],
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "overall_policy_status": report.overall_policy_status,
            "variant_count": report.variant_count,
            "file_count": report.file_count,
            "tree_set_file_count": report.tree_set_file_count,
            "posterior_sample_file_count": report.posterior_sample_file_count,
            "thinning_recommended_file_count": report.thinning_recommended_file_count,
            "thinning_required_file_count": report.thinning_required_file_count,
            "compression_recommended_file_count": report.compression_recommended_file_count,
            "compression_required_file_count": report.compression_required_file_count,
            "total_tree_count": report.total_tree_count,
            "total_tree_byte_count": report.total_tree_byte_count,
            "largest_tree_set_path": report.largest_tree_set_path,
            "largest_tree_set_tree_count": report.largest_tree_set_tree_count,
        },
        summary_metrics=[
            ("overall status", report.overall_policy_status),
            ("tree-set files", report.tree_set_file_count),
            ("posterior files", report.posterior_sample_file_count),
            ("thinning required", report.thinning_required_file_count),
            ("compression required", report.compression_required_file_count),
            ("total tree count", report.total_tree_count),
            ("total tree bytes", report.total_tree_byte_count),
        ],
        artifact_links=[
            ("tree retention checks", "slurm-tree-retention-checks.tsv", None),
            ("tree retention files", "slurm-tree-retention-files.tsv", None),
            ("tree retention summary", "slurm-tree-retention-policy.json", None),
        ],
    )

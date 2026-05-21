from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from bijux_phylogenetics.render.html import write_html_report

from .contracts import (
    RabiesMethodSensitivitySlurmMergeCheckRow,
    RabiesMethodSensitivitySlurmMergeReport,
    RabiesMethodSensitivitySlurmMergeVariantRow,
)
from .report_builder import build_rabies_method_sensitivity_slurm_merge_report
from .shared import (
    _write_tsv,
)

__all__ = [
    "RabiesMethodSensitivitySlurmMergeCheckRow",
    "RabiesMethodSensitivitySlurmMergeReport",
    "RabiesMethodSensitivitySlurmMergeVariantRow",
    "build_rabies_method_sensitivity_slurm_merge_report",
    "write_rabies_method_sensitivity_slurm_merge_checks_table",
    "write_rabies_method_sensitivity_slurm_merge_html_report",
    "write_rabies_method_sensitivity_slurm_merge_summary_json",
    "write_rabies_method_sensitivity_slurm_merge_variants_table",
]


def write_rabies_method_sensitivity_slurm_merge_checks_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write one check-level merge ledger."""
    return _write_tsv(
        path,
        fieldnames=("check_id", "surface", "status", "expected", "observed", "detail"),
        rows=[asdict(row) for row in report.checks],
    )


def write_rabies_method_sensitivity_slurm_merge_variants_table(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write one per-variant merge-decision ledger."""
    return _write_tsv(
        path,
        fieldnames=(
            "variant_id",
            "merge_status",
            "job_status",
            "output_freshness_status",
            "evidence_status",
            "included_in_merge",
            "selected_model",
            "serious_conflict_count",
            "rooted_engine_rf_distance",
            "rooted_engine_same_taxa_different_rooting",
            "issue_count",
            "issues",
            "evidence_json_path",
            "evidence_html_path",
        ),
        rows=[
            {
                "variant_id": row.variant_id,
                "merge_status": row.merge_status,
                "job_status": row.job_status,
                "output_freshness_status": row.output_freshness_status,
                "evidence_status": row.evidence_status,
                "included_in_merge": str(row.included_in_merge).lower(),
                "selected_model": row.selected_model,
                "serious_conflict_count": row.serious_conflict_count,
                "rooted_engine_rf_distance": row.rooted_engine_rf_distance,
                "rooted_engine_same_taxa_different_rooting": str(
                    row.rooted_engine_same_taxa_different_rooting
                ).lower(),
                "issue_count": row.issue_count,
                "issues": " | ".join(row.issues),
                "evidence_json_path": row.evidence_json_path,
                "evidence_html_path": row.evidence_html_path,
            }
            for row in report.variants
        ],
    )


def write_rabies_method_sensitivity_slurm_merge_summary_json(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write the structured merge summary."""
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = asdict(report)
    payload["bundle_root"] = "."
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return path


def write_rabies_method_sensitivity_slurm_merge_html_report(
    path: Path,
    report: RabiesMethodSensitivitySlurmMergeReport,
) -> Path:
    """Write the reviewer-facing global merge report."""
    sections = [
        (
            "merge-status",
            "\n".join(
                [
                    f"merge_status: {report.merge_status}",
                    f"merge_ready: {str(report.merge_ready).lower()}",
                    f"expected_variant_count: {report.expected_variant_count}",
                    f"merged_variant_count: {report.merged_variant_count}",
                    f"failed_variant_count: {report.failed_variant_count}",
                ]
            ),
        ),
        (
            "merged-science-summary",
            "\n".join(
                [
                    f"stable_clade_count: {report.stable_clade_count}",
                    f"changed_clade_count: {report.changed_clade_count}",
                    f"preprocessing_comparison_count: {report.preprocessing_comparison_count}",
                    f"conclusion_count: {report.conclusion_count}",
                    f"serious_conflict_variant_count: {report.serious_conflict_variant_count}",
                    f"rooted_engine_change_variant_count: {report.rooted_engine_change_variant_count}",
                    f"maximum_serious_conflict_count: {report.maximum_serious_conflict_count}",
                    f"selected_models: {', '.join(report.selected_models)}",
                ]
            ),
        ),
        (
            "failed-checks",
            "none"
            if report.failed_check_count == 0
            else "\n".join(
                f"{row.check_id}: {row.detail}"
                for row in report.checks
                if row.status == "failed"
            ),
        ),
        (
            "blocked-variants",
            "none"
            if report.failed_variant_count == 0
            else "\n".join(
                f"{row.variant_id}: {'; '.join(row.issues)}"
                for row in report.variants
                if not row.included_in_merge
            ),
        ),
    ]
    return write_html_report(
        title="Rabies Slurm Merge Report",
        sections=sections,
        out_path=path,
        embedded_json={
            "dataset_id": report.dataset_id,
            "workflow_prefix": report.workflow_prefix,
            "merge_status": report.merge_status,
            "merge_ready": report.merge_ready,
            "expected_variant_count": report.expected_variant_count,
            "merged_variant_count": report.merged_variant_count,
            "failed_variant_count": report.failed_variant_count,
            "failed_check_count": report.failed_check_count,
            "selected_models": list(report.selected_models),
        },
        summary_metrics=[
            ("merge status", report.merge_status),
            ("merge ready", str(report.merge_ready).lower()),
            ("merged variants", report.merged_variant_count),
            ("failed variants", report.failed_variant_count),
            ("failed checks", report.failed_check_count),
            ("stable clades", report.stable_clade_count),
            ("changed clades", report.changed_clade_count),
        ],
        artifact_links=[
            ("merge checks", "slurm-merge-checks.tsv", None),
            ("merge variants", "slurm-merge-variants.tsv", None),
            ("merge summary", "slurm-merge-report.json", None),
        ],
    )

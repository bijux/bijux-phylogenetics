from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmMergeReport
from .shared import _write_tsv


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

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmTreeRetentionReport
from .shared import _write_tsv


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

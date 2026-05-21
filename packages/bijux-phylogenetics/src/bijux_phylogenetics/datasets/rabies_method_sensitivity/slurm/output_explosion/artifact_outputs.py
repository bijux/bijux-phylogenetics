from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path

from .contracts import RabiesMethodSensitivitySlurmOutputExplosionReport
from .shared import _format_float, _write_tsv


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

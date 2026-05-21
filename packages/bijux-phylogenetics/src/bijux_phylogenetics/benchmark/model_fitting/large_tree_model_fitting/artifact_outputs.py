from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import LargeTreeModelFittingBenchmarkReport


def write_large_tree_model_fitting_summary_table(
    path: Path,
    report: LargeTreeModelFittingBenchmarkReport,
) -> Path:
    """Write the stable benchmark summary counts for one governed tier."""
    return write_taxon_rows(
        path,
        columns=[
            "tier",
            "case_count",
            "geiger_match_case_count",
            "threshold_pass_case_count",
            "too_slow_case_count",
            "unstable_case_count",
            "limitations",
        ],
        rows=[
            {
                "tier": report.tier,
                "case_count": report.case_count,
                "geiger_match_case_count": report.geiger_match_case_count,
                "threshold_pass_case_count": report.threshold_pass_case_count,
                "too_slow_case_count": report.too_slow_case_count,
                "unstable_case_count": report.unstable_case_count,
                "limitations": " | ".join(report.limitations),
            }
        ],
    )


def write_large_tree_model_fitting_observation_table(
    path: Path,
    report: LargeTreeModelFittingBenchmarkReport,
) -> Path:
    """Write one stable observation row per governed large-tree benchmark case."""
    return write_taxon_rows(
        path,
        columns=[
            "case_id",
            "tier",
            "trait_kind",
            "fit_surface",
            "taxon_count",
            "status",
            "converged",
            "stable_conclusion_supported",
            "unstable_review",
            "too_slow_review",
            "performance_threshold_passed",
            "matches_geiger_reference",
            "geiger_reference_available",
            "notes",
        ],
        rows=[
            {
                "case_id": row.case_id,
                "tier": row.tier,
                "trait_kind": row.trait_kind,
                "fit_surface": row.fit_surface,
                "taxon_count": row.taxon_count,
                "status": row.status,
                "converged": format_optional_bool(row.converged),
                "stable_conclusion_supported": format_optional_bool(
                    row.stable_conclusion_supported
                ),
                "unstable_review": row.unstable_review,
                "too_slow_review": row.too_slow_review,
                "performance_threshold_passed": format_optional_bool(
                    row.performance_threshold_passed
                ),
                "matches_geiger_reference": format_optional_bool(
                    row.matches_geiger_reference
                ),
                "geiger_reference_available": row.geiger_reference_available,
                "notes": " | ".join(row.notes),
            }
            for row in report.observations
        ],
    )


def format_optional_bool(value: bool | None) -> str:
    if value is None:
        return ""
    return "True" if value else "False"

from __future__ import annotations

from pathlib import Path
from typing import Any

from ....artifacts import preview_report_rows, section
from .shared import artifact_link, preview_payload


def build_instability_sections(
    *,
    unstable_taxa: Any,
    unstable_clades: Any,
    clade_conflicts: Any | None,
    conclusion_summary: Any | None,
    scaled_report_note: dict[str, object],
    artifact_paths: dict[str, Path],
    out_path: Path,
    preview_limit: int,
    unstable_taxa_rows: list[dict[str, Any]],
    unstable_taxa_truncated: int,
    unstable_clade_rows: list[dict[str, Any]],
    unstable_clade_truncated: int,
    conflict_rows: list[dict[str, Any]],
    conflict_truncated: int,
    robust_rows: list[dict[str, Any]],
    robust_truncated: int,
    uncertain_rows: list[dict[str, Any]],
    uncertain_truncated: int,
    conflicting_rows: list[dict[str, Any]],
    conflicting_truncated: int,
) -> list[tuple[str, object]]:
    """Build instability and conclusion sections for the tree uncertainty report."""
    return [
        section(
            "unstable-taxa",
            {
                "tree_count": unstable_taxa.tree_count,
                **preview_payload(
                    rows=unstable_taxa_rows,
                    row_count=len(unstable_taxa.taxa),
                    truncated_row_count=unstable_taxa_truncated,
                    preview_limit=preview_limit,
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="unstable_taxa",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "unstable-clades",
            {
                "tree_count": unstable_clades.tree_count,
                **preview_payload(
                    rows=unstable_clade_rows,
                    row_count=len(unstable_clades.clades),
                    truncated_row_count=unstable_clade_truncated,
                    preview_limit=preview_limit,
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="unstable_clades",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "clade-credibility-conflicts",
            {
                **(
                    {
                        "tree_count": clade_conflicts.tree_count,
                        "credibility_threshold": clade_conflicts.credibility_threshold,
                        "high_credibility_clade_count": clade_conflicts.high_credibility_clade_count,
                        **preview_payload(
                            rows=conflict_rows,
                            row_count=len(clade_conflicts.conflicts),
                            truncated_row_count=conflict_truncated,
                            preview_limit=preview_limit,
                        ),
                    }
                    if clade_conflicts is not None
                    else scaled_report_note
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="clade_credibility_conflicts",
                    out_path=out_path,
                ),
            },
        ),
        section(
            "uncertainty-aware-conclusions",
            {
                **(
                    {
                        "tree_count": conclusion_summary.tree_count,
                        "robust_clade_count": conclusion_summary.robust_clade_count,
                        "uncertain_clade_count": conclusion_summary.uncertain_clade_count,
                        "conflicting_clade_count": conclusion_summary.conflicting_clade_count,
                        "robust_rows": preview_report_rows(
                            robust_rows, limit=preview_limit
                        ),
                        "robust_truncated_row_count": robust_truncated,
                        "uncertain_rows": preview_report_rows(
                            uncertain_rows, limit=preview_limit
                        ),
                        "uncertain_truncated_row_count": uncertain_truncated,
                        "conflicting_rows": preview_report_rows(
                            conflicting_rows, limit=preview_limit
                        ),
                        "conflicting_truncated_row_count": conflicting_truncated,
                    }
                    if conclusion_summary is not None
                    else scaled_report_note
                ),
                "artifact_path": artifact_link(
                    artifact_paths=artifact_paths,
                    artifact_key="uncertainty_aware_conclusions",
                    out_path=out_path,
                ),
            },
        ),
    ]

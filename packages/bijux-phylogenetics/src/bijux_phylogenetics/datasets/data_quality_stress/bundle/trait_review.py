from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    CatarrhineDataQualityStressPanelWorkflowReport,
    TraitDuplicateResolution,
    TraitMissingObservation,
)


def _write_trait_duplicates_table(
    path: Path,
    rows: list[TraitDuplicateResolution],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "occurrence_count",
            "selected_row_number",
            "selected_non_missing_field_count",
            "discarded_row_numbers",
            "selected_reason",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "occurrence_count": str(row.occurrence_count),
                "selected_row_number": str(row.selected_row_number),
                "selected_non_missing_field_count": str(
                    row.selected_non_missing_field_count
                ),
                "discarded_row_numbers": ",".join(
                    str(value) for value in row.discarded_row_numbers
                ),
                "selected_reason": row.selected_reason,
            }
            for row in rows
        ],
    )


def _write_trait_missing_values_table(
    path: Path,
    rows: list[TraitMissingObservation],
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "row_number",
            "trait",
            "required_for_analysis",
            "action",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "row_number": str(row.row_number),
                "trait": row.trait,
                "required_for_analysis": str(row.required_for_analysis).lower(),
                "action": row.action,
            }
            for row in rows
        ],
    )


def _write_raw_trait_linkage_table(
    path: Path,
    report: CatarrhineDataQualityStressPanelWorkflowReport,
) -> Path:
    linkage = report.raw_trait_mismatch_linkage
    return write_taxon_rows(
        path,
        columns=[
            "surface",
            "tree_taxa",
            "trait_taxa",
            "linked_taxa",
            "missing_from_traits",
            "extra_trait_taxa",
            "strict_status",
            "detail",
        ],
        rows=[
            {
                "surface": "raw_trait_mismatch",
                "tree_taxa": str(linkage.tree_taxa),
                "trait_taxa": str(linkage.trait_taxa),
                "linked_taxa": str(linkage.linked_taxa),
                "missing_from_traits": ",".join(linkage.missing_from_traits),
                "extra_trait_taxa": ",".join(linkage.extra_trait_taxa),
                "strict_status": (
                    "failed"
                    if report.raw_trait_mismatch_error is not None
                    else "passed"
                ),
                "detail": report.raw_trait_mismatch_error or "raw linkage passed",
            }
        ],
    )

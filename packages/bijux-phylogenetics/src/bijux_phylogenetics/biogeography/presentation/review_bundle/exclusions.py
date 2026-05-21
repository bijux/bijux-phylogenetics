from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylogeography.geographic_map import GeographicMapReport

from .contracts import BiogeographyReportExclusionRow


def write_biogeography_report_exclusion_table(
    path: Path,
    rows: list[BiogeographyReportExclusionRow],
) -> Path:
    """Write one combined exclusion ledger for the full report package."""
    return write_taxon_rows(
        path,
        columns=[
            "surface",
            "subject_id",
            "subject_kind",
            "raw_left",
            "raw_right",
            "reason",
            "note",
        ],
        rows=[
            {
                "surface": row.surface,
                "subject_id": row.subject_id,
                "subject_kind": row.subject_kind,
                "raw_left": row.raw_left,
                "raw_right": row.raw_right,
                "reason": row.reason,
                "note": row.note,
            }
            for row in rows
        ],
    )


def build_biogeography_report_exclusion_rows(
    state_report: GeographicStateModelReport,
    map_report: GeographicMapReport,
) -> list[BiogeographyReportExclusionRow]:
    rows: list[BiogeographyReportExclusionRow] = [
        BiogeographyReportExclusionRow(
            surface="state_model",
            subject_id=row.taxon,
            subject_kind="taxon",
            raw_left=row.raw_region,
            raw_right=row.normalized_region or "",
            reason=row.reason,
            note=row.note,
        )
        for row in state_report.exclusion_rows
    ]
    rows.extend(
        BiogeographyReportExclusionRow(
            surface="map",
            subject_id=row.subject_id,
            subject_kind=row.subject_kind,
            raw_left=row.raw_left,
            raw_right=row.raw_right,
            reason=row.reason,
            note=row.note,
        )
        for row in map_report.exclusion_rows
    )
    deduplicated: list[BiogeographyReportExclusionRow] = []
    seen: set[tuple[str, str, str, str, str, str, str]] = set()
    for row in rows:
        key = (
            row.surface,
            row.subject_id,
            row.subject_kind,
            row.raw_left,
            row.raw_right,
            row.reason,
            row.note,
        )
        if key in seen:
            continue
        seen.add(key)
        deduplicated.append(row)
    return deduplicated

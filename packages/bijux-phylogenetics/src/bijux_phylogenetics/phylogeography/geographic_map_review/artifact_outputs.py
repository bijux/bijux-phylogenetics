from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import GeographicMapReport
from .shared import format_optional_float


def write_geographic_map_summary_table(
    path: Path,
    report: GeographicMapReport,
) -> Path:
    """Write one overall geographic map summary ledger."""
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "mode",
            "model",
            "analyzed_taxon_count",
            "excluded_record_count",
            "tip_marker_count",
            "internal_marker_count",
            "root_marker_count",
            "line_count",
            "visible_line_count",
            "tree_depth",
            "time_filter_applied",
            "minimum_midpoint_depth_filter",
            "maximum_midpoint_depth_filter",
            "earliest_visible_midpoint_depth",
            "latest_visible_midpoint_depth",
            "warning_count",
        ],
        rows=[
            {
                "mode": summary.mode,
                "model": summary.model,
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_record_count": str(summary.excluded_record_count),
                "tip_marker_count": str(summary.tip_marker_count),
                "internal_marker_count": str(summary.internal_marker_count),
                "root_marker_count": str(summary.root_marker_count),
                "line_count": str(summary.line_count),
                "visible_line_count": str(summary.visible_line_count),
                "tree_depth": str(summary.tree_depth),
                "time_filter_applied": str(summary.time_filter_applied).lower(),
                "minimum_midpoint_depth_filter": format_optional_float(
                    summary.minimum_midpoint_depth_filter
                ),
                "maximum_midpoint_depth_filter": format_optional_float(
                    summary.maximum_midpoint_depth_filter
                ),
                "earliest_visible_midpoint_depth": format_optional_float(
                    summary.earliest_visible_midpoint_depth
                ),
                "latest_visible_midpoint_depth": format_optional_float(
                    summary.latest_visible_midpoint_depth
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_geographic_map_marker_table(
    path: Path,
    report: GeographicMapReport,
) -> Path:
    """Write one geographic marker ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "marker_id",
            "label",
            "marker_kind",
            "latitude",
            "longitude",
            "state_label",
            "descendant_taxa",
            "confidence",
            "uncertainty_km",
            "is_tip",
            "is_root",
            "active_line_count",
        ],
        rows=[
            {
                "marker_id": row.marker_id,
                "label": row.label,
                "marker_kind": row.marker_kind,
                "latitude": str(row.latitude),
                "longitude": str(row.longitude),
                "state_label": row.state_label,
                "descendant_taxa": ",".join(row.descendant_taxa),
                "confidence": str(row.confidence),
                "uncertainty_km": format_optional_float(row.uncertainty_km),
                "is_tip": str(row.is_tip).lower(),
                "is_root": str(row.is_root).lower(),
                "active_line_count": str(row.active_line_count),
            }
            for row in report.marker_rows
        ],
    )


def write_geographic_map_line_table(
    path: Path,
    report: GeographicMapReport,
) -> Path:
    """Write one geographic line ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "line_id",
            "line_kind",
            "source_label",
            "target_label",
            "source_latitude",
            "source_longitude",
            "target_latitude",
            "target_longitude",
            "child_descendant_taxa",
            "support",
            "midpoint_depth",
            "branch_length",
            "distance_km",
            "state_transition",
            "flag_codes",
            "visible",
        ],
        rows=[
            {
                "line_id": row.line_id,
                "line_kind": row.line_kind,
                "source_label": row.source_label,
                "target_label": row.target_label,
                "source_latitude": str(row.source_latitude),
                "source_longitude": str(row.source_longitude),
                "target_latitude": str(row.target_latitude),
                "target_longitude": str(row.target_longitude),
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "support": str(row.support),
                "midpoint_depth": format_optional_float(row.midpoint_depth),
                "branch_length": format_optional_float(row.branch_length),
                "distance_km": str(row.distance_km),
                "state_transition": row.state_transition,
                "flag_codes": ",".join(row.flag_codes),
                "visible": str(row.visible).lower(),
            }
            for row in report.line_rows
        ],
    )


def write_geographic_map_exclusion_table(
    path: Path,
    report: GeographicMapReport,
) -> Path:
    """Write one geographic map exclusion ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "subject_id",
            "subject_kind",
            "raw_left",
            "raw_right",
            "reason",
            "note",
        ],
        rows=[
            {
                "subject_id": row.subject_id,
                "subject_kind": row.subject_kind,
                "raw_left": row.raw_left,
                "raw_right": row.raw_right,
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )

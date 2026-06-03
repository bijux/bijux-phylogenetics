from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.contracts import (
    PhylogeographicCoordinateReport,
)
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.movement_analysis import (
    stable_float,
)


def write_coordinate_movement_summary_table(
    path: Path,
    report: PhylogeographicCoordinateReport,
) -> Path:
    summary = report.summary
    return write_taxon_rows(
        path,
        columns=[
            "taxon_column",
            "latitude_column",
            "longitude_column",
            "model",
            "alpha",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "internal_node_count",
            "weak_support_node_count",
            "outlier_jump_count",
            "impossible_jump_count",
            "flagged_branch_count",
            "maximum_jump_km",
            "root_latitude",
            "root_longitude",
            "root_radial_standard_error_km",
            "warning_count",
        ],
        rows=[
            {
                "taxon_column": summary.taxon_column,
                "latitude_column": summary.latitude_column,
                "longitude_column": summary.longitude_column,
                "model": summary.model,
                "alpha": str(summary.alpha),
                "analyzed_taxon_count": str(summary.analyzed_taxon_count),
                "excluded_taxon_count": str(summary.excluded_taxon_count),
                "internal_node_count": str(summary.internal_node_count),
                "weak_support_node_count": str(summary.weak_support_node_count),
                "outlier_jump_count": str(summary.outlier_jump_count),
                "impossible_jump_count": str(summary.impossible_jump_count),
                "flagged_branch_count": str(summary.flagged_branch_count),
                "maximum_jump_km": str(summary.maximum_jump_km),
                "root_latitude": str(summary.root_latitude),
                "root_longitude": str(summary.root_longitude),
                "root_radial_standard_error_km": str(
                    summary.root_radial_standard_error_km
                ),
                "warning_count": str(summary.warning_count),
            }
        ],
    )


def write_coordinate_estimate_table(
    path: Path,
    report: PhylogeographicCoordinateReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "node",
            "node_name",
            "is_tip",
            "descendant_taxa",
            "latitude",
            "longitude",
            "latitude_standard_error",
            "longitude_standard_error",
            "radial_standard_error_km",
            "lower_95_latitude",
            "upper_95_latitude",
            "lower_95_longitude",
            "upper_95_longitude",
            "confidence",
            "unstable",
            "is_root",
        ],
        rows=[
            {
                "node": row.node,
                "node_name": row.node_name or "",
                "is_tip": str(row.is_tip).lower(),
                "descendant_taxa": ",".join(row.descendant_taxa),
                "latitude": str(row.latitude),
                "longitude": str(row.longitude),
                "latitude_standard_error": str(row.latitude_standard_error),
                "longitude_standard_error": str(row.longitude_standard_error),
                "radial_standard_error_km": str(row.radial_standard_error_km),
                "lower_95_latitude": str(row.lower_95_latitude),
                "upper_95_latitude": str(row.upper_95_latitude),
                "lower_95_longitude": str(row.lower_95_longitude),
                "upper_95_longitude": str(row.upper_95_longitude),
                "confidence": str(row.confidence),
                "unstable": str(row.unstable).lower(),
                "is_root": str(row.is_root).lower(),
            }
            for row in report.estimate_rows
        ],
    )


def write_coordinate_movement_branch_table(
    path: Path,
    report: PhylogeographicCoordinateReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "branch_length",
            "parent_latitude",
            "parent_longitude",
            "child_latitude",
            "child_longitude",
            "great_circle_km",
            "branch_rate_km_per_unit",
            "support",
            "impossible_jump",
            "outlier_jump",
            "flag_codes",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "branch_length": _format_optional_float(row.branch_length),
                "parent_latitude": str(row.parent_latitude),
                "parent_longitude": str(row.parent_longitude),
                "child_latitude": str(row.child_latitude),
                "child_longitude": str(row.child_longitude),
                "great_circle_km": str(row.great_circle_km),
                "branch_rate_km_per_unit": _format_optional_float(
                    row.branch_rate_km_per_unit
                ),
                "support": str(row.support),
                "impossible_jump": str(row.impossible_jump).lower(),
                "outlier_jump": str(row.outlier_jump).lower(),
                "flag_codes": ",".join(row.flag_codes),
            }
            for row in report.branch_rows
        ],
    )


def write_coordinate_movement_outlier_table(
    path: Path,
    report: PhylogeographicCoordinateReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "parent_node",
            "child_node",
            "child_descendant_taxa",
            "great_circle_km",
            "branch_rate_km_per_unit",
            "median_distance_km",
            "distance_threshold_km",
            "median_rate_km_per_unit",
            "rate_threshold_km",
            "impossible_jump",
            "outlier_jump",
            "flag_codes",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "parent_node": row.parent_node,
                "child_node": row.child_node,
                "child_descendant_taxa": ",".join(row.child_descendant_taxa),
                "great_circle_km": str(row.great_circle_km),
                "branch_rate_km_per_unit": _format_optional_float(
                    row.branch_rate_km_per_unit
                ),
                "median_distance_km": str(row.median_distance_km),
                "distance_threshold_km": str(row.distance_threshold_km),
                "median_rate_km_per_unit": _format_optional_float(
                    row.median_rate_km_per_unit
                ),
                "rate_threshold_km": _format_optional_float(row.rate_threshold_km),
                "impossible_jump": str(row.impossible_jump).lower(),
                "outlier_jump": str(row.outlier_jump).lower(),
                "flag_codes": ",".join(row.flag_codes),
            }
            for row in report.outlier_rows
        ],
    )


def write_coordinate_movement_exclusion_table(
    path: Path,
    report: PhylogeographicCoordinateReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "taxon",
            "raw_latitude",
            "raw_longitude",
            "reason",
            "note",
        ],
        rows=[
            {
                "taxon": row.taxon,
                "raw_latitude": row.raw_latitude,
                "raw_longitude": row.raw_longitude,
                "reason": row.reason,
                "note": row.note,
            }
            for row in report.exclusion_rows
        ],
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(stable_float(value))

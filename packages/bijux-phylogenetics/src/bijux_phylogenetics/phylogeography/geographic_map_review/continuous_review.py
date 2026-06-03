from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.phylogeography.continuous_coordinates import (
    summarize_continuous_phylogeography,
)

from .contracts import (
    GeographicMapExclusionRow,
    GeographicMapLineRow,
    GeographicMapMarkerRow,
    GeographicMapReport,
)
from .shared import (
    active_line_counts,
    apply_depth_filter,
    branch_midpoint_depth,
    build_map_summary,
    continuous_depths,
    validate_depth_filters,
)


def summarize_continuous_phylogeography_map(
    tree_path: Path,
    table_path: Path,
    *,
    latitude_column: str,
    longitude_column: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
    minimum_midpoint_depth: float | None = None,
    maximum_midpoint_depth: float | None = None,
) -> GeographicMapReport:
    """Build one map review surface from continuous ancestral coordinates."""
    validate_depth_filters(
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    coordinate_report = summarize_continuous_phylogeography(
        tree_path,
        table_path,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        taxon_column=taxon_column,
        model=model,
        alpha=alpha,
    )
    depth_by_node, tree_depth = continuous_depths(coordinate_report)
    line_rows = [
        GeographicMapLineRow(
            line_id=row.branch_id,
            line_kind="continuous_movement",
            source_label=row.parent_node,
            target_label=row.child_node,
            source_latitude=row.parent_latitude,
            source_longitude=row.parent_longitude,
            target_latitude=row.child_latitude,
            target_longitude=row.child_longitude,
            child_descendant_taxa=list(row.child_descendant_taxa),
            support=row.support,
            midpoint_depth=branch_midpoint_depth(
                depth_by_node.get(row.parent_node),
                depth_by_node.get(row.child_node),
            ),
            branch_length=row.branch_length,
            distance_km=row.great_circle_km,
            state_transition="",
            flag_codes=list(row.flag_codes),
            visible=False,
        )
        for row in coordinate_report.branch_rows
    ]
    line_rows = apply_depth_filter(
        line_rows,
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    visible_line_counts = active_line_counts(line_rows)
    marker_rows = [
        GeographicMapMarkerRow(
            marker_id=row.node,
            label=row.node_name or row.node,
            marker_kind=(
                "root" if row.is_root else "tip" if row.is_tip else "ancestral"
            ),
            latitude=row.latitude,
            longitude=row.longitude,
            state_label="",
            descendant_taxa=list(row.descendant_taxa),
            confidence=row.confidence,
            uncertainty_km=row.radial_standard_error_km,
            is_tip=row.is_tip,
            is_root=row.is_root,
            active_line_count=visible_line_counts.get(row.node, 0),
        )
        for row in coordinate_report.estimate_rows
    ]
    exclusion_rows = [
        GeographicMapExclusionRow(
            subject_id=row.taxon,
            subject_kind="taxon",
            raw_left=row.raw_latitude,
            raw_right=row.raw_longitude,
            reason=row.reason,
            note=row.note,
        )
        for row in coordinate_report.exclusion_rows
    ]
    warnings = list(coordinate_report.warnings)
    if minimum_midpoint_depth is not None or maximum_midpoint_depth is not None:
        warnings.append(
            "map output applies a midpoint-depth line filter; branch lines outside the selected depth window remain excluded from the visible map layer"
        )
    if not any(row.visible for row in line_rows):
        warnings.append(
            "no geographic movement lines remain visible after depth filtering"
        )
    summary = build_map_summary(
        mode="continuous",
        model=coordinate_report.model,
        analyzed_taxon_count=coordinate_report.summary.analyzed_taxon_count,
        marker_rows=marker_rows,
        line_rows=line_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
        tree_depth=tree_depth,
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    return GeographicMapReport(
        mode="continuous",
        tree_path=tree_path,
        table_path=table_path,
        summary=summary,
        marker_rows=marker_rows,
        line_rows=line_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )

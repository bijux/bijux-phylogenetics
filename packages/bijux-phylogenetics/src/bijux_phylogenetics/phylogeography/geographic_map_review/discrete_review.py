from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.biogeography import (
    summarize_geographic_migration_events,
    summarize_geographic_state_model,
)

from .contracts import (
    GeographicMapExclusionRow,
    GeographicMapMarkerRow,
    GeographicMapReport,
)
from .shared import (
    active_line_counts,
    build_discrete_line_rows,
    build_discrete_node_markers,
    build_discrete_tip_markers,
    build_map_summary,
    load_region_centroids,
    validate_depth_filters,
)


def summarize_discrete_region_map(
    tree_path: Path,
    table_path: Path,
    *,
    trait: str,
    centroids_path: Path,
    taxon_column: str | None = None,
    model: str = "er",
    region_column: str = "region",
    latitude_column: str = "latitude",
    longitude_column: str = "longitude",
    minimum_midpoint_depth: float | None = None,
    maximum_midpoint_depth: float | None = None,
) -> GeographicMapReport:
    """Build one map review surface from discrete geographic-region reconstruction."""
    validate_depth_filters(
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    state_report = summarize_geographic_state_model(
        tree_path,
        table_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    event_report = summarize_geographic_migration_events(
        tree_path,
        table_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    centroids = load_region_centroids(
        centroids_path,
        region_column=region_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
    )
    tip_markers, centroid_exclusions = build_discrete_tip_markers(
        tree_path=tree_path,
        table_path=table_path,
        taxon_column=state_report.taxon_column,
        trait=trait,
        centroids=centroids,
        excluded_taxa={row.taxon for row in state_report.exclusion_rows},
    )
    node_markers, node_exclusions = build_discrete_node_markers(
        state_report,
        centroids=centroids,
    )
    line_rows, line_exclusions = build_discrete_line_rows(
        event_report,
        centroids=centroids,
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    all_markers = [*tip_markers, *node_markers]
    visible_line_counts = active_line_counts(line_rows)
    marker_rows = [
        GeographicMapMarkerRow(
            marker_id=row.marker_id,
            label=row.label,
            marker_kind=row.marker_kind,
            latitude=row.latitude,
            longitude=row.longitude,
            state_label=row.state_label,
            descendant_taxa=list(row.descendant_taxa),
            confidence=row.confidence,
            uncertainty_km=row.uncertainty_km,
            is_tip=row.is_tip,
            is_root=row.is_root,
            active_line_count=visible_line_counts.get(row.marker_id, 0),
        )
        for row in all_markers
    ]
    exclusion_rows = [
        GeographicMapExclusionRow(
            subject_id=row.taxon,
            subject_kind="taxon",
            raw_left=row.raw_region,
            raw_right=row.normalized_region or "",
            reason=row.reason,
            note=row.note,
        )
        for row in state_report.exclusion_rows
    ]
    exclusion_rows.extend(centroid_exclusions)
    exclusion_rows.extend(node_exclusions)
    exclusion_rows.extend(line_exclusions)
    warnings = list(dict.fromkeys([*state_report.warnings, *event_report.warnings]))
    if minimum_midpoint_depth is not None or maximum_midpoint_depth is not None:
        warnings.append(
            "map output applies a midpoint-depth event filter; geographic transition lines outside the selected depth window remain excluded from the visible map layer"
        )
    if not any(row.visible for row in line_rows):
        warnings.append(
            "no geographic transition lines remain visible after depth filtering"
        )
    summary = build_map_summary(
        mode="regions",
        model=state_report.model,
        analyzed_taxon_count=state_report.summary.analyzed_taxon_count,
        marker_rows=marker_rows,
        line_rows=line_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
        tree_depth=event_report.summary.tree_depth,
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    return GeographicMapReport(
        mode="regions",
        tree_path=tree_path,
        table_path=table_path,
        summary=summary,
        marker_rows=marker_rows,
        line_rows=line_rows,
        exclusion_rows=exclusion_rows,
        warnings=warnings,
    )

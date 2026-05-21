from __future__ import annotations

from html import escape
import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature, stable_value
from bijux_phylogenetics.biogeography import (
    summarize_geographic_migration_events,
    summarize_geographic_state_model,
)
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.phylogeography.continuous_coordinates import (
    PhylogeographicCoordinateReport,
    summarize_continuous_phylogeography,
)
from bijux_phylogenetics.phylogeography.region_styles import (
    build_geographic_state_color_map,
    geographic_transition_support_colors,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError
from bijux_phylogenetics.phylogeography.geographic_map_review.contracts import (
    GeographicMapArtifact,
    GeographicMapExclusionRow,
    GeographicMapLineRow,
    GeographicMapMarkerRow,
    GeographicMapReport,
    GeographicMapSummary,
)

_EARTH_RADIUS_KM = 6371.0088


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
    _validate_depth_filters(
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
    depth_by_node, tree_depth = _continuous_depths(coordinate_report)
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
            midpoint_depth=_branch_midpoint_depth(
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
    line_rows = _apply_depth_filter(
        line_rows,
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    active_line_counts = _active_line_counts(line_rows)
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
            active_line_count=active_line_counts.get(row.node, 0),
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
    summary = _build_map_summary(
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
    _validate_depth_filters(
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
    centroids = _load_region_centroids(
        centroids_path,
        region_column=region_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
    )
    tip_markers, centroid_exclusions = _build_discrete_tip_markers(
        tree_path=tree_path,
        table_path=table_path,
        taxon_column=state_report.taxon_column,
        trait=trait,
        centroids=centroids,
        excluded_taxa={row.taxon for row in state_report.exclusion_rows},
    )
    node_markers, node_exclusions = _build_discrete_node_markers(
        state_report,
        centroids=centroids,
    )
    line_rows, line_exclusions = _build_discrete_line_rows(
        event_report,
        centroids=centroids,
        minimum_midpoint_depth=minimum_midpoint_depth,
        maximum_midpoint_depth=maximum_midpoint_depth,
    )
    all_markers = [*tip_markers, *node_markers]
    active_line_counts = _active_line_counts(line_rows)
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
            active_line_count=active_line_counts.get(row.marker_id, 0),
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
    summary = _build_map_summary(
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


def render_geographic_map_html(
    report: GeographicMapReport,
    *,
    out_path: Path,
    width: int = 1200,
    height: int = 720,
    state_colors: dict[str, str] | None = None,
) -> GeographicMapArtifact:
    """Render one self-contained HTML geographic map review artifact."""
    if out_path.suffix.lower() != ".html":
        raise ValueError("geographic map output must end in .html")
    html = _build_map_html(
        report,
        width=width,
        height=height,
        state_colors=state_colors,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(html, encoding="utf-8")
    return GeographicMapArtifact(
        output_path=out_path,
        format="html",
        width=width,
        height=height,
        visible_line_count=sum(row.visible for row in report.line_rows),
    )


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
                "minimum_midpoint_depth_filter": _format_optional_float(
                    summary.minimum_midpoint_depth_filter
                ),
                "maximum_midpoint_depth_filter": _format_optional_float(
                    summary.maximum_midpoint_depth_filter
                ),
                "earliest_visible_midpoint_depth": _format_optional_float(
                    summary.earliest_visible_midpoint_depth
                ),
                "latest_visible_midpoint_depth": _format_optional_float(
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
                "uncertainty_km": _format_optional_float(row.uncertainty_km),
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
                "midpoint_depth": _format_optional_float(row.midpoint_depth),
                "branch_length": _format_optional_float(row.branch_length),
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


def _validate_depth_filters(
    *,
    minimum_midpoint_depth: float | None,
    maximum_midpoint_depth: float | None,
) -> None:
    if minimum_midpoint_depth is not None and minimum_midpoint_depth < 0.0:
        raise ValueError("minimum midpoint depth filter must be non-negative")
    if maximum_midpoint_depth is not None and maximum_midpoint_depth < 0.0:
        raise ValueError("maximum midpoint depth filter must be non-negative")
    if (
        minimum_midpoint_depth is not None
        and maximum_midpoint_depth is not None
        and minimum_midpoint_depth > maximum_midpoint_depth
    ):
        raise ValueError(
            "minimum midpoint depth filter cannot exceed maximum midpoint depth filter"
        )


def _continuous_depths(
    report: PhylogeographicCoordinateReport,
) -> tuple[dict[str, float], float]:
    analysis_taxa = sorted(
        row.node_name
        for row in report.estimate_rows
        if row.is_tip and row.node_name is not None
    )
    tree, _ = prune_tree_to_requested_taxa(report.tree_path, analysis_taxa)
    depth_by_node: dict[str, float] = {}

    def walk(node: TreeNode, depth: float) -> None:
        depth_by_node[node_signature(node)] = stable_value(depth)
        for child in node.children:
            child_depth = depth + (child.branch_length or 0.0)
            walk(child, child_depth)

    walk(tree.root, 0.0)
    tree_depth = max(depth_by_node.values(), default=0.0)
    return depth_by_node, stable_value(tree_depth)


def _branch_midpoint_depth(
    parent_depth: float | None,
    child_depth: float | None,
) -> float | None:
    if parent_depth is None or child_depth is None:
        return None
    return stable_value((parent_depth + child_depth) / 2.0)


def _apply_depth_filter(
    rows: list[GeographicMapLineRow],
    *,
    minimum_midpoint_depth: float | None,
    maximum_midpoint_depth: float | None,
) -> list[GeographicMapLineRow]:
    reviewed: list[GeographicMapLineRow] = []
    for row in rows:
        midpoint_depth = row.midpoint_depth
        visible = True
        if midpoint_depth is not None and minimum_midpoint_depth is not None:
            visible = visible and midpoint_depth >= minimum_midpoint_depth
        if midpoint_depth is not None and maximum_midpoint_depth is not None:
            visible = visible and midpoint_depth <= maximum_midpoint_depth
        reviewed.append(
            GeographicMapLineRow(
                line_id=row.line_id,
                line_kind=row.line_kind,
                source_label=row.source_label,
                target_label=row.target_label,
                source_latitude=row.source_latitude,
                source_longitude=row.source_longitude,
                target_latitude=row.target_latitude,
                target_longitude=row.target_longitude,
                child_descendant_taxa=list(row.child_descendant_taxa),
                support=row.support,
                midpoint_depth=row.midpoint_depth,
                branch_length=row.branch_length,
                distance_km=row.distance_km,
                state_transition=row.state_transition,
                flag_codes=list(row.flag_codes),
                visible=visible,
            )
        )
    return reviewed


def _active_line_counts(rows: list[GeographicMapLineRow]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        if not row.visible:
            continue
        counts[row.source_label] = counts.get(row.source_label, 0) + 1
        counts[row.target_label] = counts.get(row.target_label, 0) + 1
    return counts


def _build_discrete_tip_markers(
    *,
    tree_path: Path,
    table_path: Path,
    taxon_column: str,
    trait: str,
    centroids: dict[str, tuple[float, float]],
    excluded_taxa: set[str],
) -> tuple[list[GeographicMapMarkerRow], list[GeographicMapExclusionRow]]:
    tree = load_tree(tree_path)
    table = load_taxon_table(table_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    markers: list[GeographicMapMarkerRow] = []
    exclusions: list[GeographicMapExclusionRow] = []
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None or taxon in excluded_taxa:
            continue
        region = row.get(trait, "").strip()
        centroid = centroids.get(region)
        if centroid is None:
            exclusions.append(
                GeographicMapExclusionRow(
                    subject_id=taxon,
                    subject_kind="taxon",
                    raw_left=region,
                    raw_right="",
                    reason="missing-region-centroid",
                    note="taxon region does not have a matching centroid row for map rendering",
                )
            )
            continue
        markers.append(
            GeographicMapMarkerRow(
                marker_id=taxon,
                label=taxon,
                marker_kind="tip",
                latitude=centroid[0],
                longitude=centroid[1],
                state_label=region,
                descendant_taxa=[taxon],
                confidence=1.0,
                uncertainty_km=None,
                is_tip=True,
                is_root=False,
                active_line_count=0,
            )
        )
    return markers, exclusions


def _build_discrete_node_markers(
    report,
    *,
    centroids: dict[str, tuple[float, float]],
) -> tuple[list[GeographicMapMarkerRow], list[GeographicMapExclusionRow]]:
    markers: list[GeographicMapMarkerRow] = []
    exclusions: list[GeographicMapExclusionRow] = []
    for row in report.node_rows:
        centroid = centroids.get(row.most_likely_region)
        if centroid is None:
            exclusions.append(
                GeographicMapExclusionRow(
                    subject_id=row.node,
                    subject_kind="node",
                    raw_left=row.most_likely_region,
                    raw_right="",
                    reason="missing-region-centroid",
                    note="ancestral region does not have a matching centroid row for map rendering",
                )
            )
            continue
        markers.append(
            GeographicMapMarkerRow(
                marker_id=row.node,
                label=row.node_name or row.node,
                marker_kind="root" if row.is_root else "ancestral",
                latitude=centroid[0],
                longitude=centroid[1],
                state_label=row.most_likely_region,
                descendant_taxa=list(row.descendant_taxa),
                confidence=row.confidence,
                uncertainty_km=None,
                is_tip=False,
                is_root=row.is_root,
                active_line_count=0,
            )
        )
    return markers, exclusions


def _build_discrete_line_rows(
    report,
    *,
    centroids: dict[str, tuple[float, float]],
    minimum_midpoint_depth: float | None,
    maximum_midpoint_depth: float | None,
) -> tuple[list[GeographicMapLineRow], list[GeographicMapExclusionRow]]:
    rows: list[GeographicMapLineRow] = []
    exclusions: list[GeographicMapExclusionRow] = []
    for row in report.event_rows:
        source_centroid = centroids.get(row.source_region)
        target_centroid = centroids.get(row.target_region)
        if source_centroid is None or target_centroid is None:
            exclusions.append(
                GeographicMapExclusionRow(
                    subject_id=row.branch_id,
                    subject_kind="event",
                    raw_left=row.source_region,
                    raw_right=row.target_region,
                    reason="missing-region-centroid",
                    note="source or target region does not have a matching centroid row for map rendering",
                )
            )
            continue
        rows.append(
            GeographicMapLineRow(
                line_id=row.branch_id,
                line_kind="region_transition",
                source_label=row.parent_node,
                target_label=row.child_node,
                source_latitude=source_centroid[0],
                source_longitude=source_centroid[1],
                target_latitude=target_centroid[0],
                target_longitude=target_centroid[1],
                child_descendant_taxa=list(row.child_descendant_taxa),
                support=row.support,
                midpoint_depth=row.midpoint_depth,
                branch_length=row.branch_length,
                distance_km=_great_circle_km(
                    source_centroid[0],
                    source_centroid[1],
                    target_centroid[0],
                    target_centroid[1],
                ),
                state_transition=f"{row.source_region}->{row.target_region}",
                flag_codes=[row.confidence_class],
                visible=False,
            )
        )
    return (
        _apply_depth_filter(
            rows,
            minimum_midpoint_depth=minimum_midpoint_depth,
            maximum_midpoint_depth=maximum_midpoint_depth,
        ),
        exclusions,
    )


def _load_region_centroids(
    centroids_path: Path,
    *,
    region_column: str,
    latitude_column: str,
    longitude_column: str,
) -> dict[str, tuple[float, float]]:
    table = load_taxon_table(centroids_path, taxon_column=region_column)
    if latitude_column not in table.columns:
        raise AncestralReconstructionError(
            f"centroid table does not contain column '{latitude_column}'"
        )
    if longitude_column not in table.columns:
        raise AncestralReconstructionError(
            f"centroid table does not contain column '{longitude_column}'"
        )
    centroids: dict[str, tuple[float, float]] = {}
    for row in table.rows:
        region = row[table.taxon_column].strip()
        if not region:
            raise AncestralReconstructionError(
                "centroid table contains an empty region identifier"
            )
        try:
            latitude = float(row[latitude_column].strip())
            longitude = float(row[longitude_column].strip())
        except ValueError as exc:
            raise AncestralReconstructionError(
                f"centroid row for region '{region}' has non-numeric coordinates"
            ) from exc
        if latitude < -90.0 or latitude > 90.0:
            raise AncestralReconstructionError(
                f"centroid row for region '{region}' has latitude outside [-90, 90]"
            )
        if longitude < -180.0 or longitude > 180.0:
            raise AncestralReconstructionError(
                f"centroid row for region '{region}' has longitude outside [-180, 180]"
            )
        centroids[region] = (stable_value(latitude), stable_value(longitude))
    return centroids


def _build_map_summary(
    *,
    mode: str,
    model: str,
    analyzed_taxon_count: int,
    marker_rows: list[GeographicMapMarkerRow],
    line_rows: list[GeographicMapLineRow],
    exclusion_rows: list[GeographicMapExclusionRow],
    warnings: list[str],
    tree_depth: float,
    minimum_midpoint_depth: float | None,
    maximum_midpoint_depth: float | None,
) -> GeographicMapSummary:
    visible_depths = [
        row.midpoint_depth
        for row in line_rows
        if row.visible and row.midpoint_depth is not None
    ]
    return GeographicMapSummary(
        mode=mode,
        model=model,
        analyzed_taxon_count=analyzed_taxon_count,
        excluded_record_count=len(exclusion_rows),
        tip_marker_count=sum(row.is_tip for row in marker_rows),
        internal_marker_count=sum(
            not row.is_tip and not row.is_root for row in marker_rows
        ),
        root_marker_count=sum(row.is_root for row in marker_rows),
        line_count=len(line_rows),
        visible_line_count=sum(row.visible for row in line_rows),
        tree_depth=stable_value(tree_depth),
        time_filter_applied=(
            minimum_midpoint_depth is not None or maximum_midpoint_depth is not None
        ),
        minimum_midpoint_depth_filter=minimum_midpoint_depth,
        maximum_midpoint_depth_filter=maximum_midpoint_depth,
        earliest_visible_midpoint_depth=(
            stable_value(min(visible_depths)) if visible_depths else None
        ),
        latest_visible_midpoint_depth=(
            stable_value(max(visible_depths)) if visible_depths else None
        ),
        warning_count=len(warnings),
    )


def _great_circle_km(
    left_latitude: float,
    left_longitude: float,
    right_latitude: float,
    right_longitude: float,
) -> float:
    left_phi = math.radians(left_latitude)
    right_phi = math.radians(right_latitude)
    delta_phi = math.radians(right_latitude - left_latitude)
    delta_lambda = math.radians(right_longitude - left_longitude)
    a = (
        math.sin(delta_phi / 2.0) ** 2
        + math.cos(left_phi) * math.cos(right_phi) * math.sin(delta_lambda / 2.0) ** 2
    )
    return stable_value(2.0 * _EARTH_RADIUS_KM * math.asin(math.sqrt(a)))


def _build_map_html(
    report: GeographicMapReport,
    *,
    width: int,
    height: int,
    state_colors: dict[str, str] | None = None,
) -> str:
    summary = report.summary
    state_colors = state_colors or build_geographic_state_color_map(
        row.state_label for row in report.marker_rows
    )
    transition_colors = geographic_transition_support_colors()
    svg = _build_map_svg(
        report,
        width=width,
        height=height,
        state_colors=state_colors,
        transition_colors=transition_colors,
    )
    warnings_markup = "".join(
        f"<li>{escape(warning)}</li>" for warning in report.warnings
    )
    title = (
        "Continuous Geographic Map Review"
        if report.mode == "continuous"
        else "Regional Transition Map Review"
    )
    filter_note = (
        f"{_format_optional_float(summary.minimum_midpoint_depth_filter)} to "
        f"{_format_optional_float(summary.maximum_midpoint_depth_filter)}"
        if summary.time_filter_applied
        else "none"
    )
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            f"  <title>{escape(title)}</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f4f7f3 0%, #e4ece8 100%); color: #163127; }",
            "    main { max-width: 1320px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 32px; }",
            "    p { margin: 0 0 12px; line-height: 1.5; }",
            "    .grid { display: grid; grid-template-columns: 340px 1fr; gap: 20px; align-items: start; }",
            "    .panel { background: rgba(255,255,255,0.78); border: 1px solid rgba(22,49,39,0.18); border-radius: 18px; padding: 18px; box-shadow: 0 18px 44px rgba(22,49,39,0.08); }",
            "    .metric { display: flex; justify-content: space-between; gap: 12px; padding: 6px 0; border-bottom: 1px solid rgba(22,49,39,0.08); }",
            "    .metric:last-child { border-bottom: none; }",
            "    .label { color: #496559; }",
            "    ul { margin: 8px 0 0 18px; padding: 0; }",
            "    .legend-block { margin-top: 16px; }",
            "    .legend-grid { display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 10px; margin-top: 8px; }",
            "    .legend-row { display: flex; align-items: center; gap: 8px; font-size: 13px; }",
            "    .legend-swatch { width: 14px; height: 14px; border-radius: 999px; border: 1px solid rgba(22,49,39,0.22); flex: 0 0 auto; }",
            "    .legend-line { width: 24px; height: 0; border-top-width: 3px; border-top-style: solid; flex: 0 0 auto; }",
            "    .map-shell { overflow: auto; }",
            "    svg { width: 100%; height: auto; display: block; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            f"<h1>{escape(title)}</h1>",
            "<p>Self-contained HTML map review over owned phylogenetic geography outputs. The world view uses a fixed latitude/longitude extent so spatial interpretation remains map-based rather than coordinate-space-only.</p>",
            '<div class="grid">',
            '<section class="panel">',
            '<div class="metric"><span class="label">mode</span><strong>'
            + escape(summary.mode)
            + "</strong></div>",
            '<div class="metric"><span class="label">model</span><strong>'
            + escape(summary.model)
            + "</strong></div>",
            '<div class="metric"><span class="label">tip markers</span><strong>'
            + str(summary.tip_marker_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">ancestral markers</span><strong>'
            + str(summary.internal_marker_count + summary.root_marker_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">visible lines</span><strong>'
            + str(summary.visible_line_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">all lines</span><strong>'
            + str(summary.line_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">tree depth</span><strong>'
            + str(summary.tree_depth)
            + "</strong></div>",
            '<div class="metric"><span class="label">midpoint-depth filter</span><strong>'
            + escape(filter_note)
            + "</strong></div>",
            '<div class="metric"><span class="label">excluded records</span><strong>'
            + str(summary.excluded_record_count)
            + "</strong></div>",
            '<div class="metric"><span class="label">warnings</span><strong>'
            + str(summary.warning_count)
            + "</strong></div>",
            "<p><strong>Legend.</strong> Geographic region colors are stable across the tree and map. Marker size still distinguishes tips, ancestral nodes, and the root. Transition line color follows support class, while pale dashed lines remain outside the active depth window.</p>",
            "<p><strong>Time filtering.</strong> When branch lengths represent dated depth, midpoint-depth filters can restrict the visible transition layer without refitting the reconstruction.</p>",
            '<section class="legend-block">',
            "  <strong>Region colors</strong>",
            _legend_grid(
                [
                    _legend_row_html(
                        swatch_color=color,
                        label=label,
                    )
                    for label, color in state_colors.items()
                ]
            ),
            "</section>",
            '<section class="legend-block">',
            "  <strong>Transition support</strong>",
            _legend_grid(
                [
                    _legend_row_html(
                        line_color=color,
                        label=f"{label} support",
                    )
                    for label, color in transition_colors.items()
                ]
            ),
            "</section>",
            "<ul>" + warnings_markup + "</ul>",
            "</section>",
            '<section class="panel map-shell">',
            svg,
            "</section>",
            "</div>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def _build_map_svg(
    report: GeographicMapReport,
    *,
    width: int,
    height: int,
    state_colors: dict[str, str],
    transition_colors: dict[str, str],
) -> str:
    margin_left = 70.0
    margin_right = 28.0
    margin_top = 28.0
    margin_bottom = 52.0
    inner_width = width - margin_left - margin_right
    inner_height = height - margin_top - margin_bottom

    def project(longitude: float, latitude: float) -> tuple[float, float]:
        x = margin_left + ((longitude + 180.0) / 360.0) * inner_width
        y = margin_top + ((90.0 - latitude) / 180.0) * inner_height
        return stable_value(x), stable_value(y)

    grid_lines: list[str] = []
    for longitude in range(-180, 181, 30):
        x, _ = project(float(longitude), 0.0)
        grid_lines.append(
            f'<line x1="{x}" y1="{margin_top}" x2="{x}" y2="{margin_top + inner_height}" stroke="#b9cfc6" stroke-width="1" stroke-dasharray="4 6" />'
        )
        grid_lines.append(
            f'<text x="{x}" y="{height - 16}" font-size="11" text-anchor="middle" fill="#557567">{longitude}</text>'
        )
    for latitude in range(-60, 91, 30):
        _, y = project(0.0, float(latitude))
        grid_lines.append(
            f'<line x1="{margin_left}" y1="{y}" x2="{margin_left + inner_width}" y2="{y}" stroke="#b9cfc6" stroke-width="1" stroke-dasharray="4 6" />'
        )
        grid_lines.append(
            f'<text x="18" y="{y + 4}" font-size="11" fill="#557567">{latitude}</text>'
        )

    line_elements: list[str] = []
    for row in report.line_rows:
        left_x, left_y = project(row.source_longitude, row.source_latitude)
        right_x, right_y = project(row.target_longitude, row.target_latitude)
        confidence_class = row.flag_codes[0] if row.flag_codes else ""
        stroke = (
            "#0c5a74"
            if row.line_kind == "continuous_movement"
            else transition_colors.get(confidence_class, "#186b3d")
        )
        opacity = "0.88" if row.visible else "0.22"
        dasharray = "" if row.visible else "6 7"
        line_elements.append(
            "\n".join(
                [
                    f'<line x1="{left_x}" y1="{left_y}" x2="{right_x}" y2="{right_y}" stroke="{stroke}" stroke-width="{3.6 if row.visible else 2.0}" opacity="{opacity}"'
                    + (f' stroke-dasharray="{dasharray}"' if dasharray else "")
                    + " >",
                    "<title>"
                    + escape(
                        f"{row.source_label} -> {row.target_label}; distance={row.distance_km} km; midpoint_depth={_format_optional_float(row.midpoint_depth)}; support={row.support}; transition={row.state_transition or 'continuous'}"
                    )
                    + "</title>",
                    "</line>",
                ]
            )
        )

    marker_elements: list[str] = []
    for row in report.marker_rows:
        x, y = project(row.longitude, row.latitude)
        fill = state_colors.get(row.state_label, "#0f6090")
        radius = 5.0 if row.is_tip else 7.0
        stroke = "#ffffff"
        stroke_width = 1.4
        if row.is_root:
            radius = 8.0
            stroke = "#17261f"
            stroke_width = 2.1
        opacity = 0.95 if row.active_line_count or row.is_tip else 0.78
        marker_elements.append(
            "\n".join(
                [
                    f'<circle cx="{x}" cy="{y}" r="{radius}" fill="{fill}" stroke="{stroke}" stroke-width="{stroke_width}" opacity="{opacity}">',
                    "<title>"
                    + escape(
                        f"{row.label}; lat={row.latitude}; lon={row.longitude}; state={row.state_label or 'continuous'}; confidence={row.confidence}; uncertainty_km={_format_optional_float(row.uncertainty_km)}"
                    )
                    + "</title>",
                    "</circle>",
                ]
            )
        )
    return "\n".join(
        [
            f'<svg viewBox="0 0 {width} {height}" role="img" aria-label="Geographic reconstruction map">',
            f'<rect x="{margin_left}" y="{margin_top}" width="{inner_width}" height="{inner_height}" rx="24" fill="#edf5f2" stroke="#7aa28f" stroke-width="1.5" />',
            *grid_lines,
            *line_elements,
            *marker_elements,
            f'<text x="{margin_left}" y="18" font-size="13" fill="#355648">longitude / latitude map extent</text>',
            "</svg>",
        ]
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(stable_value(value))


def _legend_row_html(
    *,
    label: str,
    swatch_color: str | None = None,
    line_color: str | None = None,
) -> str:
    icon = ""
    if swatch_color is not None:
        icon = f'<span class="legend-swatch" style="background:{escape(swatch_color)};"></span>'
    elif line_color is not None:
        icon = f'<span class="legend-line" style="border-top-color:{escape(line_color)};"></span>'
    return f'<div class="legend-row">{icon}<span>{escape(label)}</span></div>'


def _legend_grid(rows: list[str]) -> str:
    if not rows:
        return "<p>No legend entries.</p>"
    return '<div class="legend-grid">' + "".join(rows) + "</div>"

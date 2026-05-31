from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import node_signature, stable_value
from bijux_phylogenetics.biogeography.migration.migration_event_review import (
    GeographicMigrationEventReport,
)
from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.phylogeography.continuous_coordinates import (
    PhylogeographicCoordinateReport,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

from .contracts import (
    GeographicMapExclusionRow,
    GeographicMapLineRow,
    GeographicMapMarkerRow,
    GeographicMapSummary,
)

_EARTH_RADIUS_KM = 6371.0088


def validate_depth_filters(
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


def continuous_depths(
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


def branch_midpoint_depth(
    parent_depth: float | None,
    child_depth: float | None,
) -> float | None:
    if parent_depth is None or child_depth is None:
        return None
    return stable_value((parent_depth + child_depth) / 2.0)


def apply_depth_filter(
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


def active_line_counts(rows: list[GeographicMapLineRow]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        if not row.visible:
            continue
        counts[row.source_label] = counts.get(row.source_label, 0) + 1
        counts[row.target_label] = counts.get(row.target_label, 0) + 1
    return counts


def load_region_centroids(
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


def build_map_summary(
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


def great_circle_km(
    left_latitude: float,
    left_longitude: float,
    right_latitude: float,
    right_longitude: float,
) -> float:
    left_phi = math.radians(left_latitude)
    right_phi = math.radians(right_latitude)
    delta_phi = math.radians(right_latitude - left_latitude)
    delta_lambda = math.radians(right_longitude - left_longitude)
    distance_component = (
        math.sin(delta_phi / 2.0) ** 2
        + math.cos(left_phi) * math.cos(right_phi) * math.sin(delta_lambda / 2.0) ** 2
    )
    return stable_value(
        2.0 * _EARTH_RADIUS_KM * math.asin(math.sqrt(distance_component))
    )


def format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(stable_value(value))


def build_discrete_tip_markers(
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


def build_discrete_node_markers(
    report: GeographicStateModelReport,
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


def build_discrete_line_rows(
    report: GeographicMigrationEventReport,
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
                distance_km=great_circle_km(
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
        apply_depth_filter(
            rows,
            minimum_midpoint_depth=minimum_midpoint_depth,
            maximum_midpoint_depth=maximum_midpoint_depth,
        ),
        exclusions,
    )

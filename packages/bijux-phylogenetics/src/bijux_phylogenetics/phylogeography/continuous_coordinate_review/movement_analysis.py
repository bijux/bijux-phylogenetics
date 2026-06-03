from __future__ import annotations

import math

from bijux_phylogenetics.ancestral.continuous import ContinuousAncestralReport
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.contracts import (
    CoordinateEstimateRow,
    CoordinateMovementBranchRow,
    CoordinateMovementOutlierRow,
    CoordinateMovementSummary,
)
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.inputs import (
    PreparedCoordinateDataset,
)

_EARTH_RADIUS_KM = 6371.0088


def build_estimate_rows(
    latitude_report: ContinuousAncestralReport,
    longitude_report: ContinuousAncestralReport,
) -> list[CoordinateEstimateRow]:
    longitude_by_node = {
        estimate.node: estimate for estimate in longitude_report.estimates
    }
    root_node = latitude_report.estimates[0].node
    rows: list[CoordinateEstimateRow] = []
    for latitude_estimate in latitude_report.estimates:
        longitude_estimate = longitude_by_node[latitude_estimate.node]
        confidence = stable_float(
            min(latitude_estimate.confidence, longitude_estimate.confidence)
        )
        unstable = latitude_estimate.unstable or longitude_estimate.unstable
        rows.append(
            CoordinateEstimateRow(
                node=latitude_estimate.node,
                node_name=latitude_estimate.node_name,
                is_tip=latitude_estimate.is_tip,
                descendant_taxa=list(latitude_estimate.descendant_taxa),
                latitude=stable_float(latitude_estimate.estimate),
                longitude=stable_float(longitude_estimate.estimate),
                latitude_standard_error=stable_float(latitude_estimate.standard_error),
                longitude_standard_error=stable_float(
                    longitude_estimate.standard_error
                ),
                radial_standard_error_km=stable_float(
                    _radial_standard_error_km(
                        latitude=latitude_estimate.estimate,
                        latitude_standard_error=latitude_estimate.standard_error,
                        longitude_standard_error=longitude_estimate.standard_error,
                    )
                ),
                lower_95_latitude=stable_float(latitude_estimate.lower_95_interval),
                upper_95_latitude=stable_float(latitude_estimate.upper_95_interval),
                lower_95_longitude=stable_float(longitude_estimate.lower_95_interval),
                upper_95_longitude=stable_float(longitude_estimate.upper_95_interval),
                confidence=confidence,
                unstable=unstable,
                is_root=latitude_estimate.node == root_node,
            )
        )
    return rows


def build_branch_rows(
    latitude_report: ContinuousAncestralReport,
    longitude_report: ContinuousAncestralReport,
) -> list[CoordinateMovementBranchRow]:
    from bijux_phylogenetics.io.newick import loads_newick

    tree = loads_newick(latitude_report.analysis_tree_newick)
    estimates_by_node = {
        row.node: row for row in build_estimate_rows(latitude_report, longitude_report)
    }
    branch_rows: list[CoordinateMovementBranchRow] = []

    def visit(node) -> None:
        parent_estimate = estimates_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimates_by_node[_node_signature(child)]
            great_circle_km = stable_float(
                _haversine_km(
                    parent_estimate.latitude,
                    parent_estimate.longitude,
                    child_estimate.latitude,
                    child_estimate.longitude,
                )
            )
            branch_length = (
                None if child.branch_length is None else float(child.branch_length)
            )
            branch_rate_km_per_unit = (
                stable_float(great_circle_km / branch_length)
                if branch_length is not None and branch_length > 0.0
                else None
            )
            impossible_jump = (
                branch_length is not None
                and branch_length <= 0.0
                and great_circle_km > 1.0
            )
            branch_rows.append(
                CoordinateMovementBranchRow(
                    branch_id=child_estimate.node,
                    parent_node=parent_estimate.node,
                    child_node=child_estimate.node,
                    child_descendant_taxa=list(child_estimate.descendant_taxa),
                    branch_length=branch_length,
                    parent_latitude=parent_estimate.latitude,
                    parent_longitude=parent_estimate.longitude,
                    child_latitude=child_estimate.latitude,
                    child_longitude=child_estimate.longitude,
                    great_circle_km=great_circle_km,
                    branch_rate_km_per_unit=branch_rate_km_per_unit,
                    support=stable_float(
                        min(parent_estimate.confidence, child_estimate.confidence)
                    ),
                    impossible_jump=impossible_jump,
                    outlier_jump=False,
                    flag_codes=["zero_length_jump"] if impossible_jump else [],
                )
            )
            visit(child)

    visit(tree.root)
    return branch_rows


def review_branch_outliers(
    branch_rows: list[CoordinateMovementBranchRow],
) -> tuple[list[CoordinateMovementOutlierRow], list[CoordinateMovementBranchRow]]:
    distances = [row.great_circle_km for row in branch_rows]
    rates = [
        row.branch_rate_km_per_unit
        for row in branch_rows
        if row.branch_rate_km_per_unit is not None
    ]
    median_distance = _median(distances)
    distance_threshold = _robust_outlier_threshold(distances)
    median_rate = _median(rates) if rates else None
    rate_threshold = _robust_outlier_threshold(rates) if rates else None
    outlier_rows: list[CoordinateMovementOutlierRow] = []
    reviewed_branch_rows: list[CoordinateMovementBranchRow] = []
    for row in branch_rows:
        flag_codes = list(row.flag_codes)
        distance_outlier = row.great_circle_km > distance_threshold
        rate_outlier = (
            row.branch_rate_km_per_unit is not None
            and rate_threshold is not None
            and row.branch_rate_km_per_unit > rate_threshold
        )
        if distance_outlier:
            flag_codes.append("distance_outlier")
        if rate_outlier:
            flag_codes.append("rate_outlier")
        reviewed_row = CoordinateMovementBranchRow(
            branch_id=row.branch_id,
            parent_node=row.parent_node,
            child_node=row.child_node,
            child_descendant_taxa=row.child_descendant_taxa,
            branch_length=row.branch_length,
            parent_latitude=row.parent_latitude,
            parent_longitude=row.parent_longitude,
            child_latitude=row.child_latitude,
            child_longitude=row.child_longitude,
            great_circle_km=row.great_circle_km,
            branch_rate_km_per_unit=row.branch_rate_km_per_unit,
            support=row.support,
            impossible_jump=row.impossible_jump,
            outlier_jump=distance_outlier or rate_outlier,
            flag_codes=flag_codes,
        )
        reviewed_branch_rows.append(reviewed_row)
        if reviewed_row.impossible_jump or reviewed_row.outlier_jump:
            outlier_rows.append(
                CoordinateMovementOutlierRow(
                    branch_id=reviewed_row.branch_id,
                    parent_node=reviewed_row.parent_node,
                    child_node=reviewed_row.child_node,
                    child_descendant_taxa=reviewed_row.child_descendant_taxa,
                    great_circle_km=reviewed_row.great_circle_km,
                    branch_rate_km_per_unit=reviewed_row.branch_rate_km_per_unit,
                    median_distance_km=stable_float(median_distance),
                    distance_threshold_km=stable_float(distance_threshold),
                    median_rate_km_per_unit=(
                        None if median_rate is None else stable_float(median_rate)
                    ),
                    rate_threshold_km=(
                        None if rate_threshold is None else stable_float(rate_threshold)
                    ),
                    impossible_jump=reviewed_row.impossible_jump,
                    outlier_jump=reviewed_row.outlier_jump,
                    flag_codes=reviewed_row.flag_codes,
                )
            )
    return outlier_rows, reviewed_branch_rows


def build_summary(
    *,
    report: ContinuousAncestralReport,
    prepared: PreparedCoordinateDataset,
    latitude_column: str,
    longitude_column: str,
    estimate_rows: list[CoordinateEstimateRow],
    branch_rows: list[CoordinateMovementBranchRow],
    outlier_rows: list[CoordinateMovementOutlierRow],
    warnings: list[str],
) -> CoordinateMovementSummary:
    root_estimate = next(row for row in estimate_rows if row.is_root)
    return CoordinateMovementSummary(
        taxon_column=prepared.taxon_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        model=report.model,
        alpha=stable_float(report.alpha),
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(prepared.exclusion_rows),
        internal_node_count=sum(not row.is_tip for row in estimate_rows),
        weak_support_node_count=sum(
            (not row.is_tip) and row.confidence < 0.75 for row in estimate_rows
        ),
        outlier_jump_count=sum(row.outlier_jump for row in outlier_rows),
        impossible_jump_count=sum(row.impossible_jump for row in outlier_rows),
        flagged_branch_count=len(outlier_rows),
        maximum_jump_km=stable_float(
            max((row.great_circle_km for row in branch_rows), default=0.0)
        ),
        root_latitude=root_estimate.latitude,
        root_longitude=root_estimate.longitude,
        root_radial_standard_error_km=root_estimate.radial_standard_error_km,
        warning_count=len(warnings),
    )


def stable_float(value: float) -> float:
    normalized = round(float(value), 8)
    return 0.0 if normalized == -0.0 else normalized


def _haversine_km(
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
    return 2.0 * _EARTH_RADIUS_KM * math.asin(min(1.0, math.sqrt(a)))


def _radial_standard_error_km(
    *,
    latitude: float,
    latitude_standard_error: float,
    longitude_standard_error: float,
) -> float:
    lat_component = 111.32 * latitude_standard_error
    lon_component = 111.32 * math.cos(math.radians(latitude)) * longitude_standard_error
    return math.sqrt((lat_component**2) + (lon_component**2))


def _median(values: list[float | None]) -> float:
    numeric = sorted(float(value) for value in values if value is not None)
    if not numeric:
        return 0.0
    middle = len(numeric) // 2
    if len(numeric) % 2 == 1:
        return numeric[middle]
    return (numeric[middle - 1] + numeric[middle]) / 2.0


def _robust_outlier_threshold(values: list[float | None]) -> float:
    numeric = [float(value) for value in values if value is not None]
    if not numeric:
        return math.inf
    median_value = _median(numeric)
    deviations = [abs(value - median_value) for value in numeric]
    mad = _median(deviations)
    scale = max(mad * 1.4826, 1.0)
    return median_value + (3.0 * scale)


def _node_signature(node) -> str:
    if node.is_leaf():
        return node.name
    return "|".join(sorted(leaf.name for leaf in node.iter_leaves()))

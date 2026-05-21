from __future__ import annotations

from html import escape
import math
from pathlib import Path
import tempfile

from bijux_phylogenetics.ancestral.continuous import (
    ContinuousAncestralReport,
    reconstruct_continuous_ancestral_states,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.pruning import prune_tree_to_requested_taxa
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.contracts import (
    CoordinateEstimateRow,
    CoordinateMovementBranchRow,
    CoordinateMovementOutlierRow,
    CoordinateMovementSummary,
    CoordinateMovementVisualization,
    PhylogeographicCoordinateReport,
)
from bijux_phylogenetics.phylogeography.continuous_coordinate_review.inputs import (
    PreparedCoordinateDataset as _PreparedCoordinateDataset,
    prepare_coordinate_dataset as _prepare_coordinate_dataset,
    write_filtered_coordinate_table as _write_filtered_coordinate_table,
)
from bijux_phylogenetics.runtime.errors import AncestralReconstructionError

_EARTH_RADIUS_KM = 6371.0088


def summarize_continuous_phylogeography(
    tree_path: Path,
    table_path: Path,
    *,
    latitude_column: str,
    longitude_column: str,
    taxon_column: str | None = None,
    model: str = "brownian",
    alpha: float = 1.0,
) -> PhylogeographicCoordinateReport:
    """Reconstruct ancestral coordinates and review branchwise movement."""
    prepared = _prepare_coordinate_dataset(
        tree_path,
        table_path,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        taxon_column=taxon_column,
    )
    if len(prepared.included_taxa) < 2:
        raise AncestralReconstructionError(
            "continuous phylogeography requires at least two taxa with usable coordinates"
        )

    with tempfile.TemporaryDirectory(prefix="bijux-phylogeography-") as temp_dir:
        temp_root = Path(temp_dir)
        pruned_tree, _ = prune_tree_to_requested_taxa(tree_path, prepared.included_taxa)
        pruned_tree_path = temp_root / "analysis-tree.nwk"
        pruned_tree_path.write_text(f"{dumps_newick(pruned_tree)}\n", encoding="utf-8")
        filtered_table_path = temp_root / "coordinates.tsv"
        _write_filtered_coordinate_table(filtered_table_path, prepared)
        latitude_report = reconstruct_continuous_ancestral_states(
            pruned_tree_path,
            filtered_table_path,
            trait=latitude_column,
            taxon_column=prepared.taxon_column,
            model=model,
            alpha=alpha,
        )
        longitude_report = reconstruct_continuous_ancestral_states(
            pruned_tree_path,
            filtered_table_path,
            trait=longitude_column,
            taxon_column=prepared.taxon_column,
            model=model,
            alpha=alpha,
        )

    warnings = list(
        dict.fromkeys([*latitude_report.warnings, *longitude_report.warnings])
    )
    estimate_rows = _build_estimate_rows(latitude_report, longitude_report)
    branch_rows = _build_branch_rows(latitude_report, longitude_report)
    outlier_rows, reviewed_branch_rows = _review_branch_outliers(branch_rows)
    summary = _build_summary(
        report=latitude_report,
        prepared=prepared,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        estimate_rows=estimate_rows,
        branch_rows=reviewed_branch_rows,
        outlier_rows=outlier_rows,
        warnings=warnings,
    )
    return PhylogeographicCoordinateReport(
        tree_path=tree_path,
        table_path=table_path,
        taxon_column=prepared.taxon_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
        model=model,
        alpha=alpha,
        summary=summary,
        estimate_rows=estimate_rows,
        branch_rows=reviewed_branch_rows,
        outlier_rows=outlier_rows,
        exclusion_rows=prepared.exclusion_rows,
        warnings=warnings,
    )


def render_coordinate_movement_visualization(
    report: PhylogeographicCoordinateReport,
    *,
    out_path: Path,
    width: int = 960,
    height: int = 640,
) -> CoordinateMovementVisualization:
    """Render a coordinate-space movement visualization as SVG or HTML."""
    output_format = out_path.suffix.lower().lstrip(".")
    if output_format not in {"svg", "html"}:
        raise ValueError(
            "coordinate movement visualization output must end in .svg or .html"
        )
    svg_markup = _build_coordinate_svg(report, width=width, height=height)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if output_format == "svg":
        out_path.write_text(svg_markup, encoding="utf-8")
    else:
        out_path.write_text(
            _build_coordinate_html(report, svg_markup), encoding="utf-8"
        )
    return CoordinateMovementVisualization(
        output_path=out_path,
        format=output_format,
        width=width,
        height=height,
        highlighted_branch_count=sum(
            row.outlier_jump or row.impossible_jump for row in report.branch_rows
        ),
    )


def write_coordinate_movement_summary_table(
    path: Path,
    report: PhylogeographicCoordinateReport,
) -> Path:
    """Write one overall continuous phylogeography summary ledger."""
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
    """Write one coordinate estimate ledger for tips and internal nodes."""
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
    """Write one branchwise movement ledger."""
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
    """Write one flagged movement ledger."""
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
    """Write one excluded coordinate-row ledger."""
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


def _build_estimate_rows(
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
        confidence = _stable_float(
            min(latitude_estimate.confidence, longitude_estimate.confidence)
        )
        unstable = latitude_estimate.unstable or longitude_estimate.unstable
        rows.append(
            CoordinateEstimateRow(
                node=latitude_estimate.node,
                node_name=latitude_estimate.node_name,
                is_tip=latitude_estimate.is_tip,
                descendant_taxa=list(latitude_estimate.descendant_taxa),
                latitude=_stable_float(latitude_estimate.estimate),
                longitude=_stable_float(longitude_estimate.estimate),
                latitude_standard_error=_stable_float(latitude_estimate.standard_error),
                longitude_standard_error=_stable_float(
                    longitude_estimate.standard_error
                ),
                radial_standard_error_km=_stable_float(
                    _radial_standard_error_km(
                        latitude=latitude_estimate.estimate,
                        latitude_standard_error=latitude_estimate.standard_error,
                        longitude_standard_error=longitude_estimate.standard_error,
                    )
                ),
                lower_95_latitude=_stable_float(latitude_estimate.lower_95_interval),
                upper_95_latitude=_stable_float(latitude_estimate.upper_95_interval),
                lower_95_longitude=_stable_float(longitude_estimate.lower_95_interval),
                upper_95_longitude=_stable_float(longitude_estimate.upper_95_interval),
                confidence=confidence,
                unstable=unstable,
                is_root=latitude_estimate.node == root_node,
            )
        )
    return rows


def _build_branch_rows(
    latitude_report: ContinuousAncestralReport,
    longitude_report: ContinuousAncestralReport,
) -> list[CoordinateMovementBranchRow]:
    from bijux_phylogenetics.io.newick import loads_newick

    tree = loads_newick(latitude_report.analysis_tree_newick)
    estimates_by_node = {
        row.node: row for row in _build_estimate_rows(latitude_report, longitude_report)
    }
    branch_rows: list[CoordinateMovementBranchRow] = []

    def visit(node) -> None:
        parent_estimate = estimates_by_node[_node_signature(node)]
        for child in node.children:
            child_estimate = estimates_by_node[_node_signature(child)]
            great_circle_km = _stable_float(
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
                _stable_float(great_circle_km / branch_length)
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
                    support=_stable_float(
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


def _review_branch_outliers(
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
                    median_distance_km=_stable_float(median_distance),
                    distance_threshold_km=_stable_float(distance_threshold),
                    median_rate_km_per_unit=(
                        None if median_rate is None else _stable_float(median_rate)
                    ),
                    rate_threshold_km=(
                        None
                        if rate_threshold is None
                        else _stable_float(rate_threshold)
                    ),
                    impossible_jump=reviewed_row.impossible_jump,
                    outlier_jump=reviewed_row.outlier_jump,
                    flag_codes=reviewed_row.flag_codes,
                )
            )
    return outlier_rows, reviewed_branch_rows


def _build_summary(
    *,
    report: ContinuousAncestralReport,
    prepared: _PreparedCoordinateDataset,
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
        alpha=_stable_float(report.alpha),
        analyzed_taxon_count=report.taxon_count,
        excluded_taxon_count=len(prepared.exclusion_rows),
        internal_node_count=sum(not row.is_tip for row in estimate_rows),
        weak_support_node_count=sum(
            (not row.is_tip) and row.confidence < 0.75 for row in estimate_rows
        ),
        outlier_jump_count=sum(row.outlier_jump for row in outlier_rows),
        impossible_jump_count=sum(row.impossible_jump for row in outlier_rows),
        flagged_branch_count=len(outlier_rows),
        maximum_jump_km=_stable_float(
            max((row.great_circle_km for row in branch_rows), default=0.0)
        ),
        root_latitude=root_estimate.latitude,
        root_longitude=root_estimate.longitude,
        root_radial_standard_error_km=root_estimate.radial_standard_error_km,
        warning_count=len(warnings),
    )


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


def _build_coordinate_svg(
    report: PhylogeographicCoordinateReport,
    *,
    width: int,
    height: int,
) -> str:
    rows = report.estimate_rows
    latitudes = [row.latitude for row in rows]
    longitudes = [row.longitude for row in rows]
    min_latitude = min(latitudes)
    max_latitude = max(latitudes)
    min_longitude = min(longitudes)
    max_longitude = max(longitudes)
    padding = 48.0

    def project(longitude: float, latitude: float) -> tuple[float, float]:
        x_fraction = (
            0.5
            if max_longitude <= min_longitude
            else (longitude - min_longitude) / (max_longitude - min_longitude)
        )
        y_fraction = (
            0.5
            if max_latitude <= min_latitude
            else (latitude - min_latitude) / (max_latitude - min_latitude)
        )
        x = padding + x_fraction * (width - 2.0 * padding)
        y = padding + (1.0 - y_fraction) * (height - 2.0 * padding)
        return x, y

    estimate_by_node = {row.node: row for row in rows}
    branch_lines: list[str] = []
    for row in report.branch_rows:
        parent = estimate_by_node[row.parent_node]
        child = estimate_by_node[row.child_node]
        left_x, left_y = project(parent.longitude, parent.latitude)
        right_x, right_y = project(child.longitude, child.latitude)
        stroke = (
            "#b91c1c"
            if row.impossible_jump
            else "#c2410c"
            if row.outlier_jump
            else "#334155"
        )
        width_px = 3 if row.impossible_jump or row.outlier_jump else 1.75
        branch_lines.append(
            f'<line x1="{left_x:.1f}" y1="{left_y:.1f}" x2="{right_x:.1f}" y2="{right_y:.1f}" stroke="{stroke}" stroke-width="{width_px}" stroke-linecap="round" />'
        )

    point_marks: list[str] = []
    labels: list[str] = []
    for row in rows:
        x, y = project(row.longitude, row.latitude)
        if row.is_tip:
            fill = "#0f172a"
            radius = 4.5
        elif row.is_root:
            fill = "#a16207"
            radius = 5.0
        else:
            fill = "#0f766e"
            radius = 3.5
        point_marks.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="{radius:.1f}" fill="{fill}" opacity="0.95" />'
        )
        if row.node_name is not None:
            labels.append(
                f'<text x="{x + 6:.1f}" y="{y - 6:.1f}" font-size="12" fill="#0f172a">{escape(row.node_name)}</text>'
            )

    return "\n".join(
        [
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">',
            '<rect width="100%" height="100%" fill="#f8fafc" />',
            '<text x="24" y="30" font-size="20" fill="#0f172a">Coordinate-Space Movement Review</text>',
            f'<text x="24" y="{height - 18}" font-size="12" fill="#475569">Longitude increases left to right. Latitude increases bottom to top. This is a coordinate-space view, not a projected map.</text>',
            *branch_lines,
            *point_marks,
            *labels,
            "</svg>",
        ]
    )


def _build_coordinate_html(
    report: PhylogeographicCoordinateReport,
    svg_markup: str,
) -> str:
    summary = report.summary
    return "\n".join(
        [
            "<!doctype html>",
            '<html lang="en">',
            "<head>",
            '<meta charset="utf-8">',
            "<title>Continuous phylogeography review</title>",
            "<style>",
            "body { font-family: ui-sans-serif, system-ui, sans-serif; margin: 24px; color: #0f172a; background: #f8fafc; }",
            ".metrics { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; margin-bottom: 18px; }",
            ".card { background: white; border: 1px solid #cbd5e1; border-radius: 10px; padding: 12px; }",
            ".label { font-size: 12px; color: #475569; }",
            ".value { font-size: 18px; font-weight: 700; }",
            "svg { background: white; border: 1px solid #cbd5e1; border-radius: 10px; }",
            "</style>",
            "</head>",
            "<body>",
            "<h1>Continuous Phylogeography Review</h1>",
            '<div class="metrics">',
            f'<div class="card"><div class="label">Model</div><div class="value">{escape(summary.model)}</div></div>',
            f'<div class="card"><div class="label">Analyzed Taxa</div><div class="value">{summary.analyzed_taxon_count}</div></div>',
            f'<div class="card"><div class="label">Flagged Branches</div><div class="value">{summary.flagged_branch_count}</div></div>',
            f'<div class="card"><div class="label">Maximum Jump (km)</div><div class="value">{summary.maximum_jump_km}</div></div>',
            "</div>",
            svg_markup,
            "</body>",
            "</html>",
        ]
    )


def _format_optional_float(value: float | None) -> str:
    if value is None:
        return ""
    return str(_stable_float(value))


def _stable_float(value: float) -> float:
    normalized = round(float(value), 8)
    return 0.0 if normalized == -0.0 else normalized


def _node_signature(node) -> str:
    if node.is_leaf():
        return node.name
    return "|".join(sorted(leaf.name for leaf in node.iter_leaves()))

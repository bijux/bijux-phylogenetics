from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def _command_line_api() -> Any:
    import bijux_phylogenetics.command_line as command_line_api

    return command_line_api


def add_phylogeography_commands(subparsers: Any) -> None:
    phylogeography = subparsers.add_parser(
        get_command_spec("phylogeography").name,
        help=get_command_spec("phylogeography").summary,
    )
    phylogeography_subparsers = phylogeography.add_subparsers(
        dest="phylogeography_command",
        required=True,
    )
    phylogeography_coordinates = phylogeography_subparsers.add_parser(
        "coordinates",
        help="Reconstruct continuous geographic coordinates, review branch movement, and render coordinate-space movement output.",
    )
    phylogeography_coordinates.add_argument("tree", type=Path)
    phylogeography_coordinates.add_argument("table", type=Path)
    phylogeography_coordinates.add_argument("--latitude-column", required=True)
    phylogeography_coordinates.add_argument("--longitude-column", required=True)
    phylogeography_coordinates.add_argument("--taxon-column")
    phylogeography_coordinates.add_argument(
        "--model",
        choices=("brownian", "ou"),
        default="brownian",
    )
    phylogeography_coordinates.add_argument("--alpha", type=float, default=1.0)
    phylogeography_coordinates.add_argument("--summary-out", type=Path)
    phylogeography_coordinates.add_argument("--estimates-out", type=Path)
    phylogeography_coordinates.add_argument("--branches-out", type=Path)
    phylogeography_coordinates.add_argument("--outliers-out", type=Path)
    phylogeography_coordinates.add_argument("--exclusions-out", type=Path)
    phylogeography_coordinates.add_argument("--visualization-out", type=Path)
    phylogeography_coordinates.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogeography review as JSON.",
    )
    _add_manifest_argument(phylogeography_coordinates)

    phylogeography_coordinates_map = phylogeography_subparsers.add_parser(
        "coordinates-map",
        help="Render one HTML world map from continuous geographic coordinate reconstruction.",
    )
    phylogeography_coordinates_map.add_argument("tree", type=Path)
    phylogeography_coordinates_map.add_argument("table", type=Path)
    phylogeography_coordinates_map.add_argument("--latitude-column", required=True)
    phylogeography_coordinates_map.add_argument("--longitude-column", required=True)
    phylogeography_coordinates_map.add_argument("--taxon-column")
    phylogeography_coordinates_map.add_argument(
        "--model",
        choices=("brownian", "ou"),
        default="brownian",
    )
    phylogeography_coordinates_map.add_argument("--alpha", type=float, default=1.0)
    phylogeography_coordinates_map.add_argument(
        "--minimum-midpoint-depth",
        type=float,
    )
    phylogeography_coordinates_map.add_argument(
        "--maximum-midpoint-depth",
        type=float,
    )
    phylogeography_coordinates_map.add_argument("--summary-out", type=Path)
    phylogeography_coordinates_map.add_argument("--markers-out", type=Path)
    phylogeography_coordinates_map.add_argument("--lines-out", type=Path)
    phylogeography_coordinates_map.add_argument("--exclusions-out", type=Path)
    phylogeography_coordinates_map.add_argument("--html-out", type=Path)
    phylogeography_coordinates_map.add_argument(
        "--json",
        action="store_true",
        help="Emit the mapped phylogeography review as JSON.",
    )
    _add_manifest_argument(phylogeography_coordinates_map)

    phylogeography_regions_map = phylogeography_subparsers.add_parser(
        "regions-map",
        help="Render one HTML world map from discrete ancestral geographic region reconstruction.",
    )
    phylogeography_regions_map.add_argument("tree", type=Path)
    phylogeography_regions_map.add_argument("table", type=Path)
    phylogeography_regions_map.add_argument("--trait", required=True)
    phylogeography_regions_map.add_argument("--centroids", type=Path, required=True)
    phylogeography_regions_map.add_argument("--taxon-column")
    phylogeography_regions_map.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    phylogeography_regions_map.add_argument(
        "--region-column",
        default="region",
    )
    phylogeography_regions_map.add_argument(
        "--latitude-column",
        default="latitude",
    )
    phylogeography_regions_map.add_argument(
        "--longitude-column",
        default="longitude",
    )
    phylogeography_regions_map.add_argument(
        "--minimum-midpoint-depth",
        type=float,
    )
    phylogeography_regions_map.add_argument(
        "--maximum-midpoint-depth",
        type=float,
    )
    phylogeography_regions_map.add_argument("--summary-out", type=Path)
    phylogeography_regions_map.add_argument("--markers-out", type=Path)
    phylogeography_regions_map.add_argument("--lines-out", type=Path)
    phylogeography_regions_map.add_argument("--exclusions-out", type=Path)
    phylogeography_regions_map.add_argument("--html-out", type=Path)
    phylogeography_regions_map.add_argument(
        "--json",
        action="store_true",
        help="Emit the mapped regional reconstruction review as JSON.",
    )
    _add_manifest_argument(phylogeography_regions_map)


def run_phylogeography_command(args: Any) -> int:
    command_line_api = _command_line_api()
    if args.phylogeography_command == "coordinates":
        report = command_line_api.summarize_continuous_phylogeography(
            args.tree,
            args.table,
            latitude_column=args.latitude_column,
            longitude_column=args.longitude_column,
            taxon_column=args.taxon_column,
            model=args.model,
            alpha=args.alpha,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                command_line_api.write_coordinate_movement_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.estimates_out is not None:
            outputs.append(
                command_line_api.write_coordinate_estimate_table(
                    args.estimates_out,
                    report,
                )
            )
        if args.branches_out is not None:
            outputs.append(
                command_line_api.write_coordinate_movement_branch_table(
                    args.branches_out,
                    report,
                )
            )
        if args.outliers_out is not None:
            outputs.append(
                command_line_api.write_coordinate_movement_outlier_table(
                    args.outliers_out,
                    report,
                )
            )
        if args.exclusions_out is not None:
            outputs.append(
                command_line_api.write_coordinate_movement_exclusion_table(
                    args.exclusions_out,
                    report,
                )
            )
        if args.visualization_out is not None:
            outputs.append(
                command_line_api.render_coordinate_movement_visualization(
                    report,
                    out_path=args.visualization_out,
                ).output_path
            )
        outputs = _finalize_outputs(
            args,
            command="phylogeography",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="phylogeography",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model": report.model,
                    "analyzed_taxon_count": report.summary.analyzed_taxon_count,
                    "outlier_jump_count": report.summary.outlier_jump_count,
                    "impossible_jump_count": report.summary.impossible_jump_count,
                    "flagged_branch_count": report.summary.flagged_branch_count,
                    "maximum_jump_km": report.summary.maximum_jump_km,
                    "excluded_taxon_count": report.summary.excluded_taxon_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.phylogeography_command == "coordinates-map":
        report = command_line_api.summarize_continuous_phylogeography_map(
            args.tree,
            args.table,
            latitude_column=args.latitude_column,
            longitude_column=args.longitude_column,
            taxon_column=args.taxon_column,
            model=args.model,
            alpha=args.alpha,
            minimum_midpoint_depth=args.minimum_midpoint_depth,
            maximum_midpoint_depth=args.maximum_midpoint_depth,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                command_line_api.write_geographic_map_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.markers_out is not None:
            outputs.append(
                command_line_api.write_geographic_map_marker_table(
                    args.markers_out,
                    report,
                )
            )
        if args.lines_out is not None:
            outputs.append(
                command_line_api.write_geographic_map_line_table(
                    args.lines_out,
                    report,
                )
            )
        if args.exclusions_out is not None:
            outputs.append(
                command_line_api.write_geographic_map_exclusion_table(
                    args.exclusions_out,
                    report,
                )
            )
        if args.html_out is not None:
            outputs.append(
                command_line_api.render_geographic_map_html(
                    report,
                    out_path=args.html_out,
                ).output_path
            )
        outputs = _finalize_outputs(
            args,
            command="phylogeography",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="phylogeography",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "map_mode": report.summary.mode,
                    "model": report.summary.model,
                    "tip_marker_count": report.summary.tip_marker_count,
                    "internal_marker_count": (
                        report.summary.internal_marker_count
                        + report.summary.root_marker_count
                    ),
                    "line_count": report.summary.line_count,
                    "visible_line_count": report.summary.visible_line_count,
                    "time_filter_applied": report.summary.time_filter_applied,
                    "excluded_record_count": report.summary.excluded_record_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    report = command_line_api.summarize_discrete_region_map(
        args.tree,
        args.table,
        trait=args.trait,
        centroids_path=args.centroids,
        taxon_column=args.taxon_column,
        model=args.model,
        region_column=args.region_column,
        latitude_column=args.latitude_column,
        longitude_column=args.longitude_column,
        minimum_midpoint_depth=args.minimum_midpoint_depth,
        maximum_midpoint_depth=args.maximum_midpoint_depth,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            command_line_api.write_geographic_map_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.markers_out is not None:
        outputs.append(
            command_line_api.write_geographic_map_marker_table(
                args.markers_out,
                report,
            )
        )
    if args.lines_out is not None:
        outputs.append(
            command_line_api.write_geographic_map_line_table(
                args.lines_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            command_line_api.write_geographic_map_exclusion_table(
                args.exclusions_out,
                report,
            )
        )
    if args.html_out is not None:
        outputs.append(
            command_line_api.render_geographic_map_html(
                report,
                out_path=args.html_out,
            ).output_path
        )
    inputs = [args.tree, args.table, args.centroids]
    outputs = _finalize_outputs(
        args,
        command="phylogeography",
        inputs=inputs,
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="phylogeography",
            inputs=inputs,
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "map_mode": report.summary.mode,
                "model": report.summary.model,
                "tip_marker_count": report.summary.tip_marker_count,
                "internal_marker_count": (
                    report.summary.internal_marker_count
                    + report.summary.root_marker_count
                ),
                "line_count": report.summary.line_count,
                "visible_line_count": report.summary.visible_line_count,
                "time_filter_applied": report.summary.time_filter_applied,
                "excluded_record_count": report.summary.excluded_record_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

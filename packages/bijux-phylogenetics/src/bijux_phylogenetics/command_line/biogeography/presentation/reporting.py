from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_biogeography_report_command(biogeography_subparsers: Any) -> None:
    biogeography_report = biogeography_subparsers.add_parser(
        "report",
        help="Build a full biogeography report package with region counts, transition evidence, ancestral-region tree, and map output.",
    )
    biogeography_report.add_argument("tree", type=Path)
    biogeography_report.add_argument("table", type=Path)
    biogeography_report.add_argument("--trait", required=True)
    biogeography_report.add_argument("centroids", type=Path)
    biogeography_report.add_argument("--taxon-column")
    biogeography_report.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_report.add_argument(
        "--region-column",
        default="region",
        help="Region key column in the centroid table.",
    )
    biogeography_report.add_argument(
        "--latitude-column",
        default="latitude",
        help="Latitude column in the centroid table.",
    )
    biogeography_report.add_argument(
        "--longitude-column",
        default="longitude",
        help="Longitude column in the centroid table.",
    )
    biogeography_report.add_argument("--out-dir", required=True, type=Path)
    biogeography_report.add_argument(
        "--json", action="store_true", help="Emit the biogeography package as JSON."
    )
    _add_manifest_argument(biogeography_report)


def run_biogeography_report_command(args: Any) -> int | None:
    if args.biogeography_command != "report":
        return None

    cli_api = command_line_api()
    result = cli_api.build_biogeography_report_package(
        tree_path=args.tree,
        traits_path=args.table,
        centroids_path=args.centroids,
        trait=args.trait,
        out_dir=args.out_dir,
        taxon_column=args.taxon_column,
        model=args.model,
        region_column=args.region_column,
        latitude_column=args.latitude_column,
        longitude_column=args.longitude_column,
    )
    outputs = _finalize_outputs(
        args,
        command="biogeography",
        inputs=[args.tree, args.table, args.centroids],
        outputs=[
            result.report_path,
            result.tree_figure_path,
            result.map_path,
            result.legend_path,
            result.caption_path,
            result.summary_table_path,
            result.region_count_table_path,
            result.node_table_path,
            result.transition_matrix_path,
            result.event_table_path,
            result.map_marker_table_path,
            result.map_line_table_path,
            result.exclusion_table_path,
            result.manifest_path,
            result.reproducibility_manifest_path,
        ],
    )
    _print_result(
        build_command_result(
            command="biogeography",
            inputs=[args.tree, args.table, args.centroids],
            outputs=outputs,
            warnings=result.warnings,
            metrics={
                "report_kind": "biogeography-report-package",
                "model": result.state_report.model,
                "output_dir": str(result.output_dir),
                "artifact_count": 15,
                "observed_region_count": (
                    result.state_report.summary.observed_region_count
                ),
                "transition_rate_row_count": (
                    result.state_report.summary.transition_rate_row_count
                ),
                "event_count": result.event_report.summary.event_count,
                "visible_map_line_count": result.map_report.summary.visible_line_count,
                "publication_ready": result.audit.publication_ready,
                "legend_entry_count": result.audit.legend_entry_count,
                "caption_ready": result.audit.caption_ready,
                "rendered_internal_pie_count": (
                    result.audit.rendered_internal_pie_count
                ),
            },
            data=result,
        ),
        json_output=args.json,
    )
    return 0

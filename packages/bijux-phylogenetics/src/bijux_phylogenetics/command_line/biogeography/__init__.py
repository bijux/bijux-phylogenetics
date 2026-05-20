from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import command_line_api
from .state_models import (
    add_biogeography_state_model_commands,
    run_biogeography_state_model_command,
)


def add_biogeography_commands(subparsers: Any) -> None:
    biogeography = subparsers.add_parser(
        get_command_spec("biogeography").name,
        help=get_command_spec("biogeography").summary,
    )
    biogeography_subparsers = biogeography.add_subparsers(
        dest="biogeography_command",
        required=True,
    )
    add_biogeography_state_model_commands(biogeography_subparsers)

    biogeography_chronology = biogeography_subparsers.add_parser(
        "chronology",
        help="Place inferred geographic transitions into dated-tree age context with automatic age bins.",
    )
    biogeography_chronology.add_argument("tree", type=Path)
    biogeography_chronology.add_argument("table", type=Path)
    biogeography_chronology.add_argument("--trait", required=True)
    biogeography_chronology.add_argument("--taxon-column")
    biogeography_chronology.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_chronology.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_chronology.add_argument(
        "--time-bin-count",
        type=int,
        default=4,
        help="Number of equal-width age bins between the tips and the root age.",
    )
    biogeography_chronology.add_argument("--summary-out", type=Path)
    biogeography_chronology.add_argument("--nodes-out", type=Path)
    biogeography_chronology.add_argument("--events-out", type=Path)
    biogeography_chronology.add_argument("--bins-out", type=Path)
    biogeography_chronology.add_argument("--exclusions-out", type=Path)
    biogeography_chronology.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_chronology)

    biogeography_events = biogeography_subparsers.add_parser(
        "events",
        help="Extract inferred geographic movement events on one tree or across a retained tree set.",
    )
    biogeography_events.add_argument("tree", type=Path)
    biogeography_events.add_argument("table", type=Path)
    biogeography_events.add_argument("--trait", required=True)
    biogeography_events.add_argument("--taxon-column")
    biogeography_events.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_events.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_events.add_argument(
        "--tree-set",
        action="store_true",
        help="Interpret the input tree path as a posterior or bootstrap tree set.",
    )
    biogeography_events.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Fraction of leading trees to discard before tree-set event summarization.",
    )
    biogeography_events.add_argument("--summary-out", type=Path)
    biogeography_events.add_argument("--events-out", type=Path)
    biogeography_events.add_argument("--trees-out", type=Path)
    biogeography_events.add_argument("--event-summaries-out", type=Path)
    biogeography_events.add_argument("--exclusions-out", type=Path)
    biogeography_events.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_events)

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


def run_biogeography_command(args: Any) -> int:
    state_model_exit_code = run_biogeography_state_model_command(args)
    if state_model_exit_code is not None:
        return state_model_exit_code

    cli_api = command_line_api()
    if args.biogeography_command == "chronology":
        report = cli_api.summarize_biogeographic_transition_chronology(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            allowed_regions=cli_api._split_csv_values(args.allowed_regions)
            or None,
            time_bin_count=args.time_bin_count,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                cli_api.write_dated_biogeography_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.nodes_out is not None:
            outputs.append(
                cli_api.write_dated_biogeography_node_table(
                    args.nodes_out,
                    report,
                )
            )
        if args.events_out is not None:
            outputs.append(
                cli_api.write_dated_biogeography_event_table(
                    args.events_out,
                    report,
                )
            )
        if args.bins_out is not None:
            outputs.append(
                cli_api.write_dated_biogeography_time_bin_table(
                    args.bins_out,
                    report,
                )
            )
        if args.exclusions_out is not None:
            outputs.append(
                cli_api.write_dated_biogeography_exclusion_table(
                    args.exclusions_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="biogeography",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="biogeography",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "model": report.summary.model,
                    "tree_is_time_scaled": report.summary.tree_is_time_scaled,
                    "root_age": report.summary.root_age,
                    "event_count": report.summary.event_count,
                    "time_bin_count": report.summary.time_bin_count,
                    "high_uncertainty_bin_count": (
                        report.summary.high_uncertainty_bin_count
                    ),
                    "excluded_taxon_count": report.summary.excluded_taxon_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.biogeography_command == "events":
        outputs: list[Path | str] = []
        if args.tree_set:
            report = cli_api.summarize_geographic_migration_event_tree_set(
                args.tree,
                args.table,
                trait=args.trait,
                taxon_column=args.taxon_column,
                model=args.model,
                allowed_regions=cli_api._split_csv_values(
                    args.allowed_regions
                )
                or None,
                burnin_fraction=args.burnin_fraction,
            )
            if args.summary_out is not None:
                outputs.append(
                    cli_api.write_geographic_migration_tree_set_summary_table(
                        args.summary_out,
                        report,
                    )
                )
            if args.events_out is not None:
                outputs.append(
                    cli_api.write_geographic_migration_tree_set_event_table(
                        args.events_out,
                        report,
                    )
                )
            if args.trees_out is not None:
                outputs.append(
                    cli_api.write_geographic_migration_tree_set_tree_table(
                        args.trees_out,
                        report,
                    )
                )
            if args.event_summaries_out is not None:
                outputs.append(
                    cli_api.write_geographic_migration_tree_set_event_summary_table(
                        args.event_summaries_out,
                        report,
                    )
                )
            if args.exclusions_out is not None:
                outputs.append(
                    cli_api.write_geographic_migration_tree_set_exclusion_table(
                        args.exclusions_out,
                        report,
                    )
                )
            outputs = _finalize_outputs(
                args,
                command="biogeography",
                inputs=[args.tree, args.table],
                outputs=outputs,
            )
            _print_result(
                build_command_result(
                    command="biogeography",
                    inputs=[args.tree, args.table],
                    outputs=outputs,
                    warnings=report.warnings,
                    metrics={
                        "report_mode": "tree_set",
                        "model": report.model,
                        "kept_tree_count": report.summary.kept_tree_count,
                        "event_row_count": report.summary.event_row_count,
                        "event_summary_count": report.summary.event_summary_count,
                        "topology_sensitive_event_count": (
                            report.summary.topology_sensitive_event_count
                        ),
                        "excluded_taxon_count": report.summary.excluded_taxon_count,
                        "warning_count": report.summary.warning_count,
                    },
                    data=report,
                ),
                json_output=args.json,
            )
            return 0
        report = cli_api.summarize_geographic_migration_events(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            allowed_regions=cli_api._split_csv_values(args.allowed_regions)
            or None,
        )
        if args.summary_out is not None:
            outputs.append(
                cli_api.write_geographic_migration_event_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.events_out is not None:
            outputs.append(
                cli_api.write_geographic_migration_event_table(
                    args.events_out,
                    report,
                )
            )
        if args.exclusions_out is not None:
            outputs.append(
                cli_api.write_geographic_migration_exclusion_table(
                    args.exclusions_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="biogeography",
            inputs=[args.tree, args.table],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="biogeography",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "report_mode": "single_tree",
                    "model": report.model,
                    "event_count": report.summary.event_count,
                    "strongly_supported_event_count": (
                        report.summary.strongly_supported_event_count
                    ),
                    "mean_event_support": report.summary.mean_event_support,
                    "excluded_taxon_count": report.summary.excluded_taxon_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
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

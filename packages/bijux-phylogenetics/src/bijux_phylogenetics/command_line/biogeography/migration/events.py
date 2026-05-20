from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_biogeography_event_commands(biogeography_subparsers: Any) -> None:
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


def run_biogeography_event_command(args: Any) -> int | None:
    if args.biogeography_command != "events":
        return None

    cli_api = command_line_api()
    outputs: list[Path | str] = []
    if args.tree_set:
        report = cli_api.summarize_geographic_migration_event_tree_set(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
            model=args.model,
            allowed_regions=cli_api._split_csv_values(args.allowed_regions) or None,
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
        allowed_regions=cli_api._split_csv_values(args.allowed_regions) or None,
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

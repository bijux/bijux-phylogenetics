from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_biogeography_chronology_command(biogeography_subparsers: Any) -> None:
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


def run_biogeography_chronology_command(args: Any) -> int | None:
    if args.biogeography_command != "chronology":
        return None

    cli_api = command_line_api()
    report = cli_api.summarize_biogeographic_transition_chronology(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_regions=cli_api._split_csv_values(args.allowed_regions) or None,
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

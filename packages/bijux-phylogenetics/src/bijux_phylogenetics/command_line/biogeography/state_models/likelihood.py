from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_geographic_likelihood_command(biogeography_subparsers: Any) -> None:
    biogeography_model = biogeography_subparsers.add_parser(
        "model",
        help="Reconstruct ancestral geographic regions under an ER, SYM, or ARD transition model.",
    )
    biogeography_model.add_argument("tree", type=Path)
    biogeography_model.add_argument("table", type=Path)
    biogeography_model.add_argument("--trait", required=True)
    biogeography_model.add_argument("--taxon-column")
    biogeography_model.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_model.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_model.add_argument("--summary-out", type=Path)
    biogeography_model.add_argument("--nodes-out", type=Path)
    biogeography_model.add_argument("--rates-out", type=Path)
    biogeography_model.add_argument("--events-out", type=Path)
    biogeography_model.add_argument("--exclusions-out", type=Path)
    biogeography_model.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_model)


def run_geographic_likelihood_command(args: Any) -> int | None:
    if args.biogeography_command != "model":
        return None

    cli_api = command_line_api()
    report = cli_api.summarize_geographic_state_model(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_regions=cli_api._split_csv_values(args.allowed_regions) or None,
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            cli_api.write_geographic_state_summary_table(args.summary_out, report)
        )
    if args.nodes_out is not None:
        outputs.append(
            cli_api.write_geographic_region_probability_table(args.nodes_out, report)
        )
    if args.rates_out is not None:
        outputs.append(
            cli_api.write_geographic_transition_rate_table(args.rates_out, report)
        )
    if args.events_out is not None:
        outputs.append(
            cli_api.write_geographic_transition_event_table(args.events_out, report)
        )
    if args.exclusions_out is not None:
        outputs.append(
            cli_api.write_geographic_exclusion_table(args.exclusions_out, report)
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
                "model": report.model,
                "observed_region_count": report.summary.observed_region_count,
                "internal_node_count": report.summary.internal_node_count,
                "transition_rate_row_count": report.summary.transition_rate_row_count,
                "changed_branch_count": report.summary.changed_branch_count,
                "strongly_supported_transition_count": (
                    report.summary.strongly_supported_transition_count
                ),
                "excluded_taxon_count": report.summary.excluded_taxon_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result

from ..shared import command_line_api


def add_time_stratified_geography_command(biogeography_subparsers: Any) -> None:
    biogeography_time_stratified = biogeography_subparsers.add_parser(
        "time-stratified",
        help="Estimate interval-specific geographic transitions across explicit root-depth bins.",
    )
    biogeography_time_stratified.add_argument("tree", type=Path)
    biogeography_time_stratified.add_argument("table", type=Path)
    biogeography_time_stratified.add_argument("--trait", required=True)
    biogeography_time_stratified.add_argument("--taxon-column")
    biogeography_time_stratified.add_argument(
        "--model",
        choices=("er", "sym", "ard"),
        default="er",
    )
    biogeography_time_stratified.add_argument(
        "--allowed-regions",
        help="Comma-delimited explicit region vocabulary.",
    )
    biogeography_time_stratified.add_argument(
        "--time-bin",
        action="append",
        required=True,
        metavar="LABEL:START:END",
        help="Explicit root-depth interval in LABEL:START:END form. Repeat for multiple intervals.",
    )
    biogeography_time_stratified.add_argument("--summary-out", type=Path)
    biogeography_time_stratified.add_argument("--matrix-out", type=Path)
    biogeography_time_stratified.add_argument("--branches-out", type=Path)
    biogeography_time_stratified.add_argument("--exclusions-out", type=Path)
    biogeography_time_stratified.add_argument(
        "--json", action="store_true", help="Emit the biogeography review as JSON."
    )
    _add_manifest_argument(biogeography_time_stratified)


def run_time_stratified_geography_command(args: Any) -> int | None:
    if args.biogeography_command != "time-stratified":
        return None

    cli_api = command_line_api()
    report = cli_api.summarize_time_stratified_geographic_transitions(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_regions=cli_api._split_csv_values(args.allowed_regions) or None,
        time_bins=[
            cli_api._parse_time_bin_definition(raw_time_bin)
            for raw_time_bin in args.time_bin
        ],
    )
    outputs: list[Path | str] = []
    if args.summary_out is not None:
        outputs.append(
            cli_api.write_time_stratified_transition_summary_table(
                args.summary_out,
                report,
            )
        )
    if args.matrix_out is not None:
        outputs.append(
            cli_api.write_time_stratified_transition_matrix_table(
                args.matrix_out,
                report,
            )
        )
    if args.branches_out is not None:
        outputs.append(
            cli_api.write_time_stratified_branch_table(
                args.branches_out,
                report,
            )
        )
    if args.exclusions_out is not None:
        outputs.append(
            cli_api.write_time_stratified_exclusion_table(
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
                "model": report.model,
                "time_bin_count": report.summary.time_bin_count,
                "matrix_row_count": report.summary.matrix_row_count,
                "changed_branch_count": report.summary.changed_branch_count,
                "allocated_transition_weight_total": (
                    report.summary.allocated_transition_weight_total
                ),
                "excluded_taxon_count": report.summary.excluded_taxon_count,
                "warning_count": report.summary.warning_count,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

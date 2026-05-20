from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import write_supplementary_diversification_table
from bijux_phylogenetics.runtime.results import build_command_result


def add_diversification_supplementary_table_commands(
    report_subparsers: Any,
) -> None:
    report_supplementary_diversification_table = report_subparsers.add_parser(
        "supplementary-diversification-table",
        help="Write a supplementary diversification table with clade rates, model estimates, sampling correction, and warnings.",
    )
    report_supplementary_diversification_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_diversification_table.add_argument("--metadata", type=Path)
    report_supplementary_diversification_table.add_argument("--taxon-column")
    report_supplementary_diversification_table.add_argument("--sampling-column")
    report_supplementary_diversification_table.add_argument(
        "--clade-model",
        choices=("yule", "birth-death"),
        default="birth-death",
    )
    report_supplementary_diversification_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_diversification_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_diversification_table)


def run_diversification_supplementary_table_command(args: Any) -> int | None:
    if args.report_command != "supplementary-diversification-table":
        return None

    inputs = [args.tree]
    if args.metadata is not None:
        inputs.append(args.metadata)
    result = write_supplementary_diversification_table(
        args.out,
        tree_path=args.tree,
        metadata_path=args.metadata,
        taxon_column=args.taxon_column,
        sampling_column=args.sampling_column,
        clade_model=args.clade_model,
    )
    outputs = _finalize_outputs(
        args,
        command="report",
        inputs=inputs,
        outputs=[result.output_path],
    )
    if args.json:
        _print_result(
            build_command_result(
                command="report",
                inputs=inputs,
                outputs=outputs,
                warnings=[],
                metrics={
                    "row_count": result.row_count,
                    "better_model": result.better_model,
                    "clade_model": result.clade_model,
                    "high_clade_count": result.high_clade_count,
                    "low_clade_count": result.low_clade_count,
                    "warning_count": result.warning_count,
                    "sampling_metadata_complete": result.sampling_metadata_complete,
                },
                data=result,
            ),
            json_output=True,
        )
        return 0
    print(result.output_path)
    return 0

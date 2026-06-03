from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.reports import write_supplementary_comparative_model_table
from bijux_phylogenetics.runtime.results import build_command_result


def _parse_lambda_value(value: str | None) -> float | str:
    if value == "estimate" or value is None:
        return "estimate"
    return float(value)


def add_comparative_model_supplementary_table_commands(
    report_subparsers: Any,
) -> None:
    report_supplementary_comparative_model_table = report_subparsers.add_parser(
        "supplementary-comparative-model-table",
        help="Write a supplementary comparative-model table with coefficients, uncertainty, diagnostics, and exclusions.",
    )
    report_supplementary_comparative_model_table.add_argument(
        "--tree", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--traits", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--formula",
        dest="formulas",
        action="append",
        required=True,
        help="Add one comparative candidate formula. Repeat for each candidate model.",
    )
    report_supplementary_comparative_model_table.add_argument("--taxon-column")
    report_supplementary_comparative_model_table.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    report_supplementary_comparative_model_table.add_argument(
        "--out", required=True, type=Path
    )
    report_supplementary_comparative_model_table.add_argument(
        "--json", action="store_true", help="Emit the table write result as JSON."
    )
    _add_manifest_argument(report_supplementary_comparative_model_table)


def run_comparative_model_supplementary_table_command(args: Any) -> int | None:
    if args.report_command != "supplementary-comparative-model-table":
        return None

    result = write_supplementary_comparative_model_table(
        args.out,
        tree_path=args.tree,
        traits_path=args.traits,
        formulas=list(args.formulas),
        taxon_column=args.taxon_column,
        lambda_value=_parse_lambda_value(getattr(args, "lambda_value", None)),
    )
    inputs = [args.tree, args.traits]
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
                    "model_count": result.model_count,
                    "selected_formula": result.selected_formula,
                    "selected_criterion": result.selected_criterion,
                    "excluded_taxon_count": result.excluded_taxon_count,
                },
                data=result,
            ),
            json_output=True,
        )
        return 0
    print(result.output_path)
    return 0

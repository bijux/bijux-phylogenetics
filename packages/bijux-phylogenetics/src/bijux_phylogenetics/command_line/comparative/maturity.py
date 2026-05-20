from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.assessment import (
    assess_comparative_method_maturity,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_maturity_command(comparative_subparsers: Any) -> None:
    comparative_maturity = comparative_subparsers.add_parser(
        "maturity",
        help="Audit comparative residual diagnostics and sensitivity for one response trait workflow.",
    )
    comparative_maturity.add_argument("tree", type=Path)
    comparative_maturity.add_argument("table", type=Path)
    comparative_maturity.add_argument("--response")
    comparative_maturity.add_argument("--predictors", nargs="+")
    comparative_maturity.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_maturity.add_argument("--taxon-column")
    comparative_maturity.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_maturity.add_argument(
        "--json",
        action="store_true",
        help="Emit the maturity audit as JSON.",
    )
    _add_manifest_argument(comparative_maturity)


def run_comparative_maturity_command(args: Any) -> int | None:
    lambda_value: float | str
    if hasattr(args, "lambda_value"):
        if args.lambda_value == "estimate":
            lambda_value = "estimate"
        else:
            lambda_value = float(args.lambda_value)
    else:
        lambda_value = "estimate"

    if args.comparative_command != "maturity":
        return None

    report = assess_comparative_method_maturity(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "selected_model": report.selected_model,
                "residual_surface_count": len(report.residual_diagnostics),
                "influential_taxa": len(report.sensitivity.influential_taxa),
                "reference_validation_passed": report.reference_validation_passed,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

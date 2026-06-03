from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls import run_pgls_multiple_testing
from bijux_phylogenetics.runtime.results import build_command_result


def add_multiple_testing_pgls_commands(comparative_subparsers: Any) -> None:
    comparative_multiple_testing = comparative_subparsers.add_parser(
        "multiple-testing",
        help="Adjust PGLS coefficient p-values across many response traits.",
    )
    comparative_multiple_testing.add_argument("tree", type=Path)
    comparative_multiple_testing.add_argument("table", type=Path)
    comparative_multiple_testing.add_argument("--responses", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--predictors", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--taxon-column")
    comparative_multiple_testing.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multiple_testing.add_argument(
        "--json", action="store_true", help="Emit the correction report as JSON."
    )
    _add_manifest_argument(comparative_multiple_testing)


def run_multiple_testing_pgls_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "multiple-testing":
        return None

    report = run_pgls_multiple_testing(
        args.tree,
        args.table,
        responses=list(args.responses),
        predictors=list(args.predictors),
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
            metrics={
                "response_count": len(report.responses),
                "test_count": len(report.rows),
                "family_size": report.family_size,
                "raw_significant_count": report.raw_significant_count,
                "significant_count": sum(1 for row in report.rows if row.significant),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

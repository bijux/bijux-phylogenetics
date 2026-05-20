from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.pgls import run_pgls_multiple_testing
from bijux_phylogenetics.runtime.results import build_command_result

from .brownian import (
    add_brownian_pgls_commands,
    run_brownian_pgls_command,
)
from .covariance_audit import (
    add_covariance_audit_pgls_commands,
    run_covariance_audit_pgls_command,
)
from .estimated_lambda import (
    add_estimated_lambda_pgls_commands,
    run_estimated_lambda_pgls_command,
)
from .ou import add_ou_pgls_commands, run_ou_pgls_command


def add_comparative_pgls_commands(comparative_subparsers: Any) -> None:
    add_covariance_audit_pgls_commands(comparative_subparsers)
    add_estimated_lambda_pgls_commands(comparative_subparsers)
    add_brownian_pgls_commands(comparative_subparsers)
    add_ou_pgls_commands(comparative_subparsers)

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


def run_comparative_pgls_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    lambda_value: float | str
    if hasattr(args, "lambda_value"):
        if args.lambda_value == "estimate":
            lambda_value = "estimate"
        else:
            lambda_value = float(args.lambda_value)
    else:
        lambda_value = "estimate"
    covariance_audit_result = run_covariance_audit_pgls_command(
        args,
        lambda_value=lambda_value,
    )
    if covariance_audit_result is not None:
        return covariance_audit_result

    estimated_lambda_result = run_estimated_lambda_pgls_command(
        args,
        lambda_value=lambda_value,
    )
    if estimated_lambda_result is not None:
        return estimated_lambda_result

    brownian_result = run_brownian_pgls_command(args)
    if brownian_result is not None:
        return brownian_result

    ou_result = run_ou_pgls_command(args)
    if ou_result is not None:
        return ou_result
    if args.comparative_command == "multiple-testing":
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
    return None

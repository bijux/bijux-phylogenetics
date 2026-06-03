from __future__ import annotations

from typing import Any

from .comparisons import (
    add_regression_comparison_commands,
    run_regression_comparison_command,
)
from .logistic import (
    add_comparative_logistic_command,
    run_comparative_logistic_command,
)
from .model_selection import (
    add_regression_model_selection_command,
    run_regression_model_selection_command,
)
from .multivariate import (
    add_multivariate_regression_command,
    run_multivariate_regression_command,
)
from .pgls import add_comparative_pgls_commands, run_comparative_pgls_command
from .reporting import (
    add_regression_reporting_commands,
    run_regression_reporting_command,
)


def add_comparative_modeling_commands(comparative_subparsers: Any) -> None:
    add_regression_model_selection_command(comparative_subparsers)
    add_multivariate_regression_command(comparative_subparsers)
    add_regression_reporting_commands(comparative_subparsers)
    add_regression_comparison_commands(comparative_subparsers)


def run_comparative_modeling_command(
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

    model_selection_exit_code = run_regression_model_selection_command(
        args,
        lambda_value=lambda_value,
    )
    if model_selection_exit_code is not None:
        return model_selection_exit_code

    multivariate_exit_code = run_multivariate_regression_command(
        args,
        lambda_value=lambda_value,
    )
    if multivariate_exit_code is not None:
        return multivariate_exit_code

    reporting_exit_code = run_regression_reporting_command(
        args,
        lambda_value=lambda_value,
    )
    if reporting_exit_code is not None:
        return reporting_exit_code

    return run_regression_comparison_command(
        args,
        lambda_value=lambda_value,
    )


__all__ = [
    "add_comparative_logistic_command",
    "add_comparative_modeling_commands",
    "add_comparative_pgls_commands",
    "run_comparative_logistic_command",
    "run_comparative_modeling_command",
    "run_comparative_pgls_command",
]

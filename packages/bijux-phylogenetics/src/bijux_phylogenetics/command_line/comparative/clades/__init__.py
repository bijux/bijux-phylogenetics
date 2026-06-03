from __future__ import annotations

from typing import Any

from .posterior_trees import (
    add_posterior_tree_regression_commands,
    run_posterior_tree_regression_command,
)
from .residuals import add_clade_residual_commands, run_clade_residual_command
from .stability import add_clade_stability_commands, run_clade_stability_command


def add_comparative_clade_commands(comparative_subparsers: Any) -> None:
    add_clade_residual_commands(comparative_subparsers)
    add_clade_stability_commands(comparative_subparsers)
    add_posterior_tree_regression_commands(comparative_subparsers)


def run_comparative_clade_command(
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

    residual_exit_code = run_clade_residual_command(args, lambda_value=lambda_value)
    if residual_exit_code is not None:
        return residual_exit_code

    stability_exit_code = run_clade_stability_command(args, lambda_value=lambda_value)
    if stability_exit_code is not None:
        return stability_exit_code

    return run_posterior_tree_regression_command(args, lambda_value=lambda_value)


__all__ = [
    "add_comparative_clade_commands",
    "run_comparative_clade_command",
]

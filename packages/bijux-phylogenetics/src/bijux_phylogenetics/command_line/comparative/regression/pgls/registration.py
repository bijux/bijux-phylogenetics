from __future__ import annotations

from typing import Any

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
from .multiple_testing import (
    add_multiple_testing_pgls_commands,
    run_multiple_testing_pgls_command,
)
from .ou import add_ou_pgls_commands, run_ou_pgls_command


def add_comparative_pgls_commands(comparative_subparsers: Any) -> None:
    add_covariance_audit_pgls_commands(comparative_subparsers)
    add_estimated_lambda_pgls_commands(comparative_subparsers)
    add_brownian_pgls_commands(comparative_subparsers)
    add_ou_pgls_commands(comparative_subparsers)
    add_multiple_testing_pgls_commands(comparative_subparsers)


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
    return run_multiple_testing_pgls_command(args, lambda_value=lambda_value)

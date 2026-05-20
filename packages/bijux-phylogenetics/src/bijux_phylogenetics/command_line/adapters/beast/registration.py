from __future__ import annotations

from typing import Any

from .calibration_review import (
    add_beast_calibration_review_commands,
    run_beast_calibration_review_command,
)
from .diagnostics import (
    add_beast_diagnostic_commands,
    run_beast_diagnostic_command,
)
from .execution import (
    add_beast_execution_commands,
    run_beast_execution_command,
)
from .posterior_review import (
    add_beast_posterior_review_commands,
    run_beast_posterior_review_command,
)
from .posterior_trees import (
    add_beast_posterior_tree_commands,
    run_beast_posterior_tree_command,
)


def add_beast_adapter_commands(adapter_subparsers: Any) -> None:
    add_beast_execution_commands(adapter_subparsers)
    add_beast_calibration_review_commands(adapter_subparsers)
    add_beast_diagnostic_commands(adapter_subparsers)
    add_beast_posterior_tree_commands(adapter_subparsers)
    add_beast_posterior_review_commands(adapter_subparsers)


def run_beast_adapter_command(args: Any) -> int | None:
    if not str(args.adapter_command).startswith("beast-"):
        return None

    execution_result = run_beast_execution_command(args)
    if execution_result is not None:
        return execution_result

    calibration_review_result = run_beast_calibration_review_command(args)
    if calibration_review_result is not None:
        return calibration_review_result

    diagnostic_result = run_beast_diagnostic_command(args)
    if diagnostic_result is not None:
        return diagnostic_result

    posterior_tree_result = run_beast_posterior_tree_command(args)
    if posterior_tree_result is not None:
        return posterior_tree_result

    posterior_review_result = run_beast_posterior_review_command(args)
    if posterior_review_result is not None:
        return posterior_review_result

    return None

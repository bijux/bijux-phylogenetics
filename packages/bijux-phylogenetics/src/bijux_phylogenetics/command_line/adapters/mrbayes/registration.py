from __future__ import annotations

from typing import Any

from .diagnostics import (
    add_mrbayes_diagnostic_commands,
    run_mrbayes_diagnostic_command,
)
from .execution import (
    add_mrbayes_execution_commands,
    run_mrbayes_execution_command,
)
from .parameter_review import (
    add_mrbayes_parameter_review_commands,
    run_mrbayes_parameter_review_command,
)
from .posterior_report import (
    add_mrbayes_posterior_report_commands,
    run_mrbayes_posterior_report_command,
)
from .posterior_trees import (
    add_mrbayes_posterior_tree_commands,
    run_mrbayes_posterior_tree_command,
)


def add_mrbayes_adapter_commands(adapter_subparsers: Any) -> None:
    add_mrbayes_execution_commands(adapter_subparsers)
    add_mrbayes_posterior_tree_commands(adapter_subparsers)
    add_mrbayes_diagnostic_commands(adapter_subparsers)
    add_mrbayes_parameter_review_commands(adapter_subparsers)
    add_mrbayes_posterior_report_commands(adapter_subparsers)


def run_mrbayes_adapter_command(args: Any) -> int | None:
    if not str(args.adapter_command).startswith("mrbayes-"):
        return None

    execution_result = run_mrbayes_execution_command(args)
    if execution_result is not None:
        return execution_result

    posterior_tree_result = run_mrbayes_posterior_tree_command(args)
    if posterior_tree_result is not None:
        return posterior_tree_result

    diagnostic_result = run_mrbayes_diagnostic_command(args)
    if diagnostic_result is not None:
        return diagnostic_result

    parameter_review_result = run_mrbayes_parameter_review_command(args)
    if parameter_review_result is not None:
        return parameter_review_result

    return run_mrbayes_posterior_report_command(args)

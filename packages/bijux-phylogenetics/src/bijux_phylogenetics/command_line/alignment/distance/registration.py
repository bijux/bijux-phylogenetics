from __future__ import annotations

from typing import Any

from .diagnostics import (
    add_distance_diagnostic_commands,
    run_distance_diagnostic_command,
)
from .matrix import add_distance_matrix_command, run_distance_matrix_command
from .reporting import (
    add_distance_reporting_commands,
    run_distance_reporting_command,
)
from .sensitivity import (
    add_distance_sensitivity_commands,
    run_distance_sensitivity_command,
)
from .support import add_distance_support_commands, run_distance_support_command
from .trees import add_distance_tree_commands, run_distance_tree_command


def add_alignment_distance_commands(alignment_subparsers: Any) -> None:
    add_distance_matrix_command(alignment_subparsers)
    add_distance_diagnostic_commands(alignment_subparsers)
    add_distance_tree_commands(alignment_subparsers)
    add_distance_support_commands(alignment_subparsers)
    add_distance_sensitivity_commands(alignment_subparsers)
    add_distance_reporting_commands(alignment_subparsers)


def run_alignment_distance_command(args: Any) -> int | None:
    matrix_result = run_distance_matrix_command(args)
    if matrix_result is not None:
        return matrix_result

    diagnostic_result = run_distance_diagnostic_command(args)
    if diagnostic_result is not None:
        return diagnostic_result

    tree_result = run_distance_tree_command(args)
    if tree_result is not None:
        return tree_result

    support_result = run_distance_support_command(args)
    if support_result is not None:
        return support_result

    sensitivity_result = run_distance_sensitivity_command(args)
    if sensitivity_result is not None:
        return sensitivity_result

    return run_distance_reporting_command(args)

from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec

from .comparison import (
    add_tree_set_comparison_commands,
    run_tree_set_comparison_command,
)
from .presentation import (
    add_tree_set_presentation_commands,
    run_tree_set_presentation_command,
)
from .structure import (
    add_tree_set_structure_commands,
    run_tree_set_structure_command,
)
from .summary import (
    add_tree_set_summary_commands,
    run_tree_set_summary_command,
)
from .uncertainty import (
    add_tree_set_uncertainty_commands,
    run_tree_set_uncertainty_command,
)


def add_tree_set_commands(subparsers: Any) -> None:
    tree_set = subparsers.add_parser(
        get_command_spec("tree-set").name,
        help=get_command_spec("tree-set").summary,
    )
    tree_set_subparsers = tree_set.add_subparsers(
        dest="tree_set_command", required=True
    )
    add_tree_set_summary_commands(tree_set_subparsers)
    add_tree_set_structure_commands(tree_set_subparsers)
    add_tree_set_uncertainty_commands(tree_set_subparsers)
    add_tree_set_comparison_commands(tree_set_subparsers)
    add_tree_set_presentation_commands(tree_set_subparsers)


def run_tree_set_command(args: Any) -> int:
    summary_result = run_tree_set_summary_command(args)
    if summary_result is not None:
        return summary_result

    structure_result = run_tree_set_structure_command(args)
    if structure_result is not None:
        return structure_result

    uncertainty_result = run_tree_set_uncertainty_command(args)
    if uncertainty_result is not None:
        return uncertainty_result

    comparison_result = run_tree_set_comparison_command(args)
    if comparison_result is not None:
        return comparison_result

    presentation_result = run_tree_set_presentation_command(args)
    if presentation_result is not None:
        return presentation_result

    raise ValueError(f"unsupported tree-set command: {args.tree_set_command}")

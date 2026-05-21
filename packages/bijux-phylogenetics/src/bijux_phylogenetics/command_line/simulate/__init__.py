from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs

from .alignments import (
    add_simulate_alignment_commands,
    run_simulate_alignment_command,
)
from .continuous_traits import (
    add_simulate_continuous_trait_commands,
    run_simulate_continuous_trait_command,
)
from .correlated_traits import (
    add_simulate_correlated_trait_commands,
    run_simulate_correlated_trait_command,
)
from .discrete_histories import (
    add_simulate_discrete_history_commands,
    run_simulate_discrete_history_command,
)
from .discrete_traits import (
    add_simulate_discrete_trait_commands,
    run_simulate_discrete_trait_command,
)
from .reference_validation import (
    add_simulate_reference_validation_commands,
    run_simulate_reference_validation_command,
)
from .trees import add_simulate_tree_commands, run_simulate_tree_command


def add_simulate_command(subparsers: Any) -> None:
    simulate = subparsers.add_parser(
        get_command_spec("simulate").name, help=get_command_spec("simulate").summary
    )
    simulate_subparsers = simulate.add_subparsers(
        dest="simulate_command", required=True
    )
    add_simulate_tree_commands(simulate_subparsers)
    add_simulate_continuous_trait_commands(simulate_subparsers)
    add_simulate_correlated_trait_commands(simulate_subparsers)
    add_simulate_discrete_trait_commands(simulate_subparsers)
    add_simulate_discrete_history_commands(simulate_subparsers)
    add_simulate_alignment_commands(simulate_subparsers)
    add_simulate_reference_validation_commands(simulate_subparsers)


def run_simulate_command(args: Any, *, parser: Any) -> int:
    tree_exit_code = run_simulate_tree_command(args)
    if tree_exit_code is not None:
        return tree_exit_code
    continuous_trait_exit_code = run_simulate_continuous_trait_command(args)
    if continuous_trait_exit_code is not None:
        return continuous_trait_exit_code
    correlated_trait_exit_code = run_simulate_correlated_trait_command(
        args,
        parser=parser,
    )
    if correlated_trait_exit_code is not None:
        return correlated_trait_exit_code
    discrete_trait_exit_code = run_simulate_discrete_trait_command(args)
    if discrete_trait_exit_code is not None:
        return discrete_trait_exit_code
    discrete_history_exit_code = run_simulate_discrete_history_command(args)
    if discrete_history_exit_code is not None:
        return discrete_history_exit_code
    alignment_exit_code = run_simulate_alignment_command(args)
    if alignment_exit_code is not None:
        return alignment_exit_code
    reference_validation_exit_code = run_simulate_reference_validation_command(args)
    if reference_validation_exit_code is not None:
        return reference_validation_exit_code
    return 0

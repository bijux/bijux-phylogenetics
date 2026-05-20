from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec

from .modeling import add_modeling_commands, run_modeling_command
from .presentation import add_presentation_commands, run_presentation_command
from .stochastic_maps import add_stochastic_map_commands, run_stochastic_map_command
from .validation import add_validation_commands, run_validation_command


def add_discrete_evolution_commands(subparsers: Any) -> None:
    discrete_evolution = subparsers.add_parser(
        get_command_spec("discrete-evolution").name,
        help=get_command_spec("discrete-evolution").summary,
    )
    discrete_evolution_subparsers = discrete_evolution.add_subparsers(
        dest="discrete_evolution_command",
        required=True,
    )
    add_validation_commands(discrete_evolution_subparsers)
    add_modeling_commands(discrete_evolution_subparsers)
    add_stochastic_map_commands(discrete_evolution_subparsers)
    add_presentation_commands(discrete_evolution_subparsers)


def run_discrete_evolution_command(args: Any) -> int:
    validation_result = run_validation_command(args)
    if validation_result is not None:
        return validation_result

    modeling_result = run_modeling_command(args)
    if modeling_result is not None:
        return modeling_result

    stochastic_map_result = run_stochastic_map_command(args)
    if stochastic_map_result is not None:
        return stochastic_map_result

    presentation_result = run_presentation_command(args)
    if presentation_result is not None:
        return presentation_result

    raise ValueError(
        f"unsupported discrete-evolution command: {args.discrete_evolution_command}"
    )

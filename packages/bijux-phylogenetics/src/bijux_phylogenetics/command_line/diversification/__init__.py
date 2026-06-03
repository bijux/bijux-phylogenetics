from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec

from .clades import (
    add_diversification_clade_command,
    run_diversification_clade_command,
)
from .exclusions import (
    add_diversification_exclusion_commands,
    run_diversification_exclusion_command,
)
from .inspection import (
    add_diversification_inspection_commands,
    run_diversification_inspection_command,
)
from .modeling import (
    add_diversification_modeling_commands,
    run_diversification_modeling_command,
)
from .presentation import (
    add_diversification_presentation_commands,
    run_diversification_presentation_command,
)
from .trait_dependence import (
    add_diversification_trait_dependence_command,
    run_diversification_trait_dependence_command,
)


def add_diversification_commands(subparsers: Any) -> None:
    diversification = subparsers.add_parser(
        get_command_spec("diversification").name,
        help=get_command_spec("diversification").summary,
    )
    diversification_subparsers = diversification.add_subparsers(
        dest="diversification_command",
        required=True,
    )
    add_diversification_inspection_commands(diversification_subparsers)
    add_diversification_modeling_commands(diversification_subparsers)
    add_diversification_clade_command(diversification_subparsers)
    add_diversification_trait_dependence_command(diversification_subparsers)
    add_diversification_presentation_commands(diversification_subparsers)
    add_diversification_exclusion_commands(diversification_subparsers)


def run_diversification_command(args: Any) -> int:
    inspection_exit_code = run_diversification_inspection_command(args)
    if inspection_exit_code is not None:
        return inspection_exit_code

    modeling_exit_code = run_diversification_modeling_command(args)
    if modeling_exit_code is not None:
        return modeling_exit_code

    clade_exit_code = run_diversification_clade_command(args)
    if clade_exit_code is not None:
        return clade_exit_code

    trait_dependence_exit_code = run_diversification_trait_dependence_command(args)
    if trait_dependence_exit_code is not None:
        return trait_dependence_exit_code

    presentation_exit_code = run_diversification_presentation_command(args)
    if presentation_exit_code is not None:
        return presentation_exit_code

    exclusion_exit_code = run_diversification_exclusion_command(args)
    if exclusion_exit_code is not None:
        return exclusion_exit_code

    raise NotImplementedError(
        f"unsupported diversification command: {args.diversification_command}"
    )


__all__ = [
    "add_diversification_commands",
    "run_diversification_command",
]

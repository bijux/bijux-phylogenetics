from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.ancestral.discrete_diagnostics import (
    add_discrete_diagnostic_ancestral_commands,
    run_discrete_diagnostic_ancestral_command,
)
from bijux_phylogenetics.command_line.ancestral.presentation import (
    add_presentation_ancestral_commands,
    run_presentation_ancestral_command,
)
from bijux_phylogenetics.command_line.ancestral.reconstruction import (
    add_reconstruction_ancestral_commands,
    run_reconstruction_ancestral_command,
)
from bijux_phylogenetics.command_line.ancestral.sensitivity import (
    add_sensitivity_ancestral_commands,
    run_sensitivity_ancestral_command,
)
from bijux_phylogenetics.command_line.ancestral.stability import (
    add_stability_ancestral_commands,
    run_stability_ancestral_command,
)
from bijux_phylogenetics.command_line.registry import get_command_spec


def add_ancestral_commands(subparsers: Any) -> None:
    ancestral = subparsers.add_parser(
        get_command_spec("ancestral").name,
        help=get_command_spec("ancestral").summary,
    )
    ancestral_subparsers = ancestral.add_subparsers(
        dest="ancestral_command", required=True
    )
    add_reconstruction_ancestral_commands(ancestral_subparsers)
    add_stability_ancestral_commands(ancestral_subparsers)
    add_discrete_diagnostic_ancestral_commands(ancestral_subparsers)
    add_sensitivity_ancestral_commands(ancestral_subparsers)
    add_presentation_ancestral_commands(ancestral_subparsers)


def run_ancestral_command(args: Any, *, parser: Any) -> int:
    reconstruction_exit_code = run_reconstruction_ancestral_command(
        args,
        parser=parser,
    )
    if reconstruction_exit_code is not None:
        return reconstruction_exit_code
    stability_exit_code = run_stability_ancestral_command(args, parser=parser)
    if stability_exit_code is not None:
        return stability_exit_code
    discrete_diagnostic_exit_code = run_discrete_diagnostic_ancestral_command(
        args,
        parser=parser,
    )
    if discrete_diagnostic_exit_code is not None:
        return discrete_diagnostic_exit_code
    sensitivity_exit_code = run_sensitivity_ancestral_command(args, parser=parser)
    if sensitivity_exit_code is not None:
        return sensitivity_exit_code
    presentation_exit_code = run_presentation_ancestral_command(args, parser=parser)
    if presentation_exit_code is not None:
        return presentation_exit_code
    raise ValueError(f"unsupported ancestral command: {args.ancestral_command}")

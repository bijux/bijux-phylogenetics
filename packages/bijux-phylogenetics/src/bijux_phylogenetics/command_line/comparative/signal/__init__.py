from __future__ import annotations

from typing import Any

from .contrasts import add_signal_contrast_commands, run_signal_contrast_command
from .phylogenetic_signal import (
    add_phylogenetic_signal_commands,
    run_phylogenetic_signal_command,
)
from .readiness import add_signal_readiness_commands, run_signal_readiness_command


def add_comparative_signal_commands(comparative_subparsers: Any) -> None:
    add_signal_readiness_commands(comparative_subparsers)
    add_signal_contrast_commands(comparative_subparsers)
    add_phylogenetic_signal_commands(comparative_subparsers)


def run_comparative_signal_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    readiness_exit_code = run_signal_readiness_command(args)
    if readiness_exit_code is not None:
        return readiness_exit_code

    contrasts_exit_code = run_signal_contrast_command(args, parser=parser)
    if contrasts_exit_code is not None:
        return contrasts_exit_code

    return run_phylogenetic_signal_command(args)


__all__ = [
    "add_comparative_signal_commands",
    "run_comparative_signal_command",
]

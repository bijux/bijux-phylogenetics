from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec

from .assessment import (
    add_comparative_maturity_command,
    run_comparative_maturity_command,
)
from .clades import (
    add_comparative_clade_commands,
    run_comparative_clade_command,
)
from .evolution import (
    add_comparative_evolution_commands,
    run_comparative_evolution_command,
)
from .regression import (
    add_comparative_logistic_command,
    add_comparative_modeling_commands,
    add_comparative_pgls_commands,
    run_comparative_logistic_command,
    run_comparative_modeling_command,
    run_comparative_pgls_command,
)
from .signal import (
    add_comparative_signal_commands,
    run_comparative_signal_command,
)
from .traits import (
    add_comparative_trait_commands,
    run_comparative_trait_command,
)


def add_comparative_commands(subparsers: Any) -> None:
    comparative = subparsers.add_parser(
        get_command_spec("comparative").name,
        help=get_command_spec("comparative").summary,
    )
    comparative_subparsers = comparative.add_subparsers(
        dest="comparative_command",
        required=True,
    )
    add_comparative_signal_commands(comparative_subparsers)
    add_comparative_evolution_commands(comparative_subparsers)
    add_comparative_trait_commands(comparative_subparsers)
    add_comparative_maturity_command(comparative_subparsers)
    add_comparative_pgls_commands(comparative_subparsers)
    add_comparative_logistic_command(comparative_subparsers)
    add_comparative_clade_commands(comparative_subparsers)
    add_comparative_modeling_commands(comparative_subparsers)


def run_comparative_command(
    args: Any,
    *,
    parser: Any,
) -> int:
    signal_exit_code = run_comparative_signal_command(args, parser=parser)
    if signal_exit_code is not None:
        return signal_exit_code

    evolution_exit_code = run_comparative_evolution_command(args, parser=parser)
    if evolution_exit_code is not None:
        return evolution_exit_code

    review_exit_code = run_comparative_trait_command(args, parser=parser)
    if review_exit_code is not None:
        return review_exit_code

    maturity_exit_code = run_comparative_maturity_command(args)
    if maturity_exit_code is not None:
        return maturity_exit_code

    pgls_exit_code = run_comparative_pgls_command(args, parser=parser)
    if pgls_exit_code is not None:
        return pgls_exit_code

    logistic_exit_code = run_comparative_logistic_command(args)
    if logistic_exit_code is not None:
        return logistic_exit_code

    support_exit_code = run_comparative_clade_command(args, parser=parser)
    if support_exit_code is not None:
        return support_exit_code

    modeling_exit_code = run_comparative_modeling_command(args, parser=parser)
    if modeling_exit_code is not None:
        return modeling_exit_code

    raise NotImplementedError(
        f"unsupported comparative command: {args.comparative_command}"
    )


__all__ = [
    "add_comparative_commands",
    "add_comparative_evolution_commands",
    "add_comparative_signal_commands",
    "add_comparative_logistic_command",
    "add_comparative_maturity_command",
    "add_comparative_modeling_commands",
    "add_comparative_pgls_commands",
    "add_comparative_clade_commands",
    "add_comparative_trait_commands",
    "run_comparative_command",
    "run_comparative_evolution_command",
    "run_comparative_signal_command",
    "run_comparative_logistic_command",
    "run_comparative_maturity_command",
    "run_comparative_modeling_command",
    "run_comparative_pgls_command",
    "run_comparative_clade_command",
    "run_comparative_trait_command",
]

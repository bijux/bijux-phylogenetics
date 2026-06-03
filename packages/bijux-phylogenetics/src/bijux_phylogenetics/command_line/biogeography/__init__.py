from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.registry import get_command_spec

from .migration import (
    add_biogeography_migration_commands,
    run_biogeography_migration_command,
)
from .presentation import (
    add_biogeography_report_command,
    run_biogeography_report_command,
)
from .state_models import (
    add_biogeography_state_model_commands,
    run_biogeography_state_model_command,
)


def add_biogeography_commands(subparsers: Any) -> None:
    biogeography = subparsers.add_parser(
        get_command_spec("biogeography").name,
        help=get_command_spec("biogeography").summary,
    )
    biogeography_subparsers = biogeography.add_subparsers(
        dest="biogeography_command",
        required=True,
    )
    add_biogeography_state_model_commands(biogeography_subparsers)
    add_biogeography_migration_commands(biogeography_subparsers)
    add_biogeography_report_command(biogeography_subparsers)


def run_biogeography_command(args: Any) -> int:
    state_model_exit_code = run_biogeography_state_model_command(args)
    if state_model_exit_code is not None:
        return state_model_exit_code

    migration_exit_code = run_biogeography_migration_command(args)
    if migration_exit_code is not None:
        return migration_exit_code

    report_exit_code = run_biogeography_report_command(args)
    if report_exit_code is not None:
        return report_exit_code

    raise NotImplementedError(
        f"unsupported biogeography command: {args.biogeography_command}"
    )


__all__ = [
    "add_biogeography_commands",
    "run_biogeography_command",
]

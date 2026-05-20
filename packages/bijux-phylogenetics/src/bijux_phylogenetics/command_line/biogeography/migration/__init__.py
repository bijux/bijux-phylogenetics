from __future__ import annotations

from typing import Any

from .chronology import (
    add_biogeography_chronology_command,
    run_biogeography_chronology_command,
)
from .events import (
    add_biogeography_event_commands,
    run_biogeography_event_command,
)


def add_biogeography_migration_commands(biogeography_subparsers: Any) -> None:
    add_biogeography_chronology_command(biogeography_subparsers)
    add_biogeography_event_commands(biogeography_subparsers)


def run_biogeography_migration_command(args: Any) -> int | None:
    chronology_exit_code = run_biogeography_chronology_command(args)
    if chronology_exit_code is not None:
        return chronology_exit_code

    return run_biogeography_event_command(args)


__all__ = [
    "add_biogeography_migration_commands",
    "run_biogeography_migration_command",
]

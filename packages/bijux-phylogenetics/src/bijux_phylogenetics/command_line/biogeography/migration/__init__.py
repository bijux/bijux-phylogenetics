from __future__ import annotations

from typing import Any

from .chronology import (
    add_biogeography_chronology_command,
    run_biogeography_chronology_command,
)


def add_biogeography_migration_commands(biogeography_subparsers: Any) -> None:
    add_biogeography_chronology_command(biogeography_subparsers)


def run_biogeography_migration_command(args: Any) -> int | None:
    return run_biogeography_chronology_command(args)


__all__ = [
    "add_biogeography_migration_commands",
    "run_biogeography_migration_command",
]

from __future__ import annotations

from typing import Any

from .likelihood import (
    add_geographic_likelihood_command,
    run_geographic_likelihood_command,
)


def add_biogeography_state_model_commands(biogeography_subparsers: Any) -> None:
    add_geographic_likelihood_command(biogeography_subparsers)


def run_biogeography_state_model_command(args: Any) -> int | None:
    return run_geographic_likelihood_command(args)


__all__ = [
    "add_biogeography_state_model_commands",
    "run_biogeography_state_model_command",
]

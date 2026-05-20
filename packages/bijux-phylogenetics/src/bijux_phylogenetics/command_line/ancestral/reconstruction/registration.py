from __future__ import annotations

from typing import Any

from .comparison import (
    add_reconstruction_comparison_commands,
    run_reconstruction_comparison_command,
)
from .continuous import (
    add_continuous_reconstruction_commands,
    run_continuous_reconstruction_command,
)
from .discrete import (
    add_discrete_reconstruction_commands,
    run_discrete_reconstruction_command,
)
from .reference_validation import (
    add_reconstruction_reference_validation_commands,
    run_reconstruction_reference_validation_command,
)


def add_reconstruction_ancestral_commands(ancestral_subparsers: Any) -> None:
    add_continuous_reconstruction_commands(ancestral_subparsers)
    add_discrete_reconstruction_commands(ancestral_subparsers)
    add_reconstruction_reference_validation_commands(ancestral_subparsers)
    add_reconstruction_comparison_commands(ancestral_subparsers)


def run_reconstruction_ancestral_command(args: Any, *, parser: Any) -> int | None:
    continuous_result = run_continuous_reconstruction_command(args, parser=parser)
    if continuous_result is not None:
        return continuous_result
    discrete_result = run_discrete_reconstruction_command(args, parser=parser)
    if discrete_result is not None:
        return discrete_result
    reference_validation_result = run_reconstruction_reference_validation_command(args)
    if reference_validation_result is not None:
        return reference_validation_result
    comparison_result = run_reconstruction_comparison_command(args)
    if comparison_result is not None:
        return comparison_result

    return None

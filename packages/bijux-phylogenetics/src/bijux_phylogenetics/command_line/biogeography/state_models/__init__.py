from __future__ import annotations

from typing import Any

from .constrained import (
    add_constrained_geography_command,
    run_constrained_geography_command,
)
from .likelihood import (
    add_geographic_likelihood_command,
    run_geographic_likelihood_command,
)
from .sampling_bias import (
    add_sampling_bias_geography_command,
    run_sampling_bias_geography_command,
)
from .stratified import (
    add_time_stratified_geography_command,
    run_time_stratified_geography_command,
)


def add_biogeography_state_model_commands(biogeography_subparsers: Any) -> None:
    add_geographic_likelihood_command(biogeography_subparsers)
    add_constrained_geography_command(biogeography_subparsers)
    add_time_stratified_geography_command(biogeography_subparsers)
    add_sampling_bias_geography_command(biogeography_subparsers)


def run_biogeography_state_model_command(args: Any) -> int | None:
    likelihood_exit_code = run_geographic_likelihood_command(args)
    if likelihood_exit_code is not None:
        return likelihood_exit_code

    constrained_exit_code = run_constrained_geography_command(args)
    if constrained_exit_code is not None:
        return constrained_exit_code

    stratified_exit_code = run_time_stratified_geography_command(args)
    if stratified_exit_code is not None:
        return stratified_exit_code

    return run_sampling_bias_geography_command(args)


__all__ = [
    "add_biogeography_state_model_commands",
    "run_biogeography_state_model_command",
]

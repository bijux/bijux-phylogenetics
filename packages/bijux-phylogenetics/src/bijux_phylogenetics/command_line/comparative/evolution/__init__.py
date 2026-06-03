from __future__ import annotations

from typing import Any

from .continuous_models import (
    add_continuous_model_comparative_evolution_commands,
    run_continuous_model_comparative_evolution_command,
)
from .discrete_models import (
    add_discrete_model_comparative_evolution_commands,
    run_discrete_model_comparative_evolution_command,
)
from .disparity_workflows import (
    add_disparity_workflow_comparative_evolution_commands,
    run_disparity_workflow_comparative_evolution_command,
)
from .model_assessment import (
    add_model_assessment_comparative_evolution_commands,
    run_model_assessment_comparative_evolution_command,
)
from .regime_workflows import (
    add_regime_workflow_comparative_evolution_commands,
    run_regime_workflow_comparative_evolution_command,
)
from .trait_dependence import (
    add_trait_dependence_comparative_evolution_commands,
    run_trait_dependence_comparative_evolution_command,
)


def add_comparative_evolution_commands(comparative_subparsers: Any) -> None:
    add_discrete_model_comparative_evolution_commands(comparative_subparsers)
    add_trait_dependence_comparative_evolution_commands(comparative_subparsers)
    add_continuous_model_comparative_evolution_commands(comparative_subparsers)
    add_regime_workflow_comparative_evolution_commands(comparative_subparsers)
    add_disparity_workflow_comparative_evolution_commands(comparative_subparsers)
    add_model_assessment_comparative_evolution_commands(comparative_subparsers)


def run_comparative_evolution_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    discrete_model_exit_code = run_discrete_model_comparative_evolution_command(
        args,
        parser=parser,
    )
    if discrete_model_exit_code is not None:
        return discrete_model_exit_code

    trait_dependence_exit_code = run_trait_dependence_comparative_evolution_command(
        args,
        parser=parser,
    )
    if trait_dependence_exit_code is not None:
        return trait_dependence_exit_code

    continuous_model_exit_code = run_continuous_model_comparative_evolution_command(
        args,
        parser=parser,
    )
    if continuous_model_exit_code is not None:
        return continuous_model_exit_code

    regime_workflow_exit_code = run_regime_workflow_comparative_evolution_command(
        args,
        parser=parser,
    )
    if regime_workflow_exit_code is not None:
        return regime_workflow_exit_code

    disparity_workflow_exit_code = run_disparity_workflow_comparative_evolution_command(
        args,
        parser=parser,
    )
    if disparity_workflow_exit_code is not None:
        return disparity_workflow_exit_code

    return run_model_assessment_comparative_evolution_command(args, parser=parser)


__all__ = [
    "add_comparative_evolution_commands",
    "run_comparative_evolution_command",
]

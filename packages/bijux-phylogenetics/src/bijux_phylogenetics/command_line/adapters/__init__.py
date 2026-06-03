from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_external_adapter_execution_arguments,
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.errors import EngineUnavailableError
from bijux_phylogenetics.runtime.results import build_command_result

from .alignment_workflows import (
    add_adapter_alignment_workflow_commands,
    run_adapter_alignment_workflow_command,
)
from .bayesian import (
    add_bayesian_adapter_commands,
    run_bayesian_adapter_command,
)
from .beast import (
    add_beast_adapter_commands,
    run_beast_adapter_command,
)
from .fasta_to_tree import (
    add_adapter_fasta_to_tree_commands,
    run_adapter_fasta_to_tree_command,
)
from .inference import (
    add_inference_adapter_commands,
    run_inference_adapter_command,
)
from .inspection import (
    add_adapter_inspection_commands,
    run_adapter_inspection_command,
)
from .maximum_likelihood import (
    add_adapter_maximum_likelihood_commands,
    run_adapter_maximum_likelihood_command,
)
from .mrbayes import (
    add_mrbayes_adapter_commands,
    run_mrbayes_adapter_command,
)
from .reporting import (
    add_adapter_reporting_commands,
    run_adapter_reporting_command,
)
from .support_estimation import (
    add_adapter_support_estimation_commands,
    run_adapter_support_estimation_command,
)


def add_adapter_commands(subparsers: Any) -> None:
    adapter = subparsers.add_parser(
        get_command_spec("adapter").name, help=get_command_spec("adapter").summary
    )
    adapter_subparsers = adapter.add_subparsers(dest="adapter_command", required=True)
    add_adapter_inspection_commands(adapter_subparsers)
    add_adapter_reporting_commands(adapter_subparsers)

    add_adapter_alignment_workflow_commands(adapter_subparsers)
    add_adapter_maximum_likelihood_commands(adapter_subparsers)
    add_adapter_support_estimation_commands(adapter_subparsers)
    add_adapter_fasta_to_tree_commands(adapter_subparsers)

    add_inference_adapter_commands(adapter_subparsers)
    add_mrbayes_adapter_commands(adapter_subparsers)
    add_beast_adapter_commands(adapter_subparsers)
    add_bayesian_adapter_commands(adapter_subparsers)


def run_adapter_command(args: Any) -> int | None:
    if args.command != "adapter":
        return None
    inspection_exit_code = run_adapter_inspection_command(args)
    if inspection_exit_code is not None:
        return inspection_exit_code
    reporting_exit_code = run_adapter_reporting_command(args)
    if reporting_exit_code is not None:
        return reporting_exit_code
    alignment_exit_code = run_adapter_alignment_workflow_command(args)
    if alignment_exit_code is not None:
        return alignment_exit_code
    maximum_likelihood_exit_code = run_adapter_maximum_likelihood_command(args)
    if maximum_likelihood_exit_code is not None:
        return maximum_likelihood_exit_code
    support_estimation_exit_code = run_adapter_support_estimation_command(args)
    if support_estimation_exit_code is not None:
        return support_estimation_exit_code
    fasta_to_tree_exit_code = run_adapter_fasta_to_tree_command(args)
    if fasta_to_tree_exit_code is not None:
        return fasta_to_tree_exit_code
    mrbayes_exit_code = run_mrbayes_adapter_command(args)
    if mrbayes_exit_code is not None:
        return mrbayes_exit_code
    beast_exit_code = run_beast_adapter_command(args)
    if beast_exit_code is not None:
        return beast_exit_code
    bayesian_exit_code = run_bayesian_adapter_command(args)
    if bayesian_exit_code is not None:
        return bayesian_exit_code
    inference_exit_code = run_inference_adapter_command(args)
    if inference_exit_code is not None:
        return inference_exit_code
    raise EngineUnavailableError(f"unsupported adapter command: {args.adapter_command}")

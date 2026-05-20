from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.discrete_evolution import (
    detect_state_imbalance_problems,
    validate_discrete_state_coding,
    validate_discrete_transition_reference_examples,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import COMMAND_NAME, allowed_states, model_inputs, ordered_states


def add_validation_commands(discrete_evolution_subparsers: Any) -> None:
    discrete_validate = discrete_evolution_subparsers.add_parser(
        "validate-coding",
        help="Validate discrete-state labels against tree-overlapping taxa.",
    )
    discrete_validate.add_argument("tree", type=Path)
    discrete_validate.add_argument("table", type=Path)
    discrete_validate.add_argument("--trait", required=True)
    discrete_validate.add_argument("--taxon-column")
    discrete_validate.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, accept any single token state label.",
    )
    discrete_validate.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_validate.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_validate.add_argument(
        "--json", action="store_true", help="Emit the validation report as JSON."
    )
    _add_manifest_argument(discrete_validate)

    discrete_imbalance = discrete_evolution_subparsers.add_parser(
        "imbalance",
        help="Detect rare, dominant, or degenerate state balance problems.",
    )
    discrete_imbalance.add_argument("tree", type=Path)
    discrete_imbalance.add_argument("table", type=Path)
    discrete_imbalance.add_argument("--trait", required=True)
    discrete_imbalance.add_argument("--taxon-column")
    discrete_imbalance.add_argument(
        "--json", action="store_true", help="Emit the imbalance report as JSON."
    )
    _add_manifest_argument(discrete_imbalance)

    discrete_reference = discrete_evolution_subparsers.add_parser(
        "reference",
        help="Validate deterministic discrete-state transition examples against built-in reference expectations.",
    )
    discrete_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the reference-validation report as JSON.",
    )
    _add_manifest_argument(discrete_reference)


def run_validation_command(args: Any) -> int | None:
    if args.discrete_evolution_command == "validate-coding":
        return _run_validate_coding(args)
    if args.discrete_evolution_command == "imbalance":
        return _run_imbalance(args)
    if args.discrete_evolution_command == "reference":
        return _run_reference(args)
    return None


def _run_validate_coding(args: Any) -> int:
    report = validate_discrete_state_coding(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
    )
    outputs = _finalize_outputs(
        args,
        command=COMMAND_NAME,
        inputs=model_inputs(args),
    )
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=model_inputs(args),
            outputs=outputs,
            metrics={
                "valid": report.valid,
                "issue_count": len(report.issues),
                "observed_state_count": len(report.observed_states),
                "state_ordering": report.state_ordering,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_imbalance(args: Any) -> int:
    report = detect_state_imbalance_problems(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
    )
    outputs = _finalize_outputs(
        args,
        command=COMMAND_NAME,
        inputs=model_inputs(args),
    )
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=model_inputs(args),
            outputs=outputs,
            warnings=[warning.message for warning in report.warnings],
            metrics={
                "taxon_count": report.taxon_count,
                "observed_state_count": len(report.observed_states),
                "warning_count": len(report.warnings),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_reference(args: Any) -> int:
    report = validate_discrete_transition_reference_examples()
    outputs = _finalize_outputs(args, command=COMMAND_NAME, inputs=[])
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=[],
            outputs=outputs,
            metrics={
                "case_count": report.case_count,
                "all_passed": report.all_passed,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

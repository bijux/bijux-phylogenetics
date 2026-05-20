from __future__ import annotations

from typing import Any

from bijux_phylogenetics.ancestral.discrete.review import (
    validate_discrete_ancestral_reference_examples,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_reconstruction_reference_validation_commands(
    ancestral_subparsers: Any,
) -> None:
    ancestral_discrete_reference = ancestral_subparsers.add_parser(
        "discrete-reference",
        help="Validate built-in discrete ancestral reference examples.",
    )
    ancestral_discrete_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the reference validation report as JSON.",
    )


def run_reconstruction_reference_validation_command(args: Any) -> int | None:
    if args.ancestral_command != "discrete-reference":
        return None

    report = validate_discrete_ancestral_reference_examples()
    outputs = _finalize_outputs(args, command="ancestral", inputs=[])
    _print_result(
        build_command_result(
            command="ancestral",
            inputs=[],
            outputs=outputs,
            metrics={
                "case_count": report.case_count,
                "external_case_count": report.external_case_count,
                "all_passed": report.all_passed,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

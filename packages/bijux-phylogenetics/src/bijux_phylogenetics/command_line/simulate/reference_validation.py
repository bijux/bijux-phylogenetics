from __future__ import annotations

from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    validate_geiger_sim_char_reference_examples,
)


def add_simulate_reference_validation_commands(simulate_subparsers: Any) -> None:
    simulate_validate_sim_char_reference = simulate_subparsers.add_parser(
        "validate-sim-char-reference",
        help="Validate governed geiger::sim.char summary envelopes.",
    )
    simulate_validate_sim_char_reference.add_argument(
        "--json",
        action="store_true",
        help="Emit the governed validation report as JSON.",
    )
    _add_manifest_argument(simulate_validate_sim_char_reference)


def run_simulate_reference_validation_command(args: Any) -> int | None:
    if args.simulate_command != "validate-sim-char-reference":
        return None

    report = validate_geiger_sim_char_reference_examples()
    _print_result(
        build_command_result(
            command="simulate",
            inputs=[],
            outputs=[],
            metrics={
                "case_count": report.case_count,
                "all_passed": report.all_passed,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

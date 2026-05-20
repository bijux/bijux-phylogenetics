from __future__ import annotations

from pathlib import Path
from typing import Any

from .continuous_traits import (
    add_simulate_continuous_trait_commands,
    run_simulate_continuous_trait_command,
)
from .correlated_traits import (
    add_simulate_correlated_trait_commands,
    run_simulate_correlated_trait_command,
)
from .discrete_histories import (
    add_simulate_discrete_history_commands,
    run_simulate_discrete_history_command,
)
from .discrete_traits import (
    add_simulate_discrete_trait_commands,
    run_simulate_discrete_trait_command,
)
from .trees import add_simulate_tree_commands, run_simulate_tree_command
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_dna_alignment,
    simulate_protein_alignment,
    validate_geiger_sim_char_reference_examples,
    write_simulated_alignment,
)


def add_simulate_command(subparsers: Any) -> None:
    simulate = subparsers.add_parser(
        get_command_spec("simulate").name, help=get_command_spec("simulate").summary
    )
    simulate_subparsers = simulate.add_subparsers(
        dest="simulate_command", required=True
    )
    add_simulate_tree_commands(simulate_subparsers)
    add_simulate_continuous_trait_commands(simulate_subparsers)
    add_simulate_correlated_trait_commands(simulate_subparsers)
    add_simulate_discrete_trait_commands(simulate_subparsers)
    add_simulate_discrete_history_commands(simulate_subparsers)

    simulate_dna = simulate_subparsers.add_parser(
        "alignment-dna",
        help="Simulate a DNA alignment along a rooted tree.",
    )
    simulate_dna.add_argument("tree", type=Path)
    simulate_dna.add_argument("--sequence-length", type=int, required=True)
    simulate_dna.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_dna.add_argument("--seed", type=int, default=1)
    simulate_dna.add_argument("--out", required=True, type=Path)
    simulate_dna.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_dna)

    simulate_protein = simulate_subparsers.add_parser(
        "alignment-protein",
        help="Simulate a protein alignment along a rooted tree.",
    )
    simulate_protein.add_argument("tree", type=Path)
    simulate_protein.add_argument("--sequence-length", type=int, required=True)
    simulate_protein.add_argument("--substitution-rate", type=float, default=1.0)
    simulate_protein.add_argument("--seed", type=int, default=1)
    simulate_protein.add_argument("--out", required=True, type=Path)
    simulate_protein.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_protein)

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


def run_simulate_command(args: Any, *, parser: Any) -> int:
    tree_exit_code = run_simulate_tree_command(args)
    if tree_exit_code is not None:
        return tree_exit_code
    continuous_trait_exit_code = run_simulate_continuous_trait_command(args)
    if continuous_trait_exit_code is not None:
        return continuous_trait_exit_code
    correlated_trait_exit_code = run_simulate_correlated_trait_command(
        args,
        parser=parser,
    )
    if correlated_trait_exit_code is not None:
        return correlated_trait_exit_code
    discrete_trait_exit_code = run_simulate_discrete_trait_command(args)
    if discrete_trait_exit_code is not None:
        return discrete_trait_exit_code
    discrete_history_exit_code = run_simulate_discrete_history_command(args)
    if discrete_history_exit_code is not None:
        return discrete_history_exit_code
    if args.simulate_command == "alignment-dna":
        report = simulate_dna_alignment(
            args.tree,
            sequence_length=args.sequence_length,
            substitution_rate=args.substitution_rate,
            seed=args.seed,
        )
        output_path = write_simulated_alignment(args.out, report)
        outputs = _finalize_outputs(
            args, command="simulate", inputs=[args.tree], outputs=[output_path]
        )
        _print_result(
            build_command_result(
                command="simulate",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tip_count": report.tip_count,
                    "sequence_length": report.sequence_length,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.simulate_command == "validate-sim-char-reference":
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
    report = simulate_protein_alignment(
        args.tree,
        sequence_length=args.sequence_length,
        substitution_rate=args.substitution_rate,
        seed=args.seed,
    )
    output_path = write_simulated_alignment(args.out, report)
    outputs = _finalize_outputs(
        args, command="simulate", inputs=[args.tree], outputs=[output_path]
    )
    _print_result(
        build_command_result(
            command="simulate",
            inputs=[args.tree],
            outputs=outputs,
            metrics={
                "tip_count": report.tip_count,
                "sequence_length": report.sequence_length,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

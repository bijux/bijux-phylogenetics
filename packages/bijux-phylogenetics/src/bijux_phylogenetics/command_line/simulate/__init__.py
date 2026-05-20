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
from .trees import add_simulate_tree_commands, run_simulate_tree_command
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_probability_assignments,
    _parse_rate_rows,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_discrete_histories,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_protein_alignment,
    validate_geiger_sim_char_reference_examples,
    write_discrete_history_branch_truth_table,
    write_discrete_history_event_table,
    write_discrete_history_node_truth_table,
    write_discrete_history_segment_table,
    write_discrete_history_summary_table,
    write_discrete_history_tip_truth_table,
    write_discrete_trait_table,
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

    simulate_discrete = simulate_subparsers.add_parser(
        "traits-discrete",
        help="Simulate a discrete tip trait under a symmetric jump process.",
    )
    simulate_discrete.add_argument("tree", type=Path)
    simulate_discrete.add_argument("--states", nargs="+", required=True)
    simulate_discrete.add_argument("--transition-rate", type=float, default=1.0)
    simulate_discrete.add_argument("--root-state")
    simulate_discrete.add_argument("--seed", type=int, default=1)
    simulate_discrete.add_argument("--out", required=True, type=Path)
    simulate_discrete.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_discrete)

    simulate_history_discrete = simulate_subparsers.add_parser(
        "history-discrete",
        help="Simulate true discrete histories on a fixed tree from one explicit rate matrix.",
    )
    simulate_history_discrete.add_argument("tree", type=Path)
    simulate_history_discrete.add_argument("--states", nargs="+", required=True)
    simulate_history_discrete.add_argument(
        "--rate",
        action="append",
        required=True,
        help="One SOURCE->TARGET=RATE entry. Repeat to build the rate matrix.",
    )
    simulate_history_discrete.add_argument("--root-state")
    simulate_history_discrete.add_argument(
        "--root-probability",
        action="append",
        default=[],
        help="One STATE=PROBABILITY entry. Repeat to define the root prior.",
    )
    simulate_history_discrete.add_argument(
        "--transform",
        choices=("lambda", "kappa", "delta", "early-burst"),
        help="Optional discrete branch-length transform applied before history simulation.",
    )
    simulate_history_discrete.add_argument(
        "--transform-parameter-value",
        type=float,
        help="Transform parameter value used for the discrete branch-length transform.",
    )
    simulate_history_discrete.add_argument("--replicates", type=int, default=1)
    simulate_history_discrete.add_argument("--seed", type=int, default=1)
    simulate_history_discrete.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Tip-state truth table output path.",
    )
    simulate_history_discrete.add_argument("--nodes-out", type=Path)
    simulate_history_discrete.add_argument("--branches-out", type=Path)
    simulate_history_discrete.add_argument("--events-out", type=Path)
    simulate_history_discrete.add_argument("--segments-out", type=Path)
    simulate_history_discrete.add_argument("--summary-out", type=Path)
    simulate_history_discrete.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_history_discrete)

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
    if args.simulate_command == "traits-discrete":
        report = simulate_discrete_traits(
            args.tree,
            states=args.states,
            transition_rate=args.transition_rate,
            root_state=args.root_state,
            seed=args.seed,
        )
        output_path = write_discrete_trait_table(args.out, report)
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
                    "trait_count": len(report.traits),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.simulate_command == "history-discrete":
        root_probability_rows = _parse_probability_assignments(args.root_probability)
        report = simulate_discrete_histories(
            args.tree,
            states=args.states,
            rate_rows=_parse_rate_rows(args.rate),
            root_state=args.root_state,
            root_state_probabilities=(root_probability_rows or None),
            transform=args.transform,
            transform_parameter_value=args.transform_parameter_value,
            replicates=args.replicates,
            seed=args.seed,
        )
        outputs_to_finalize = [write_discrete_history_tip_truth_table(args.out, report)]
        if args.nodes_out is not None:
            outputs_to_finalize.append(
                write_discrete_history_node_truth_table(args.nodes_out, report)
            )
        if args.branches_out is not None:
            outputs_to_finalize.append(
                write_discrete_history_branch_truth_table(args.branches_out, report)
            )
        if args.events_out is not None:
            outputs_to_finalize.append(
                write_discrete_history_event_table(args.events_out, report)
            )
        if args.segments_out is not None:
            outputs_to_finalize.append(
                write_discrete_history_segment_table(args.segments_out, report)
            )
        if args.summary_out is not None:
            outputs_to_finalize.append(
                write_discrete_history_summary_table(args.summary_out, report)
            )
        outputs = _finalize_outputs(
            args,
            command="simulate",
            inputs=[args.tree],
            outputs=outputs_to_finalize,
        )
        _print_result(
            build_command_result(
                command="simulate",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tip_count": report.tip_count,
                    "branch_count": report.branch_count,
                    "replicate_count": report.replicate_count,
                    "state_count": len(report.states),
                    "mean_total_transition_count": report.mean_total_transition_count,
                    "transform_name": report.transform_name,
                    "transform_parameter_value": report.transform_parameter_value,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
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

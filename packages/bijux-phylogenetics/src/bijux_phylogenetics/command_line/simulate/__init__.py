from __future__ import annotations

from pathlib import Path
from typing import Any

from .trees import add_simulate_tree_commands, run_simulate_tree_command
from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_float_csv_row,
    _parse_probability_assignments,
    _parse_rate_rows,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_brownian_traits,
    simulate_correlated_brownian_trait_collection,
    simulate_discrete_histories,
    simulate_discrete_traits,
    simulate_dna_alignment,
    simulate_early_burst_traits,
    simulate_ou_traits,
    simulate_protein_alignment,
    simulate_speciational_traits,
    validate_geiger_sim_char_reference_examples,
    write_continuous_trait_table,
    write_correlated_continuous_trait_collection_summary_table,
    write_correlated_continuous_trait_collection_table,
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

    simulate_brownian = simulate_subparsers.add_parser(
        "traits-brownian",
        help="Simulate a continuous tip trait under Brownian motion.",
    )
    simulate_brownian.add_argument("tree", type=Path)
    simulate_brownian.add_argument("--root-state", type=float, default=0.0)
    simulate_brownian.add_argument("--sigma", type=float)
    simulate_brownian.add_argument("--sigma-squared", type=float)
    simulate_brownian.add_argument("--seed", type=int, default=1)
    simulate_brownian.add_argument("--out", required=True, type=Path)
    simulate_brownian.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_brownian)

    simulate_speciational = simulate_subparsers.add_parser(
        "traits-speciational",
        help="Simulate a continuous tip trait under geiger-style speciational Brownian motion.",
    )
    simulate_speciational.add_argument("tree", type=Path)
    simulate_speciational.add_argument("--root-state", type=float, default=0.0)
    simulate_speciational.add_argument("--sigma", type=float)
    simulate_speciational.add_argument("--sigma-squared", type=float)
    simulate_speciational.add_argument("--seed", type=int, default=1)
    simulate_speciational.add_argument("--out", required=True, type=Path)
    simulate_speciational.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_speciational)

    simulate_brownian_correlated = simulate_subparsers.add_parser(
        "traits-brownian-correlated",
        help="Simulate correlated continuous tip traits under multivariate Brownian motion.",
    )
    simulate_brownian_correlated.add_argument("tree", type=Path)
    simulate_brownian_correlated.add_argument(
        "--trait",
        action="append",
        required=True,
        help="One trait name. Repeat to define the multivariate trait order.",
    )
    simulate_brownian_correlated.add_argument(
        "--root-state",
        action="append",
        type=float,
        default=[],
        help="One root state value aligned to the declared trait order.",
    )
    covariance_group = simulate_brownian_correlated.add_mutually_exclusive_group(
        required=True
    )
    covariance_group.add_argument(
        "--covariance-row",
        action="append",
        dest="covariance_rows",
        help="One comma-separated covariance-matrix row. Repeat to build the full matrix.",
    )
    covariance_group.add_argument(
        "--correlation-row",
        action="append",
        dest="correlation_rows",
        help="One comma-separated correlation-matrix row. Repeat to build the full matrix.",
    )
    simulate_brownian_correlated.add_argument(
        "--trait-standard-deviation",
        action="append",
        type=float,
        default=[],
        help="One trait standard deviation aligned to the declared trait order when using --correlation-row.",
    )
    simulate_brownian_correlated.add_argument("--replicates", type=int, default=128)
    simulate_brownian_correlated.add_argument("--seed", type=int, default=1)
    simulate_brownian_correlated.add_argument("--out", required=True, type=Path)
    simulate_brownian_correlated.add_argument("--summary-out", type=Path)
    simulate_brownian_correlated.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_brownian_correlated)

    simulate_ou = simulate_subparsers.add_parser(
        "traits-ou",
        help="Simulate a continuous tip trait under an OU process.",
    )
    simulate_ou.add_argument("tree", type=Path)
    simulate_ou.add_argument("--root-state", type=float, default=0.0)
    simulate_ou.add_argument("--sigma", type=float, default=1.0)
    simulate_ou.add_argument("--alpha", type=float, default=1.0)
    simulate_ou.add_argument("--theta", type=float, default=0.0)
    simulate_ou.add_argument("--seed", type=int, default=1)
    simulate_ou.add_argument("--out", required=True, type=Path)
    simulate_ou.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_ou)

    simulate_early_burst = simulate_subparsers.add_parser(
        "traits-early-burst",
        help="Simulate a continuous tip trait under an early-burst branch-rate process.",
    )
    simulate_early_burst.add_argument("tree", type=Path)
    simulate_early_burst.add_argument("--root-state", type=float, default=0.0)
    simulate_early_burst.add_argument("--sigma", type=float, default=1.0)
    simulate_early_burst.add_argument("--rate-change", type=float, default=1.0)
    simulate_early_burst.add_argument("--seed", type=int, default=1)
    simulate_early_burst.add_argument("--out", required=True, type=Path)
    simulate_early_burst.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_early_burst)

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
    if args.simulate_command == "traits-brownian":
        report = simulate_brownian_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            sigma_squared=args.sigma_squared,
            seed=args.seed,
        )
        output_path = write_continuous_trait_table(args.out, report)
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
                    "sigma_squared": report.sigma_squared,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.simulate_command == "traits-speciational":
        report = simulate_speciational_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            sigma_squared=args.sigma_squared,
            seed=args.seed,
        )
        output_path = write_continuous_trait_table(args.out, report)
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
                    "sigma_squared": report.sigma_squared,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.simulate_command == "traits-brownian-correlated":
        if args.covariance_rows is None and args.correlation_rows is None:
            parser.error(
                "correlated Brownian simulation requires either --covariance-row or --correlation-row"
            )
        if args.covariance_rows is not None and args.trait_standard_deviation:
            parser.error(
                "trait standard deviations can only be used with --correlation-row"
            )
        if args.correlation_rows is not None and not args.trait_standard_deviation:
            parser.error(
                "correlated Brownian simulation requires --trait-standard-deviation with --correlation-row"
            )
        report = simulate_correlated_brownian_trait_collection(
            args.tree,
            trait_names=args.trait,
            evolutionary_covariance_matrix=(
                None
                if args.covariance_rows is None
                else [_parse_float_csv_row(raw_row) for raw_row in args.covariance_rows]
            ),
            evolutionary_correlation_matrix=(
                None
                if args.correlation_rows is None
                else [_parse_float_csv_row(raw_row) for raw_row in args.correlation_rows]
            ),
            trait_standard_deviations=(
                None if not args.trait_standard_deviation else args.trait_standard_deviation
            ),
            root_states=None if not args.root_state else args.root_state,
            replicates=args.replicates,
            seed=args.seed,
        )
        outputs: list[Path | str] = [
            write_correlated_continuous_trait_collection_table(args.out, report)
        ]
        if args.summary_out is not None:
            outputs.append(
                write_correlated_continuous_trait_collection_summary_table(
                    args.summary_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="simulate",
            inputs=[args.tree],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="simulate",
                inputs=[args.tree],
                outputs=outputs,
                metrics={
                    "tip_count": report.tip_count,
                    "trait_count": len(report.trait_names),
                    "replicate_count": report.replicate_count,
                    "summary_row_count": len(report.rows),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.simulate_command == "traits-ou":
        report = simulate_ou_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            alpha=args.alpha,
            theta=args.theta,
            seed=args.seed,
        )
        output_path = write_continuous_trait_table(args.out, report)
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
    if args.simulate_command == "traits-early-burst":
        report = simulate_early_burst_traits(
            args.tree,
            root_state=args.root_state,
            sigma=args.sigma,
            rate_change=args.rate_change,
            seed=args.seed,
        )
        output_path = write_continuous_trait_table(args.out, report)
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
                    "rate_change": report.rate_change,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
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

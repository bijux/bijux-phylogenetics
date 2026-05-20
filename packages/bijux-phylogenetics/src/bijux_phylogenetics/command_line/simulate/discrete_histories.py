from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_manifest_argument,
    _parse_probability_assignments,
    _parse_rate_rows,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_discrete_histories,
    write_discrete_history_branch_truth_table,
    write_discrete_history_event_table,
    write_discrete_history_node_truth_table,
    write_discrete_history_segment_table,
    write_discrete_history_summary_table,
    write_discrete_history_tip_truth_table,
)


def add_simulate_discrete_history_commands(simulate_subparsers: Any) -> None:
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


def run_simulate_discrete_history_command(args: Any) -> int | None:
    if args.simulate_command != "history-discrete":
        return None

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

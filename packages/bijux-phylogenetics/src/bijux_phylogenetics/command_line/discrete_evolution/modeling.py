from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.discrete_evolution import (
    compare_discrete_state_models,
    estimate_ancestral_geographic_states,
    write_discrete_model_comparison_table,
    write_node_state_probability_table,
    write_transition_summary_table,
)
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import COMMAND_NAME, allowed_states, model_inputs, ordered_states


def add_modeling_commands(discrete_evolution_subparsers: Any) -> None:
    discrete_model = discrete_evolution_subparsers.add_parser(
        "model",
        help="Run one discrete-state transition model and export node or branch summaries.",
    )
    discrete_model.add_argument("tree", type=Path)
    discrete_model.add_argument("table", type=Path)
    discrete_model.add_argument("--trait", required=True)
    discrete_model.add_argument("--taxon-column")
    discrete_model.add_argument(
        "--model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_model.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_model.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_model.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_model.add_argument(
        "--node-table-out", type=Path, help="Write node-state probabilities as TSV."
    )
    discrete_model.add_argument(
        "--transitions-out", type=Path, help="Write branch transition summaries as TSV."
    )
    discrete_model.add_argument(
        "--json", action="store_true", help="Emit the model report as JSON."
    )
    _add_manifest_argument(discrete_model)

    discrete_compare = discrete_evolution_subparsers.add_parser(
        "compare-models",
        help="Compare two supported discrete-state evolution models node by node.",
    )
    discrete_compare.add_argument("tree", type=Path)
    discrete_compare.add_argument("table", type=Path)
    discrete_compare.add_argument("--trait", required=True)
    discrete_compare.add_argument("--taxon-column")
    discrete_compare.add_argument(
        "--left-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="equal-rates",
    )
    discrete_compare.add_argument(
        "--right-model",
        choices=("equal-rates", "symmetric", "all-rates-different", "meristic"),
        default="all-rates-different",
    )
    discrete_compare.add_argument(
        "--state-ordering", choices=("unordered", "ordered"), default="unordered"
    )
    discrete_compare.add_argument(
        "--ordered-states", help="Comma-delimited explicit ordered state vocabulary."
    )
    discrete_compare.add_argument(
        "--allowed-states",
        help="Comma-delimited allowed state vocabulary. When omitted, infer observed states from the table.",
    )
    discrete_compare.add_argument(
        "--table-out", type=Path, help="Write node-wise model differences as TSV."
    )
    discrete_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(discrete_compare)


def run_modeling_command(args: Any) -> int | None:
    if args.discrete_evolution_command == "model":
        return _run_model(args)
    if args.discrete_evolution_command == "compare-models":
        return _run_compare_models(args)
    return None


def _run_model(args: Any) -> int:
    report = estimate_ancestral_geographic_states(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        model=args.model,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
    )
    outputs: list[Path | str] = []
    if args.node_table_out is not None:
        outputs.append(write_node_state_probability_table(args.node_table_out, report))
    if args.transitions_out is not None:
        outputs.append(write_transition_summary_table(args.transitions_out, report))
    outputs = _finalize_outputs(
        args,
        command=COMMAND_NAME,
        inputs=model_inputs(args),
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=model_inputs(args),
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "taxon_count": report.taxon_count,
                "observed_state_count": len(report.observed_states),
                "transition_count": report.transition_summary.transition_count,
                "strongly_supported_transition_count": (
                    report.transition_summary.strongly_supported_transition_count
                ),
                "model": report.model,
                "state_ordering": report.state_ordering,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0


def _run_compare_models(args: Any) -> int:
    comparison = compare_discrete_state_models(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
        left_model=args.left_model,
        right_model=args.right_model,
        allowed_states=allowed_states(args),
        state_ordering=args.state_ordering,
        ordered_states=ordered_states(args),
    )
    outputs: list[Path | str] = []
    if args.table_out is not None:
        outputs.append(
            write_discrete_model_comparison_table(args.table_out, comparison)
        )
    outputs = _finalize_outputs(
        args,
        command=COMMAND_NAME,
        inputs=model_inputs(args),
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command=COMMAND_NAME,
            inputs=model_inputs(args),
            outputs=outputs,
            metrics={
                "better_model": comparison.better_model,
                "model_count": len(comparison.rows),
                "differing_node_count": sum(
                    1 for row in comparison.node_differences if row.differs
                ),
                "state_ordering": args.state_ordering,
            },
            data=comparison,
        ),
        json_output=args.json,
    )
    return 0

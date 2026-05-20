from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.simulation import (
    simulate_birth_death_trees,
    simulate_coalescent_trees,
    simulate_random_trees,
    write_tree_set,
    write_tree_simulation_envelope_table,
    write_tree_simulation_record_table,
)


def add_simulate_tree_commands(simulate_subparsers: Any) -> None:
    simulate_birth_death = simulate_subparsers.add_parser(
        "tree-birth-death",
        help="Simulate one or more trees under a birth-death process.",
    )
    simulate_birth_death.add_argument("--tree-count", type=int, default=1)
    simulate_birth_death.add_argument("--tip-count", type=int, required=True)
    simulate_birth_death.add_argument("--birth-rate", type=float, default=1.0)
    simulate_birth_death.add_argument("--death-rate", type=float, default=0.25)
    simulate_birth_death.add_argument("--seed", type=int, default=1)
    simulate_birth_death.add_argument("--out", required=True, type=Path)
    simulate_birth_death.add_argument("--record-table-out", type=Path)
    simulate_birth_death.add_argument("--envelope-table-out", type=Path)
    simulate_birth_death.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_birth_death)

    simulate_random_tree = simulate_subparsers.add_parser(
        "tree-random",
        help="Simulate one or more rooted random trees with uniform branch lengths.",
    )
    simulate_random_tree.add_argument("--tree-count", type=int, default=1)
    simulate_random_tree.add_argument("--tip-count", type=int, required=True)
    simulate_random_tree.add_argument("--seed", type=int, default=1)
    simulate_random_tree.add_argument("--out", required=True, type=Path)
    simulate_random_tree.add_argument("--record-table-out", type=Path)
    simulate_random_tree.add_argument("--envelope-table-out", type=Path)
    simulate_random_tree.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_random_tree)

    simulate_coalescent = simulate_subparsers.add_parser(
        "tree-coalescent",
        help="Simulate one or more trees under a coalescent model.",
    )
    simulate_coalescent.add_argument("--tree-count", type=int, default=1)
    simulate_coalescent.add_argument("--tip-count", type=int, required=True)
    simulate_coalescent.add_argument("--population-size", type=float, default=1.0)
    simulate_coalescent.add_argument("--seed", type=int, default=1)
    simulate_coalescent.add_argument("--out", required=True, type=Path)
    simulate_coalescent.add_argument("--record-table-out", type=Path)
    simulate_coalescent.add_argument("--envelope-table-out", type=Path)
    simulate_coalescent.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_coalescent)


def run_simulate_tree_command(args: Any) -> int | None:
    if args.simulate_command == "tree-birth-death":
        trees, report = simulate_birth_death_trees(
            tree_count=args.tree_count,
            tip_count=args.tip_count,
            birth_rate=args.birth_rate,
            death_rate=args.death_rate,
            seed=args.seed,
        )
        _print_tree_result(args, trees=trees, report=report)
        return 0
    if args.simulate_command == "tree-random":
        trees, report = simulate_random_trees(
            tree_count=args.tree_count,
            tip_count=args.tip_count,
            seed=args.seed,
        )
        _print_tree_result(
            args,
            trees=trees,
            report=report,
            extra_metrics={"branch_length_model": report.branch_length_model},
        )
        return 0
    if args.simulate_command != "tree-coalescent":
        return None

    trees, report = simulate_coalescent_trees(
        tree_count=args.tree_count,
        tip_count=args.tip_count,
        population_size=args.population_size,
        seed=args.seed,
    )
    _print_tree_result(args, trees=trees, report=report)
    return 0


def _print_tree_result(
    args: Any,
    *,
    trees: list[Any],
    report: Any,
    extra_metrics: dict[str, Any] | None = None,
) -> None:
    output_path = write_tree_set(args.out, trees)
    outputs_to_finalize = [output_path]
    if args.record_table_out is not None:
        outputs_to_finalize.append(
            write_tree_simulation_record_table(args.record_table_out, report)
        )
    if args.envelope_table_out is not None:
        outputs_to_finalize.append(
            write_tree_simulation_envelope_table(args.envelope_table_out, report)
        )
    outputs = _finalize_outputs(
        args,
        command="simulate",
        inputs=[],
        outputs=outputs_to_finalize,
    )
    metrics = {
        "tree_count": report.tree_count,
        "tip_count": report.tip_count,
        "pooled_branch_count": report.pooled_branch_count,
        "envelope_metric_count": len(report.envelope_metrics),
    }
    if extra_metrics is not None:
        metrics.update(extra_metrics)
    _print_result(
        build_command_result(
            command="simulate",
            inputs=[],
            outputs=outputs,
            metrics=metrics,
            data=report,
        ),
        json_output=args.json,
    )

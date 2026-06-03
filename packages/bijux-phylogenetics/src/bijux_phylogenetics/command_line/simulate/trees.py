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
    simulate_multispecies_coalescent_gene_tree,
    simulate_random_trees,
    write_coalescent_skyline_table,
    write_coalescent_waiting_time_table,
    write_multispecies_coalescent_branch_table,
    write_multispecies_coalescent_event_table,
    write_simulated_tree,
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
    simulate_coalescent.add_argument(
        "--waiting-time-tolerance", type=float, default=0.2
    )
    simulate_coalescent.add_argument("--seed", type=int, default=1)
    simulate_coalescent.add_argument("--out", required=True, type=Path)
    simulate_coalescent.add_argument("--record-table-out", type=Path)
    simulate_coalescent.add_argument("--envelope-table-out", type=Path)
    simulate_coalescent.add_argument("--waiting-time-table-out", type=Path)
    simulate_coalescent.add_argument("--skyline-table-out", type=Path)
    simulate_coalescent.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_coalescent)

    simulate_multispecies_coalescent = simulate_subparsers.add_parser(
        "gene-tree-multispecies-coalescent",
        help="Simulate one gene tree inside a rooted ultrametric species tree.",
    )
    simulate_multispecies_coalescent.add_argument("tree", type=Path)
    simulate_multispecies_coalescent.add_argument(
        "--population-size", type=float, default=1.0
    )
    simulate_multispecies_coalescent.add_argument("--population-size-table", type=Path)
    simulate_multispecies_coalescent.add_argument("--sample-count-table", type=Path)
    simulate_multispecies_coalescent.add_argument("--seed", type=int, default=1)
    simulate_multispecies_coalescent.add_argument("--out", required=True, type=Path)
    simulate_multispecies_coalescent.add_argument("--event-table-out", type=Path)
    simulate_multispecies_coalescent.add_argument("--branch-table-out", type=Path)
    simulate_multispecies_coalescent.add_argument(
        "--json", action="store_true", help="Emit the simulation report as JSON."
    )
    _add_manifest_argument(simulate_multispecies_coalescent)


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
    if args.simulate_command == "tree-coalescent":
        trees, report = simulate_coalescent_trees(
            tree_count=args.tree_count,
            tip_count=args.tip_count,
            population_size=args.population_size,
            waiting_time_tolerance=args.waiting_time_tolerance,
            seed=args.seed,
        )
        _print_tree_result(
            args,
            trees=trees,
            report=report,
            extra_outputs=(
                (
                    []
                    if args.waiting_time_table_out is None
                    else [
                        write_coalescent_waiting_time_table(
                            args.waiting_time_table_out,
                            report,
                        )
                    ]
                )
                + (
                    []
                    if args.skyline_table_out is None
                    else [
                        write_coalescent_skyline_table(
                            args.skyline_table_out,
                            report,
                        )
                    ]
                )
            ),
            extra_metrics={
                "waiting_time_lineage_count": len(report.coalescent_waiting_time_rows),
                "waiting_time_within_tolerance_count": sum(
                    1
                    for row in report.coalescent_waiting_time_rows
                    if row.within_tolerance
                ),
                "waiting_time_all_within_tolerance": all(
                    row.within_tolerance for row in report.coalescent_waiting_time_rows
                ),
                "skyline_interval_count": len(report.coalescent_skyline_rows),
                "skyline_high_uncertainty_count": sum(
                    1
                    for row in report.coalescent_skyline_rows
                    if row.uncertainty_flag == "high"
                ),
            },
        )
        return 0
    if args.simulate_command != "gene-tree-multispecies-coalescent":
        return None

    tree, report = simulate_multispecies_coalescent_gene_tree(
        args.tree,
        default_population_size=args.population_size,
        sample_count_table_path=args.sample_count_table,
        population_size_table_path=args.population_size_table,
        seed=args.seed,
    )
    _print_multispecies_coalescent_result(
        args,
        tree=tree,
        report=report,
    )
    return 0


def _print_tree_result(
    args: Any,
    *,
    trees: list[Any],
    report: Any,
    extra_outputs: list[Path] | None = None,
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
    outputs_to_finalize.extend(extra_outputs or [])
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


def _print_multispecies_coalescent_result(
    args: Any,
    *,
    tree: Any,
    report: Any,
) -> None:
    output_path = write_simulated_tree(args.out, tree)
    outputs_to_finalize = [output_path]
    if args.event_table_out is not None:
        outputs_to_finalize.append(
            write_multispecies_coalescent_event_table(args.event_table_out, report)
        )
    if args.branch_table_out is not None:
        outputs_to_finalize.append(
            write_multispecies_coalescent_branch_table(args.branch_table_out, report)
        )
    inputs = [args.tree]
    if args.sample_count_table is not None:
        inputs.append(args.sample_count_table)
    if args.population_size_table is not None:
        inputs.append(args.population_size_table)
    outputs = _finalize_outputs(
        args,
        command="simulate",
        inputs=inputs,
        outputs=outputs_to_finalize,
    )
    _print_result(
        build_command_result(
            command="simulate",
            inputs=inputs,
            outputs=outputs,
            metrics={
                "species_tip_count": report.species_tip_count,
                "gene_tip_count": report.gene_tip_count,
                "coalescent_event_count": len(report.event_rows),
                "species_branch_count": len(report.branch_rows),
                "deep_coalescence_total": report.deep_coalescence_total,
            },
            data=report,
        ),
        json_output=args.json,
    )

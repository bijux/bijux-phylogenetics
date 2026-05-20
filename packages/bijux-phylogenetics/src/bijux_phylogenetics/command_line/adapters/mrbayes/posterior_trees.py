from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    parse_mrbayes_posterior_tree_samples,
    subsample_mrbayes_posterior_tree_set,
    summarize_mrbayes_posterior_trees,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_mrbayes_posterior_tree_commands(adapter_subparsers: Any) -> None:
    adapter_mrbayes_summarize = adapter_subparsers.add_parser(
        "mrbayes-summarize",
        help="Summarize MrBayes posterior trees after burn-in removal.",
    )
    adapter_mrbayes_summarize.add_argument("input_path", type=Path)
    adapter_mrbayes_summarize.add_argument(
        "--burnin-fraction", type=float, default=0.25
    )
    adapter_mrbayes_summarize.add_argument(
        "--json", action="store_true", help="Emit the posterior summary as JSON."
    )
    _add_manifest_argument(adapter_mrbayes_summarize)

    adapter_mrbayes_trees = adapter_subparsers.add_parser(
        "mrbayes-trees",
        help="Parse a MrBayes posterior tree set into sampled trees.",
    )
    adapter_mrbayes_trees.add_argument("input_path", type=Path)
    adapter_mrbayes_trees.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree set report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_trees)

    adapter_mrbayes_subsample = adapter_subparsers.add_parser(
        "mrbayes-subsample",
        help="Subsample MrBayes posterior trees while preserving generation metadata.",
    )
    adapter_mrbayes_subsample.add_argument("input_path", type=Path)
    adapter_mrbayes_subsample.add_argument(
        "--method",
        required=True,
        choices=("evenly-spaced", "random"),
        help="Select evenly spaced thinning or a seeded random retained subset.",
    )
    adapter_mrbayes_subsample.add_argument("--burnin-fraction", type=float, default=0.0)
    adapter_mrbayes_subsample.add_argument("--thinning-interval", type=int)
    adapter_mrbayes_subsample.add_argument("--sample-count", type=int)
    adapter_mrbayes_subsample.add_argument("--seed", type=int)
    adapter_mrbayes_subsample.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_mrbayes_subsample.add_argument(
        "--sample-table-out",
        type=Path,
        help="Write a TSV ledger of retained posterior-tree metadata.",
    )
    adapter_mrbayes_subsample.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior subsampling report as JSON.",
    )
    _add_manifest_argument(adapter_mrbayes_subsample)


def run_mrbayes_posterior_tree_command(args: Any) -> int | None:
    if args.adapter_command == "mrbayes-summarize":
        consensus_tree, report = summarize_mrbayes_posterior_trees(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=[report.filtered_tree_set_path],
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "kept_tree_count": report.kept_tree_count,
                    "rooted_topology_count": report.rooted_topology_count,
                    "tip_count": consensus_tree.tip_count,
                },
                data={
                    "summary": report,
                    "consensus_newick": report.consensus_newick,
                },
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-trees":
        report = parse_mrbayes_posterior_tree_samples(args.input_path)
        outputs = _finalize_outputs(args, command="adapter", inputs=[args.input_path])
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "rooted_tree_count": report.rooted_tree_count,
                    "sampled_generation_count": len(report.sampled_generations),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command == "mrbayes-subsample":
        report = subsample_mrbayes_posterior_tree_set(
            args.input_path,
            method=args.method,
            thinning_interval=args.thinning_interval,
            sample_count=args.sample_count,
            burnin_fraction=args.burnin_fraction,
            random_seed=args.seed,
        )
        outputs: list[Path | str] = []
        if args.tree_set_out is not None:
            outputs.append(write_posterior_tree_subsample(args.tree_set_out, report))
        if args.sample_table_out is not None:
            outputs.append(
                write_posterior_tree_subsample_table(
                    args.sample_table_out,
                    report,
                )
            )
        outputs = _finalize_outputs(
            args,
            command="adapter",
            inputs=[args.input_path],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="adapter",
                inputs=[args.input_path],
                outputs=outputs,
                metrics={
                    "total_tree_count": report.total_tree_count,
                    "burnin_tree_count": report.burnin_tree_count,
                    "pre_subsampling_tree_count": report.pre_subsampling_tree_count,
                    "retained_tree_count": report.retained_tree_count,
                    "selection_method": report.selection_method,
                    "retained_generation_count": len(
                        [tree for tree in report.trees if tree.generation is not None]
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None

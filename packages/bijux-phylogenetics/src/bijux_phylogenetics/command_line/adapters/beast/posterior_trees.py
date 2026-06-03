from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.bayesian import (
    parse_beast_posterior_tree_samples,
    subsample_beast_posterior_tree_set,
    write_beast_posterior_tree_set,
    write_posterior_tree_subsample,
    write_posterior_tree_subsample_table,
)
from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result


def add_beast_posterior_tree_commands(adapter_subparsers: Any) -> None:
    adapter_beast_trees = adapter_subparsers.add_parser(
        "beast-trees",
        help="Parse a BEAST posterior tree set into state-tagged normalized trees.",
    )
    adapter_beast_trees.add_argument("input_path", type=Path)
    adapter_beast_trees.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early sampled trees before reporting summaries.",
    )
    adapter_beast_trees.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_trees.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior tree set report as JSON.",
    )
    _add_manifest_argument(adapter_beast_trees)

    adapter_beast_subsample = adapter_subparsers.add_parser(
        "beast-subsample",
        help="Subsample BEAST posterior trees while preserving state metadata.",
    )
    adapter_beast_subsample.add_argument("input_path", type=Path)
    adapter_beast_subsample.add_argument(
        "--method",
        required=True,
        choices=("evenly-spaced", "random"),
        help="Select evenly spaced thinning or a seeded random retained subset.",
    )
    adapter_beast_subsample.add_argument(
        "--burnin-fraction",
        type=float,
        default=0.0,
        help="Discard this fraction of early sampled trees before subsampling.",
    )
    adapter_beast_subsample.add_argument("--thinning-interval", type=int)
    adapter_beast_subsample.add_argument("--sample-count", type=int)
    adapter_beast_subsample.add_argument("--seed", type=int)
    adapter_beast_subsample.add_argument(
        "--tree-set-out",
        type=Path,
        help="Write the retained posterior tree set as normalized Newick.",
    )
    adapter_beast_subsample.add_argument(
        "--sample-table-out",
        type=Path,
        help="Write a TSV ledger of retained posterior-tree metadata.",
    )
    adapter_beast_subsample.add_argument(
        "--json",
        action="store_true",
        help="Emit the posterior subsampling report as JSON.",
    )
    _add_manifest_argument(adapter_beast_subsample)


def run_beast_posterior_tree_command(args: Any) -> int | None:
    if args.adapter_command == "beast-trees":
        report = parse_beast_posterior_tree_samples(
            args.input_path,
            burnin_fraction=args.burnin_fraction,
        )
        outputs: list[Path | str] = []
        if args.tree_set_out is not None:
            outputs.append(write_beast_posterior_tree_set(args.tree_set_out, report))
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
                    "kept_tree_count": report.kept_tree_count,
                    "rooted_tree_count": report.rooted_tree_count,
                    "burnin_fraction": report.burnin_fraction,
                    "clade_count": len(report.clades),
                    "sampled_state_count": len(report.sampled_states),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.adapter_command != "beast-subsample":
        return None

    report = subsample_beast_posterior_tree_set(
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
                "retained_state_count": len(
                    [tree for tree in report.trees if tree.state is not None]
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

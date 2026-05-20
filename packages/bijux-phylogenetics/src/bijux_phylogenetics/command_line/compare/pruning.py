from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    prune_trees_to_shared_taxa,
    write_shared_taxa_pruning_table,
    write_shared_taxa_removed_taxa_table,
    write_tree_comparison_table,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_pruning_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "prune":
        return None

    if args.third is None:
        parser.exit(status=2, message="compare prune requires two tree paths\n")
    if args.out is None:
        parser.exit(
            status=2,
            message="compare prune requires --out as an output directory\n",
        )

    left_path = Path(args.right)
    right_path = Path(args.third)
    pruned_left, pruned_right, report = prune_trees_to_shared_taxa(
        left_path, right_path
    )
    args.out.mkdir(parents=True, exist_ok=True)
    left_out = write_newick(args.out / "left-shared.nwk", pruned_left)
    right_out = write_newick(args.out / "right-shared.nwk", pruned_right)
    pruning_out = write_shared_taxa_pruning_table(
        args.out / "shared-taxa-pruning.tsv",
        left_path,
        right_path,
    )
    removed_out = write_shared_taxa_removed_taxa_table(
        args.out / "shared-taxa-removed.tsv",
        left_path,
        right_path,
    )
    comparison_out = write_tree_comparison_table(
        args.out / "shared-taxa-comparison.tsv",
        left_out,
        right_out,
    )
    outputs = _finalize_outputs(
        args,
        command="compare",
        inputs=[left_path, right_path],
        outputs=[left_out, right_out, pruning_out, removed_out, comparison_out],
    )
    _print_result(
        build_command_result(
            command="compare",
            inputs=[left_path, right_path],
            outputs=outputs,
            metrics={
                "shared_taxa": len(report.shared_taxa),
                "left_removed_taxa": len(report.left_pruning.removed_taxa),
                "right_removed_taxa": len(report.right_pruning.removed_taxa),
                "topology_equal_after_pruning": (
                    report.post_pruning_comparison.topology_equal
                ),
                "post_pruning_robinson_foulds_distance": (
                    report.post_pruning_comparison.robinson_foulds_distance
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

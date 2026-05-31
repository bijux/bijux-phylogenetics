from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    approximate_maximum_agreement_subtree,
    prune_trees_to_agreement_subtree,
    prune_trees_to_shared_taxa,
    write_agreement_subtree_pruning_table,
    write_agreement_subtree_removed_taxa_table,
    write_agreement_subtree_search_table,
    write_maximum_agreement_subtree_pruning_table,
    write_maximum_agreement_subtree_removed_taxa_table,
    write_maximum_agreement_subtree_search_table,
    write_shared_taxa_pruning_table,
    write_shared_taxa_removed_taxa_table,
    write_tree_comparison_table,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_pruning_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left == "maximum-agreement-subtree":
        return _run_compare_maximum_agreement_subtree_command(args, parser=parser)
    if args.left == "agreement-subtree":
        return _run_compare_agreement_subtree_command(args, parser=parser)
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


def _run_compare_agreement_subtree_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int:
    if args.third is None:
        parser.exit(
            status=2,
            message="compare agreement-subtree requires two tree paths\n",
        )
    if args.out is None:
        parser.exit(
            status=2,
            message="compare agreement-subtree requires --out as an output directory\n",
        )

    left_path = Path(args.right)
    right_path = Path(args.third)
    pruned_left, pruned_right, report = prune_trees_to_agreement_subtree(
        left_path,
        right_path,
        rf_mode=args.rf_mode,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    left_out = write_newick(args.out / "left-agreement.nwk", pruned_left)
    right_out = write_newick(args.out / "right-agreement.nwk", pruned_right)
    pruning_out = write_agreement_subtree_pruning_table(
        args.out / "agreement-subtree-pruning.tsv",
        left_path,
        right_path,
        rf_mode=args.rf_mode,
    )
    removed_out = write_agreement_subtree_removed_taxa_table(
        args.out / "agreement-subtree-removed.tsv",
        left_path,
        right_path,
        rf_mode=args.rf_mode,
    )
    search_out = write_agreement_subtree_search_table(
        args.out / "agreement-subtree-search.tsv",
        left_path,
        right_path,
        rf_mode=args.rf_mode,
    )
    comparison_out = write_tree_comparison_table(
        args.out / "agreement-subtree-comparison.tsv",
        left_out,
        right_out,
    )
    outputs = _finalize_outputs(
        args,
        command="compare",
        inputs=[left_path, right_path],
        outputs=[
            left_out,
            right_out,
            pruning_out,
            removed_out,
            search_out,
            comparison_out,
        ],
    )
    _print_result(
        build_command_result(
            command="compare",
            inputs=[left_path, right_path],
            outputs=outputs,
            metrics={
                "shared_taxa": len(report.shared_taxa),
                "retained_taxa": len(report.retained_taxa),
                "agreement_removed_taxa": len(report.agreement_removed_taxa),
                "evaluated_candidate_count": report.evaluated_candidate_count,
                "rf_mode": report.rf_mode,
                "topology_equal_after_pruning": (
                    report.post_pruning_comparison.robinson_foulds_distance == 0
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


def _run_compare_maximum_agreement_subtree_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int:
    if args.third is None:
        parser.exit(
            status=2,
            message="compare maximum-agreement-subtree requires two tree paths\n",
        )
    if args.out is None:
        parser.exit(
            status=2,
            message=(
                "compare maximum-agreement-subtree requires --out as an output directory\n"
            ),
        )
    if args.max_evaluated_candidates is None:
        parser.exit(
            status=2,
            message=(
                "compare maximum-agreement-subtree requires --max-evaluated-candidates\n"
            ),
        )

    left_path = Path(args.right)
    right_path = Path(args.third)
    pruned_left, pruned_right, report = approximate_maximum_agreement_subtree(
        left_path,
        right_path,
        rf_mode=args.rf_mode,
        max_evaluated_candidate_count=args.max_evaluated_candidates,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    left_out = write_newick(args.out / "left-maximum-agreement.nwk", pruned_left)
    right_out = write_newick(args.out / "right-maximum-agreement.nwk", pruned_right)
    pruning_out = write_maximum_agreement_subtree_pruning_table(
        args.out / "maximum-agreement-subtree-pruning.tsv",
        left_path,
        right_path,
        rf_mode=args.rf_mode,
        max_evaluated_candidate_count=args.max_evaluated_candidates,
    )
    removed_out = write_maximum_agreement_subtree_removed_taxa_table(
        args.out / "maximum-agreement-subtree-removed.tsv",
        left_path,
        right_path,
        rf_mode=args.rf_mode,
        max_evaluated_candidate_count=args.max_evaluated_candidates,
    )
    search_out = write_maximum_agreement_subtree_search_table(
        args.out / "maximum-agreement-subtree-search.tsv",
        left_path,
        right_path,
        rf_mode=args.rf_mode,
        max_evaluated_candidate_count=args.max_evaluated_candidates,
    )
    comparison_out = write_tree_comparison_table(
        args.out / "maximum-agreement-subtree-comparison.tsv",
        left_out,
        right_out,
    )
    outputs = _finalize_outputs(
        args,
        command="compare",
        inputs=[left_path, right_path],
        outputs=[
            left_out,
            right_out,
            pruning_out,
            removed_out,
            search_out,
            comparison_out,
        ],
    )
    _print_result(
        build_command_result(
            command="compare",
            inputs=[left_path, right_path],
            outputs=outputs,
            metrics={
                "shared_taxa": len(report.shared_taxa),
                "retained_taxa": len(report.retained_taxa),
                "approximation_removed_taxa": len(report.approximation_removed_taxa),
                "evaluated_candidate_count": report.evaluated_candidate_count,
                "max_evaluated_candidate_count": (report.max_evaluated_candidate_count),
                "rf_mode": report.rf_mode,
                "approximation_status": report.approximation_status,
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

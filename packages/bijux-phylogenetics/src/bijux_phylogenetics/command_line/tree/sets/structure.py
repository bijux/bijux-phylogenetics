from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    analyze_tree_set_branch_lengths,
    cluster_trees_by_topology,
    compute_tree_distance_matrix,
    summarize_posterior_topology_diversity,
    summarize_tree_set_shapes,
    write_branch_length_table,
    write_tree_distance_distribution_table,
    write_tree_distance_matrix,
    write_tree_shape_table,
)


def add_tree_set_structure_commands(tree_set_subparsers: Any) -> None:
    tree_set_shape = tree_set_subparsers.add_parser(
        "shape",
        help="Summarize tree balance, ladderization, and height across a tree set.",
    )
    tree_set_shape.add_argument("tree_set", type=Path)
    tree_set_shape.add_argument(
        "--out", type=Path, help="Write one shape-summary row per tree as TSV."
    )
    tree_set_shape.add_argument(
        "--json", action="store_true", help="Emit the tree-shape report as JSON."
    )
    _add_manifest_argument(tree_set_shape)

    tree_set_branch_lengths = tree_set_subparsers.add_parser(
        "branch-lengths",
        help="Summarize branch-length distributions across a tree set.",
    )
    tree_set_branch_lengths.add_argument("tree_set", type=Path)
    tree_set_branch_lengths.add_argument(
        "--out", type=Path, help="Write one row per branch as TSV."
    )
    tree_set_branch_lengths.add_argument(
        "--json",
        action="store_true",
        help="Emit the branch-length distribution report as JSON.",
    )
    _add_manifest_argument(tree_set_branch_lengths)

    tree_set_distances = tree_set_subparsers.add_parser(
        "distance-matrix",
        help="Compute pairwise RF distances across a tree set.",
    )
    tree_set_distances.add_argument("tree_set", type=Path)
    tree_set_distances.add_argument(
        "--out", type=Path, help="Write the pairwise distance table as TSV."
    )
    tree_set_distances.add_argument(
        "--json", action="store_true", help="Emit the distance report as JSON."
    )
    _add_manifest_argument(tree_set_distances)

    tree_set_diversity_summary = tree_set_subparsers.add_parser(
        "diversity",
        help="Summarize topology diversity and RF distribution across one tree set.",
    )
    tree_set_diversity_summary.add_argument("tree_set", type=Path)
    tree_set_diversity_summary.add_argument(
        "--out",
        type=Path,
        help="Write the RF-distribution table as TSV.",
    )
    tree_set_diversity_summary.add_argument(
        "--json", action="store_true", help="Emit the diversity report as JSON."
    )
    _add_manifest_argument(tree_set_diversity_summary)

    tree_set_clusters = tree_set_subparsers.add_parser(
        "cluster",
        help="Cluster trees by identical rooted topology signatures.",
    )
    tree_set_clusters.add_argument("tree_set", type=Path)
    tree_set_clusters.add_argument(
        "--json", action="store_true", help="Emit the cluster report as JSON."
    )
    _add_manifest_argument(tree_set_clusters)


def run_tree_set_structure_command(args: Any) -> int | None:
    if args.tree_set_command == "shape":
        report = summarize_tree_set_shapes(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_tree_shape_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "balanced_tree_count": report.aggregate.balanced_tree_count,
                    "ladderized_tree_count": report.aggregate.ladderized_tree_count,
                    "star_like_tree_count": report.aggregate.star_like_tree_count,
                    "comb_like_tree_count": report.aggregate.comb_like_tree_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "branch-lengths":
        report = analyze_tree_set_branch_lengths(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_branch_length_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "branch_count": report.aggregate.branch_count,
                    "zero_length_branch_count": report.aggregate.zero_length_branch_count,
                    "negative_branch_count": report.aggregate.negative_branch_count,
                    "long_outlier_count": report.aggregate.long_outlier_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "distance-matrix":
        report = compute_tree_distance_matrix(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_tree_distance_matrix(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "pair_count": len(report.pairs),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "diversity":
        report = summarize_posterior_topology_diversity(args.tree_set)
        outputs: list[Path] = []
        if args.out is not None:
            outputs.append(write_tree_distance_distribution_table(args.out, report))
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.tree_set],
            outputs=outputs,
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "rooted_topology_count": report.rooted_topology_count,
                    "pair_count": report.pair_count,
                    "rf_bucket_count": len(report.rf_distribution),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "cluster":
        report = cluster_trees_by_topology(args.tree_set)
        outputs = _finalize_outputs(args, command="tree-set", inputs=[args.tree_set])
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.tree_set],
                outputs=outputs,
                metrics={
                    "tree_count": report.tree_count,
                    "runtime_seconds": report.processing.runtime_seconds,
                    "peak_memory_bytes": report.processing.peak_memory_bytes,
                    "skipped_malformed_tree_count": (
                        report.processing.skipped_malformed_tree_count
                    ),
                    "cluster_count": len(report.clusters),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None

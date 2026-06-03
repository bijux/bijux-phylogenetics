from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.runtime.results import build_command_result
from bijux_phylogenetics.trees import (
    compare_posterior_topological_diversity,
    compare_posterior_tree_sets,
)


def add_tree_set_comparison_commands(tree_set_subparsers: Any) -> None:
    tree_set_compare = tree_set_subparsers.add_parser(
        "compare",
        help="Compare two posterior tree sets over clade support and topology distance.",
    )
    tree_set_compare.add_argument("left", type=Path)
    tree_set_compare.add_argument("right", type=Path)
    tree_set_compare.add_argument(
        "--json", action="store_true", help="Emit the comparison report as JSON."
    )
    _add_manifest_argument(tree_set_compare)

    tree_set_diversity = tree_set_subparsers.add_parser(
        "diversity-compare",
        help="Compare posterior topological diversity across two analyses.",
    )
    tree_set_diversity.add_argument("left", type=Path)
    tree_set_diversity.add_argument("right", type=Path)
    tree_set_diversity.add_argument(
        "--json", action="store_true", help="Emit the diversity report as JSON."
    )
    _add_manifest_argument(tree_set_diversity)


def run_tree_set_comparison_command(args: Any) -> int | None:
    if args.tree_set_command == "compare":
        report = compare_posterior_tree_sets(args.left, args.right)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.left, args.right],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.left, args.right],
                outputs=outputs,
                metrics={
                    "left_tree_count": report.left_tree_count,
                    "right_tree_count": report.right_tree_count,
                    "shared_rooted_topology_count": report.shared_rooted_topology_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.tree_set_command == "diversity-compare":
        report = compare_posterior_topological_diversity(args.left, args.right)
        outputs = _finalize_outputs(
            args,
            command="tree-set",
            inputs=[args.left, args.right],
        )
        _print_result(
            build_command_result(
                command="tree-set",
                inputs=[args.left, args.right],
                outputs=outputs,
                metrics={
                    "left_rooted_topology_count": (
                        report.left_summary.rooted_topology_count
                    ),
                    "right_rooted_topology_count": (
                        report.right_summary.rooted_topology_count
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None

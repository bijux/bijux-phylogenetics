from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import (
    _add_distance_tree_method_argument,
    _add_manifest_argument,
    _add_missing_distance_policy_argument,
)
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.distance import (
    build_distance_tree,
    compare_distance_tree_to_reference_tree,
    compare_distance_tree_topologies,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.results import build_command_result

from .shared import (
    add_ambiguity_policy_option,
    add_distance_model_option,
    add_gap_handling_option,
)


def add_distance_tree_commands(alignment_subparsers: Any) -> None:
    alignment_build_tree = alignment_subparsers.add_parser(
        "build-tree",
        help="Build a neighbor-joining, BIONJ, UPGMA, WPGMA, single-linkage, or complete-linkage tree from a DNA distance matrix.",
    )
    alignment_build_tree.add_argument("alignment", type=Path)
    _add_distance_tree_method_argument(alignment_build_tree)
    _add_missing_distance_policy_argument(alignment_build_tree)
    add_distance_model_option(alignment_build_tree)
    add_gap_handling_option(alignment_build_tree)
    add_ambiguity_policy_option(alignment_build_tree)
    alignment_build_tree.add_argument("--out", required=True, type=Path)
    alignment_build_tree.add_argument(
        "--json", action="store_true", help="Emit the build report as JSON."
    )
    _add_manifest_argument(alignment_build_tree)

    alignment_compare_distance_trees = alignment_subparsers.add_parser(
        "compare-distance-trees",
        help="Compare NJ and UPGMA topologies built from the same DNA alignment.",
    )
    alignment_compare_distance_trees.add_argument("alignment", type=Path)
    _add_missing_distance_policy_argument(alignment_compare_distance_trees)
    add_distance_model_option(alignment_compare_distance_trees)
    add_gap_handling_option(alignment_compare_distance_trees)
    add_ambiguity_policy_option(alignment_compare_distance_trees)
    alignment_compare_distance_trees.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(alignment_compare_distance_trees)

    alignment_compare_distance_reference = alignment_subparsers.add_parser(
        "compare-distance-to-tree",
        help="Compare one built distance tree against an external inferred or reviewer reference tree.",
    )
    alignment_compare_distance_reference.add_argument("alignment", type=Path)
    alignment_compare_distance_reference.add_argument("reference_tree", type=Path)
    _add_distance_tree_method_argument(alignment_compare_distance_reference)
    _add_missing_distance_policy_argument(alignment_compare_distance_reference)
    add_distance_model_option(alignment_compare_distance_reference)
    add_gap_handling_option(alignment_compare_distance_reference)
    add_ambiguity_policy_option(alignment_compare_distance_reference)
    alignment_compare_distance_reference.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(alignment_compare_distance_reference)


def run_distance_tree_command(args: Any) -> int | None:
    if args.alignment_command == "build-tree":
        tree, report = build_distance_tree(
            args.alignment,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            missing_distance_policy=args.missing_distance_policy,
        )
        output_path = write_newick(args.out, tree)
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "taxon_count": report.taxon_count,
                    "pair_count": report.pair_count,
                    "method": report.method,
                    "ambiguity_policy": report.ambiguity_policy,
                    "missing_distance_policy": report.missing_distance_policy_report.policy,
                    "imputed_pair_count": len(
                        report.missing_distance_policy_report.imputed_rows
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "compare-distance-trees":
        report = compare_distance_tree_topologies(
            args.alignment,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            missing_distance_policy=args.missing_distance_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment],
                outputs=outputs,
                metrics={
                    "shared_taxa": len(report.shared_taxa),
                    "robinson_foulds_distance": report.robinson_foulds_distance,
                    "same_unrooted_topology": report.same_unrooted_topology,
                    "ambiguity_policy": report.ambiguity_policy,
                    "missing_distance_policy": args.missing_distance_policy,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.alignment_command == "compare-distance-to-tree":
        report = compare_distance_tree_to_reference_tree(
            args.alignment,
            args.reference_tree,
            method=args.method,
            model=args.model,
            gap_handling=args.gap_handling,
            ambiguity_policy=args.ambiguity_policy,
            missing_distance_policy=args.missing_distance_policy,
        )
        outputs = _finalize_outputs(
            args,
            command="alignment",
            inputs=[args.alignment, args.reference_tree],
        )
        _print_result(
            build_command_result(
                command="alignment",
                inputs=[args.alignment, args.reference_tree],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "topology_equal": report.topology.topology_equal,
                    "same_unrooted_topology": report.topology.same_unrooted_topology,
                    "shared_taxa": len(report.topology.shared_taxa),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None

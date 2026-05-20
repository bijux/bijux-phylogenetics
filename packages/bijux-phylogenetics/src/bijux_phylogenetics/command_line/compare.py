from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.reports import build_tree_comparison_report
from bijux_phylogenetics.compare.taxon_influence import (
    analyze_taxon_influence,
    write_taxon_influence_table,
)
from bijux_phylogenetics.compare.topology import (
    compare_branch_lengths,
    compare_clade_overlap,
    compare_support_values,
    compare_topology_distance,
    compare_tree_paths,
    detect_clade_changes,
    prune_trees_to_shared_taxa,
    write_clade_overlap_table,
    write_shared_taxa_pruning_table,
    write_shared_taxa_removed_taxa_table,
    write_support_comparison_table,
    write_topology_distance_split_table,
    write_tree_comparison_table,
)
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.runtime.results import build_command_result


def add_compare_command(subparsers: Any) -> None:
    compare = subparsers.add_parser(
        get_command_spec("compare").name, help=get_command_spec("compare").summary
    )
    compare.add_argument("left")
    compare.add_argument("right")
    compare.add_argument("third", nargs="?")
    compare.add_argument(
        "--tree",
        dest="extra_trees",
        action="append",
        type=Path,
        help="Add another tree path for compare clades.",
    )
    compare.add_argument("--out", type=Path)
    compare.add_argument(
        "--split-table-out",
        type=Path,
        help="Write the clades or splits used by the topology-distance comparison as TSV.",
    )
    compare.add_argument(
        "--rf-mode",
        choices=("rooted", "unrooted"),
        default="rooted",
        help="Compute Robinson-Foulds distance on rooted clades or unrooted bipartitions.",
    )
    compare.add_argument(
        "--taxon-overlap-policy",
        choices=("prune-to-shared", "require-identical"),
        default="prune-to-shared",
        help="Either prune both trees to shared taxa or require identical taxon sets.",
    )
    compare.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(compare)


def run_compare_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    if args.left == "report":
        if args.third is None:
            parser.exit(status=2, message="compare report requires two tree paths\n")
        if args.out is None:
            parser.exit(status=2, message="compare report requires --out\n")
        left_path = Path(args.right)
        right_path = Path(args.third)
        report = build_tree_comparison_report(left_path, right_path, out_path=args.out)
        outputs = _finalize_outputs(
            args,
            command="compare",
            inputs=[left_path, right_path],
            outputs=[args.out],
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=[left_path, right_path],
                outputs=outputs,
                metrics={"shared_taxa": len(report.topology.shared_taxa)},
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.left == "support":
        if args.third is None:
            parser.exit(status=2, message="compare support requires two tree paths\n")
        left_path = Path(args.right)
        right_path = Path(args.third)
        report = compare_support_values(left_path, right_path)
        output_paths: list[Path | str] = []
        if args.out is not None:
            output_paths.append(
                write_support_comparison_table(args.out, left_path, right_path)
            )
        outputs = _finalize_outputs(
            args,
            command="compare",
            inputs=[left_path, right_path],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=[left_path, right_path],
                outputs=outputs,
                metrics={
                    "shared_clades": len(report.shared_clades),
                    "support_disagreements": sum(
                        1
                        for row in report.shared_clades
                        if row.support_disagreement
                    ),
                    "high_support_conflicts": sum(
                        1
                        for row in report.conflicting_clades
                        if row.conflict_classification == "high_support_conflict"
                    ),
                    "low_support_disagreements": sum(
                        1
                        for row in report.conflicting_clades
                        if row.conflict_classification == "low_support_disagreement"
                    ),
                    "moderate_support_disagreements": sum(
                        1
                        for row in report.conflicting_clades
                        if row.conflict_classification
                        == "moderate_support_disagreement"
                    ),
                    "support_unavailable_conflicts": sum(
                        1
                        for row in report.conflicting_clades
                        if row.conflict_classification == "support_unavailable"
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.left == "influence":
        if args.third is None:
            parser.exit(status=2, message="compare influence requires two tree paths\n")
        left_path = Path(args.right)
        right_path = Path(args.third)
        report = analyze_taxon_influence(left_path, right_path)
        output_paths: list[Path | str] = []
        if args.out is not None:
            output_paths.append(
                write_taxon_influence_table(args.out, left_path, right_path)
            )
        outputs = _finalize_outputs(
            args,
            command="compare",
            inputs=[left_path, right_path],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=[left_path, right_path],
                outputs=outputs,
                metrics={
                    "shared_taxa": len(report.shared_taxa),
                    "top_influential_taxon": (
                        report.rows[0].taxon if report.rows else None
                    ),
                    "taxa_with_topology_change": sum(
                        1 for row in report.rows if row.topology_changed
                    ),
                    "taxa_with_support_change": sum(
                        1 for row in report.rows if row.support_changed
                    ),
                    "maximum_influence_score": (
                        report.rows[0].influence_score if report.rows else 0.0
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.left == "clades":
        if args.third is None:
            parser.exit(status=2, message="compare clades requires two tree paths\n")
        tree_paths = [Path(args.right), Path(args.third)]
        if args.extra_trees:
            tree_paths.extend(args.extra_trees)
        report = compare_clade_overlap(tree_paths)
        output_paths: list[Path | str] = []
        if args.out is not None:
            output_paths.append(write_clade_overlap_table(args.out, tree_paths))
        outputs = _finalize_outputs(
            args,
            command="compare",
            inputs=tree_paths,
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=tree_paths,
                outputs=outputs,
                metrics={
                    "shared_clades": len(report.shared_clades),
                    "conflicting_clades": len(report.conflicting_clades),
                    "tree_count": len(report.tree_paths),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.left == "prune":
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
            outputs=[
                left_out,
                right_out,
                pruning_out,
                removed_out,
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

    if args.left == "changes":
        if args.third is None:
            parser.exit(status=2, message="compare changes requires two tree paths\n")
        left_path = Path(args.right)
        right_path = Path(args.third)
        report = detect_clade_changes(left_path, right_path)
        outputs = _finalize_outputs(
            args, command="compare", inputs=[left_path, right_path]
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=[left_path, right_path],
                outputs=outputs,
                metrics={
                    "lost_clades": len(report.lost_clades),
                    "gained_clades": len(report.gained_clades),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.left == "branch-lengths":
        if args.third is None:
            parser.exit(
                status=2,
                message="compare branch-lengths requires two tree paths\n",
            )
        left_path = Path(args.right)
        right_path = Path(args.third)
        report = compare_branch_lengths(
            left_path,
            right_path,
            taxon_overlap_policy=args.taxon_overlap_policy,
        )
        outputs = _finalize_outputs(
            args, command="compare", inputs=[left_path, right_path]
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=[left_path, right_path],
                outputs=outputs,
                metrics={
                    "shared_taxa": len(report.shared_taxa),
                    "same_taxon_set": report.same_taxon_set,
                    "shared_splits": len(report.shared_splits),
                    "branch_score_distance": report.branch_score.branch_score_distance,
                    "missing_length_splits": report.branch_score.missing_length_split_count,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.left == "table":
        if args.third is None:
            parser.exit(status=2, message="compare table requires two tree paths\n")
        if args.out is None:
            parser.exit(status=2, message="compare table requires --out\n")
        left_path = Path(args.right)
        right_path = Path(args.third)
        output_path = write_tree_comparison_table(args.out, left_path, right_path)
        outputs = _finalize_outputs(
            args,
            command="compare",
            inputs=[left_path, right_path],
            outputs=[output_path],
        )
        _print_result(
            build_command_result(
                command="compare",
                inputs=[left_path, right_path],
                outputs=outputs,
                metrics={
                    "table_rows": sum(
                        1
                        for _ in output_path.read_text(
                            encoding="utf-8"
                        ).splitlines()[1:]
                    )
                },
                data={"table_path": output_path},
            ),
            json_output=args.json,
        )
        return 0

    left_path = Path(args.left)
    right_path = Path(args.right)
    report = compare_tree_paths(
        left_path,
        right_path,
        rf_mode=args.rf_mode,
        taxon_overlap_policy=args.taxon_overlap_policy,
    )
    outputs: list[Path | str] = []
    if args.split_table_out is not None:
        outputs.append(
            write_topology_distance_split_table(
                args.split_table_out,
                left_path,
                right_path,
                rf_mode=args.rf_mode,
                taxon_overlap_policy=args.taxon_overlap_policy,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="compare",
        inputs=[left_path, right_path],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="compare",
            inputs=[left_path, right_path],
            outputs=outputs,
            metrics={
                "shared_taxa": len(report.shared_taxa),
                "robinson_foulds_distance": report.robinson_foulds_distance,
                "rf_mode": report.rf_mode,
                "taxon_overlap_policy": report.taxon_overlap_policy,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

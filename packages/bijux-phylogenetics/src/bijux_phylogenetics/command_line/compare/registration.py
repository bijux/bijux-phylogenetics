from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.registry import get_command_spec
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import compare_branch_lengths
from bijux_phylogenetics.runtime.results import build_command_result

from .clades import run_compare_clade_command
from .influence import run_compare_influence_command
from .presentation import run_compare_presentation_command
from .pruning import run_compare_pruning_command
from .support import run_compare_support_command
from .topology_distance import run_compare_topology_distance_command


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
    presentation_result = run_compare_presentation_command(args, parser=parser)
    if presentation_result is not None:
        return presentation_result
    support_result = run_compare_support_command(args, parser=parser)
    if support_result is not None:
        return support_result
    influence_result = run_compare_influence_command(args, parser=parser)
    if influence_result is not None:
        return influence_result
    clade_result = run_compare_clade_command(args, parser=parser)
    if clade_result is not None:
        return clade_result
    pruning_result = run_compare_pruning_command(args, parser=parser)
    if pruning_result is not None:
        return pruning_result
    topology_distance_result = run_compare_topology_distance_command(
        args, parser=parser
    )
    if topology_distance_result is not None:
        return topology_distance_result

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

    parser.error("compare requires a supported workflow or two tree paths")

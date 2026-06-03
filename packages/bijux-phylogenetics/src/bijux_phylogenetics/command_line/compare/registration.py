from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.registry import get_command_spec

from .branch_lengths import run_compare_branch_lengths_command
from .clade_ages import run_compare_clade_ages_command
from .clades import run_compare_clade_command
from .coalescence import run_compare_deep_coalescence_command
from .influence import run_compare_influence_command
from .presentation import run_compare_presentation_command
from .pruning import run_compare_pruning_command
from .reconciliation import run_compare_duplication_loss_transfer_command
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
        "--taxon-map",
        type=Path,
        help="Map gene-tree tip labels onto species-tree taxa for deep-coalescence workflows.",
    )
    compare.add_argument(
        "--mapping-out",
        type=Path,
        help="Write the resolved gene-to-species taxon mapping as TSV.",
    )
    compare.add_argument(
        "--duplication-cost",
        type=float,
        default=2.0,
        help="Per-event cost assigned to duplication steps in DLT reconciliation.",
    )
    compare.add_argument(
        "--loss-cost",
        type=float,
        default=1.0,
        help="Per-branch cost assigned to loss steps in DLT reconciliation.",
    )
    compare.add_argument(
        "--transfer-cost",
        type=float,
        default=3.0,
        help="Per-event cost assigned to transfer or switch steps in DLT reconciliation.",
    )
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
        "--max-evaluated-candidates",
        type=int,
        help=(
            "Explicit candidate budget for heuristic compare workflows that search retained taxon subsets."
        ),
    )
    compare.add_argument("--json", action="store_true", help="Emit the report as JSON.")
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
    clade_age_result = run_compare_clade_ages_command(args, parser=parser)
    if clade_age_result is not None:
        return clade_age_result
    dlt_result = run_compare_duplication_loss_transfer_command(args, parser=parser)
    if dlt_result is not None:
        return dlt_result
    deep_coalescence_result = run_compare_deep_coalescence_command(args, parser=parser)
    if deep_coalescence_result is not None:
        return deep_coalescence_result
    pruning_result = run_compare_pruning_command(args, parser=parser)
    if pruning_result is not None:
        return pruning_result
    branch_length_result = run_compare_branch_lengths_command(args, parser=parser)
    if branch_length_result is not None:
        return branch_length_result
    topology_distance_result = run_compare_topology_distance_command(
        args, parser=parser
    )
    if topology_distance_result is not None:
        return topology_distance_result

    parser.error("compare requires a supported workflow or two tree paths")
    return 2

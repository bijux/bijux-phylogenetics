from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import compare_branch_lengths
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_branch_lengths_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "branch-lengths":
        return None

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
    outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
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

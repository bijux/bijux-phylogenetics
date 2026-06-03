from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    compare_clade_overlap,
    detect_clade_changes,
    write_clade_overlap_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_clade_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left == "clades":
        return _run_compare_clades_command(args, parser=parser)
    if args.left == "changes":
        return _run_compare_changes_command(args, parser=parser)
    return None


def _run_compare_clades_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
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


def _run_compare_changes_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
    if args.third is None:
        parser.exit(status=2, message="compare changes requires two tree paths\n")

    left_path = Path(args.right)
    right_path = Path(args.third)
    report = detect_clade_changes(left_path, right_path)
    outputs = _finalize_outputs(args, command="compare", inputs=[left_path, right_path])
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

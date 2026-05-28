from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    compare_tree_paths,
    write_topology_distance_split_table,
    write_tree_comparison_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_topology_distance_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left == "table":
        return _run_compare_table_command(args, parser=parser)

    if args.left in {
        "clades",
        "prune",
        "agreement-subtree",
        "maximum-agreement-subtree",
        "changes",
        "branch-lengths",
        "clade-ages",
        "deep-coalescence",
        "duplication-loss-transfer",
    }:
        return None

    return _run_compare_topology_distance(args)


def _run_compare_table_command(args: Any, *, parser: argparse.ArgumentParser) -> int:
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
                    1 for _ in output_path.read_text(encoding="utf-8").splitlines()[1:]
                )
            },
            data={"table_path": output_path},
        ),
        json_output=args.json,
    )
    return 0


def _run_compare_topology_distance(args: Any) -> int:
    left_path = Path(args.left)
    right_path = Path(args.right)
    report = compare_tree_paths(
        left_path,
        right_path,
        rf_mode=args.rf_mode,
        taxon_overlap_policy=args.taxon_overlap_policy,
    )
    output_paths: list[Path | str] = []
    if args.split_table_out is not None:
        output_paths.append(
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
        outputs=output_paths,
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

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    compare_clade_ages,
    write_date_aware_tree_comparison_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_clade_ages_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "clade-ages":
        return None
    if args.third is None:
        parser.exit(status=2, message="compare clade-ages requires two tree paths\n")

    left_path = Path(args.right)
    right_path = Path(args.third)
    report = compare_clade_ages(
        left_path,
        right_path,
        taxon_overlap_policy=args.taxon_overlap_policy,
    )
    output_paths: list[Path | str] = []
    if args.out is not None:
        output_paths.append(
            write_date_aware_tree_comparison_table(
                args.out,
                left_path,
                right_path,
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
                "matched_clades": report.matched_clade_count,
                "age_rmse": report.age_rmse,
                "unstable_clades": report.unstable_clade_count,
                "comparison_scope": report.comparison_scope,
                "topology_equal": report.topology.topology_equal,
                "robinson_foulds_distance": (report.topology.robinson_foulds_distance),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

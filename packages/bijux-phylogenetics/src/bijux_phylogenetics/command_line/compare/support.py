from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.topology import (
    compare_support_values,
    write_support_comparison_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_support_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "support":
        return None

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
                    1 for row in report.shared_clades if row.support_disagreement
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
                    if row.conflict_classification == "moderate_support_disagreement"
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

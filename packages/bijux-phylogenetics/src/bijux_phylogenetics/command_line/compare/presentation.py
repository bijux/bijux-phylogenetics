from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.presentation import build_tree_comparison_report
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_presentation_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "report":
        return None

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

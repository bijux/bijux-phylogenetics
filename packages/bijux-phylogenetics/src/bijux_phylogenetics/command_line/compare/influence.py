from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.compare.influence import (
    analyze_taxon_influence,
    write_taxon_influence_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def run_compare_influence_command(
    args: Any, *, parser: argparse.ArgumentParser
) -> int | None:
    if args.left != "influence":
        return None

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
                "top_influential_taxon": report.rows[0].taxon if report.rows else None,
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

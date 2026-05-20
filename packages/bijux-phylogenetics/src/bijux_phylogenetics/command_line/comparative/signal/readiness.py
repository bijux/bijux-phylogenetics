from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.common import (
    summarize_numeric_trait,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_signal_readiness_commands(comparative_subparsers: Any) -> None:
    comparative_readiness = comparative_subparsers.add_parser(
        "readiness",
        help="Check whether a rooted tree and numeric trait are ready for comparative analysis.",
    )
    comparative_readiness.add_argument("tree", type=Path)
    comparative_readiness.add_argument("table", type=Path)
    comparative_readiness.add_argument("--trait", required=True)
    comparative_readiness.add_argument("--taxon-column")
    comparative_readiness.add_argument(
        "--json", action="store_true", help="Emit the report as JSON."
    )
    _add_manifest_argument(comparative_readiness)

    comparative_summarize = comparative_subparsers.add_parser(
        "summarize",
        help="Summarize a numeric trait after pruning to overlapping phylogenetic taxa.",
    )
    comparative_summarize.add_argument("tree", type=Path)
    comparative_summarize.add_argument("table", type=Path)
    comparative_summarize.add_argument("--trait", required=True)
    comparative_summarize.add_argument("--taxon-column")
    comparative_summarize.add_argument(
        "--json", action="store_true", help="Emit the summary as JSON."
    )
    _add_manifest_argument(comparative_summarize)


def run_signal_readiness_command(args: Any) -> int | None:
    if args.comparative_command == "readiness":
        report = summarize_numeric_trait_readiness(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                warnings=report.warnings,
                metrics={
                    "tree_taxa": report.tree_taxa,
                    "analysis_taxa": len(report.analysis_taxa),
                    "ready": report.ready,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    if args.comparative_command == "summarize":
        report = summarize_numeric_trait(
            args.tree,
            args.table,
            trait=args.trait,
            taxon_column=args.taxon_column,
        )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "taxon_count": report.taxon_count,
                    "mean": report.mean,
                    "variance": report.variance,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0

    return None

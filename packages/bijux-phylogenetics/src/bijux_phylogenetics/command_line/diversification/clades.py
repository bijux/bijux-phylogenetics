from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative import (
    detect_diversification_outlier_clades,
    write_clade_diversification_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_diversification_clade_command(diversification_subparsers: Any) -> None:
    diversification_clades = diversification_subparsers.add_parser(
        "clades",
        help="Detect clades with unusually high or low diversification.",
    )
    diversification_clades.add_argument("tree", type=Path)
    diversification_clades.add_argument(
        "--model", choices=("yule", "birth-death"), default="birth-death"
    )
    diversification_clades.add_argument("--min-tip-count", type=int, default=2)
    diversification_clades.add_argument(
        "--out", type=Path, help="Write the clade diversification table as TSV."
    )
    diversification_clades.add_argument(
        "--json", action="store_true", help="Emit the clade scan report as JSON."
    )
    _add_manifest_argument(diversification_clades)


def run_diversification_clade_command(args: Any) -> int | None:
    if args.diversification_command != "clades":
        return None

    report = detect_diversification_outlier_clades(
        args.tree,
        min_tip_count=args.min_tip_count,
        model=args.model,
    )
    outputs: list[Path | str] = []
    if args.out is not None:
        outputs.append(write_clade_diversification_table(args.out, report))
    outputs = _finalize_outputs(
        args,
        command="diversification",
        inputs=[args.tree],
        outputs=outputs,
    )
    _print_result(
        build_command_result(
            command="diversification",
            inputs=[args.tree],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "global_rate": report.global_rate,
                "high_clade_count": len(report.high_diversification_clades),
                "low_clade_count": len(report.low_diversification_clades),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

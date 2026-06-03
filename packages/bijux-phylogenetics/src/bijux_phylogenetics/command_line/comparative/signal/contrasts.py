from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.signal import (
    compute_phylogenetic_independent_contrasts,
    summarize_independent_contrast_regression,
    write_independent_contrast_regression_table,
    write_independent_contrast_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_signal_contrast_commands(comparative_subparsers: Any) -> None:
    comparative_contrasts = comparative_subparsers.add_parser(
        "contrasts",
        help="Compute phylogenetic independent contrasts for one numeric trait.",
    )
    comparative_contrasts.add_argument("tree", type=Path)
    comparative_contrasts.add_argument("table", type=Path)
    comparative_contrasts.add_argument("--trait", required=True)
    comparative_contrasts.add_argument(
        "--predictor-trait",
        help="Optional second numeric trait for regression through the origin on matched contrasts.",
    )
    comparative_contrasts.add_argument("--taxon-column")
    comparative_contrasts.add_argument(
        "--contrasts-out",
        type=Path,
        help="Write one flat contrast ledger as TSV or CSV.",
    )
    comparative_contrasts.add_argument(
        "--regression-out",
        type=Path,
        help="Write one regression-through-origin contrast ledger as TSV or CSV.",
    )
    comparative_contrasts.add_argument(
        "--json", action="store_true", help="Emit the contrast report as JSON."
    )
    _add_manifest_argument(comparative_contrasts)


def run_signal_contrast_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    if args.comparative_command != "contrasts":
        return None

    if args.regression_out and not args.predictor_trait:
        parser.error(
            "--regression-out requires --predictor-trait for regression-through-origin output"
        )

    report = compute_phylogenetic_independent_contrasts(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
    )
    regression_report = None
    if args.contrasts_out:
        write_independent_contrast_table(args.contrasts_out, report)
    if args.predictor_trait:
        regression_report = summarize_independent_contrast_regression(
            args.tree,
            args.table,
            response_trait=args.trait,
            predictor_trait=args.predictor_trait,
            taxon_column=args.taxon_column,
        )
        if args.regression_out:
            write_independent_contrast_regression_table(
                args.regression_out,
                regression_report,
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
                "contrast_count": len(report.contrasts),
                "regression_row_count": (
                    len(regression_report.rows) if regression_report is not None else 0
                ),
                "regression_slope": (
                    regression_report.slope if regression_report is not None else None
                ),
                "regression_p_value": (
                    regression_report.p_value if regression_report is not None else None
                ),
            },
            data={
                "contrast_report": report,
                "regression": regression_report,
            },
        ),
        json_output=args.json,
    )
    return 0

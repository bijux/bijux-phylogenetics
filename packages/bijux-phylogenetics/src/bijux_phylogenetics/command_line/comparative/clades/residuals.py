from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.clades.residuals import (
    analyze_comparative_residual_clades,
    write_comparative_residual_clade_table,
    write_comparative_residual_taxon_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_clade_residual_commands(comparative_subparsers: Any) -> None:
    comparative_clade_residuals = comparative_subparsers.add_parser(
        "clade-residuals",
        help="Aggregate comparative model residuals across internal clades.",
    )
    comparative_clade_residuals.add_argument("tree", type=Path)
    comparative_clade_residuals.add_argument("table", type=Path)
    comparative_clade_residuals.add_argument("--response")
    comparative_clade_residuals.add_argument("--predictors", nargs="+")
    comparative_clade_residuals.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass + habitat'.",
    )
    comparative_clade_residuals.add_argument("--taxon-column")
    comparative_clade_residuals.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response residual aggregation requires a numeric value.",
    )
    comparative_clade_residuals.add_argument(
        "--taxa-out",
        type=Path,
        help="Write the analyzed taxon residual ledger as TSV or CSV.",
    )
    comparative_clade_residuals.add_argument(
        "--clades-out",
        type=Path,
        help="Write the internal clade residual aggregation ledger as TSV or CSV.",
    )
    comparative_clade_residuals.add_argument(
        "--json", action="store_true", help="Emit the clade residual report as JSON."
    )
    _add_manifest_argument(comparative_clade_residuals)


def run_clade_residual_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "clade-residuals":
        return None

    report = analyze_comparative_residual_clades(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs: list[Path | str] = []
    if args.taxa_out is not None:
        outputs.append(
            write_comparative_residual_taxon_table(
                args.taxa_out,
                report,
            )
        )
    if args.clades_out is not None:
        outputs.append(
            write_comparative_residual_clade_table(
                args.clades_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    top_clade = (
        min(report.clade_rows, key=lambda row: row.rank).clade_id
        if report.clade_rows
        else None
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "model_family": report.model_family,
                "taxon_count": len(report.taxon_rows),
                "clade_count": len(report.clade_rows),
                "residual_heavy_clade_count": len(report.residual_heavy_clades),
                "top_influential_clade": top_clade,
                "standardized_residual_method": report.standardized_residual_method,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

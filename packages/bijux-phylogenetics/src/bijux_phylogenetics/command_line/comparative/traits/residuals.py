from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_residuals,
    write_phylogenetic_residual_coefficient_table,
    write_phylogenetic_residual_exclusion_table,
    write_phylogenetic_residual_summary_table,
    write_phylogenetic_residual_taxon_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_phylogenetic_residual_review_command(comparative_subparsers: Any) -> None:
    comparative_phylogenetic_residuals = comparative_subparsers.add_parser(
        "phylogenetic-residuals",
        help="Summarize tree-aware fitted values and residuals for one continuous response and predictor.",
    )
    comparative_phylogenetic_residuals.add_argument("tree", type=Path)
    comparative_phylogenetic_residuals.add_argument("table", type=Path)
    comparative_phylogenetic_residuals.add_argument("--response", required=True)
    comparative_phylogenetic_residuals.add_argument("--predictor", required=True)
    comparative_phylogenetic_residuals.add_argument("--taxon-column")
    comparative_phylogenetic_residuals.add_argument(
        "--method",
        choices=("brownian", "lambda"),
        default="lambda",
        help="Use fixed Brownian covariance or estimate Pagel lambda before computing residuals.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--summary-out",
        type=Path,
        help="Write one phylogenetic-residual summary ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--residuals-out",
        type=Path,
        help="Write one taxon-level phylogenetic-residual ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write one phylogenetic-residual coefficient ledger as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for phylogenetic residual review as TSV or CSV.",
    )
    comparative_phylogenetic_residuals.add_argument(
        "--json",
        action="store_true",
        help="Emit the phylogenetic-residual review as JSON.",
    )
    _add_manifest_argument(comparative_phylogenetic_residuals)


def run_phylogenetic_residual_review_command(args: Any) -> int | None:
    if args.comparative_command != "phylogenetic-residuals":
        return None

    report = summarize_phylogenetic_residuals(
        args.tree,
        args.table,
        response=args.response,
        predictor=args.predictor,
        taxon_column=args.taxon_column,
        method=args.method,
    )
    if args.summary_out:
        write_phylogenetic_residual_summary_table(args.summary_out, report)
    if args.residuals_out:
        write_phylogenetic_residual_taxon_table(args.residuals_out, report)
    if args.coefficients_out:
        write_phylogenetic_residual_coefficient_table(
            args.coefficients_out,
            report,
        )
    if args.excluded_taxa_out:
        write_phylogenetic_residual_exclusion_table(
            args.excluded_taxa_out,
            report,
        )
    outputs = _finalize_outputs(
        args, command="comparative", inputs=[args.tree, args.table]
    )
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "method": report.method,
                "outlier_count": len(report.outlier_taxa),
                "top_outlier_taxon": (
                    None
                    if not report.taxon_rows
                    else max(
                        report.taxon_rows,
                        key=lambda row: row.abs_standardized_residual,
                    ).taxon
                ),
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

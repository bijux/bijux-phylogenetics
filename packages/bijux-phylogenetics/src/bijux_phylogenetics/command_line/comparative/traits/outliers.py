from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.traits.outliers import (
    summarize_trait_outliers,
    write_trait_outlier_exclusion_table,
    write_trait_outlier_summary_table,
    write_trait_outlier_taxon_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_trait_outlier_review_command(comparative_subparsers: Any) -> None:
    comparative_trait_outliers = comparative_subparsers.add_parser(
        "trait-outliers",
        help="Rank continuous-trait taxa by conditional phylogenetic residual size.",
    )
    comparative_trait_outliers.add_argument("tree", type=Path)
    comparative_trait_outliers.add_argument("table", type=Path)
    comparative_trait_outliers.add_argument("--trait", required=True)
    comparative_trait_outliers.add_argument("--taxon-column")
    comparative_trait_outliers.add_argument(
        "--summary-out",
        type=Path,
        help="Write one trait-outlier summary ledger as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--outliers-out",
        type=Path,
        help="Write one ranked taxon outlier ledger as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for trait outlier review as TSV or CSV.",
    )
    comparative_trait_outliers.add_argument(
        "--json",
        action="store_true",
        help="Emit the trait outlier review as JSON.",
    )
    _add_manifest_argument(comparative_trait_outliers)


def run_trait_outlier_review_command(args: Any) -> int | None:
    if args.comparative_command != "trait-outliers":
        return None

    report = summarize_trait_outliers(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
    )
    if args.summary_out:
        write_trait_outlier_summary_table(args.summary_out, report)
    if args.outliers_out:
        write_trait_outlier_taxon_table(args.outliers_out, report)
    if args.excluded_taxa_out:
        write_trait_outlier_exclusion_table(
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
                "selected_model": report.selected_model,
                "outlier_count": len(report.outlier_taxa),
                "top_outlier_taxon": report.top_outlier_taxon,
                "top_abs_standardized_residual": report.top_abs_standardized_residual,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

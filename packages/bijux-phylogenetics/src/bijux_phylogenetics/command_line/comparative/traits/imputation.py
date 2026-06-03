from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.traits.imputation import (
    summarize_trait_imputation,
    write_trait_imputation_exclusion_table,
    write_trait_imputation_holdout_table,
    write_trait_imputation_summary_table,
    write_trait_imputation_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_trait_imputation_review_command(comparative_subparsers: Any) -> None:
    comparative_trait_imputation = comparative_subparsers.add_parser(
        "trait-imputation",
        help="Impute missing continuous-trait values under a Brownian phylogenetic model.",
    )
    comparative_trait_imputation.add_argument("tree", type=Path)
    comparative_trait_imputation.add_argument("table", type=Path)
    comparative_trait_imputation.add_argument("--trait", required=True)
    comparative_trait_imputation.add_argument("--taxon-column")
    comparative_trait_imputation.add_argument(
        "--summary-out",
        type=Path,
        help="Write one Brownian trait-imputation summary ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--imputations-out",
        type=Path,
        help="Write one imputed-value ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--holdout-out",
        type=Path,
        help="Write one leave-one-observed-out validation ledger as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write one excluded-taxa ledger for trait imputation as TSV or CSV.",
    )
    comparative_trait_imputation.add_argument(
        "--json",
        action="store_true",
        help="Emit the trait-imputation review as JSON.",
    )
    _add_manifest_argument(comparative_trait_imputation)


def run_trait_imputation_review_command(args: Any) -> int | None:
    if args.comparative_command != "trait-imputation":
        return None

    report = summarize_trait_imputation(
        args.tree,
        args.table,
        trait=args.trait,
        taxon_column=args.taxon_column,
    )
    if args.summary_out:
        write_trait_imputation_summary_table(args.summary_out, report)
    if args.imputations_out:
        write_trait_imputation_table(args.imputations_out, report)
    if args.holdout_out:
        write_trait_imputation_holdout_table(args.holdout_out, report)
    if args.excluded_taxa_out:
        write_trait_imputation_exclusion_table(
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
                "observed_taxon_count": report.observed_taxon_count,
                "imputed_taxon_count": len(report.imputation_rows),
                "excluded_taxon_count": len(report.excluded_taxa),
                "holdout_validation_status": report.holdout_validation_status,
                "holdout_count": len(report.holdout_rows),
                "holdout_mean_absolute_error": report.holdout_mean_absolute_error,
                "holdout_interval_coverage": report.holdout_interval_coverage,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

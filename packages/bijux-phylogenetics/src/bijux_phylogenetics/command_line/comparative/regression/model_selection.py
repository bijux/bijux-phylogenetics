from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.regression import (
    compare_comparative_regression_models,
    write_comparative_regression_excluded_taxa_table,
    write_comparative_regression_model_ranking_table,
    write_comparative_regression_pairwise_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_regression_model_selection_command(comparative_subparsers: Any) -> None:
    comparative_model_selection = comparative_subparsers.add_parser(
        "model-selection",
        help="Rank competing comparative regression formulas on one shared taxon set.",
    )
    comparative_model_selection.add_argument("tree", type=Path)
    comparative_model_selection.add_argument("table", type=Path)
    comparative_model_selection.add_argument(
        "--formula",
        action="append",
        required=True,
        dest="formulas",
        help="Candidate comparative formula. Repeat this option once per model.",
    )
    comparative_model_selection.add_argument("--taxon-column")
    comparative_model_selection.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1. Binary-response model selection requires a numeric value.",
    )
    comparative_model_selection.add_argument(
        "--ranking-out",
        type=Path,
        help="Write the ranked comparative model table as TSV or CSV.",
    )
    comparative_model_selection.add_argument(
        "--pairwise-out",
        type=Path,
        help="Write nested-versus-non-nested pairwise comparison rows as TSV or CSV.",
    )
    comparative_model_selection.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the shared-complete-case excluded-taxa ledger as TSV or CSV.",
    )
    comparative_model_selection.add_argument(
        "--json", action="store_true", help="Emit the model-selection report as JSON."
    )
    _add_manifest_argument(comparative_model_selection)


def run_regression_model_selection_command(
    args: Any,
    *,
    lambda_value: float | str,
) -> int | None:
    if args.comparative_command != "model-selection":
        return None

    report = compare_comparative_regression_models(
        args.tree,
        args.table,
        formulas=list(args.formulas),
        taxon_column=args.taxon_column,
        lambda_value=lambda_value,
    )
    outputs: list[Path | str] = []
    if args.ranking_out is not None:
        outputs.append(
            write_comparative_regression_model_ranking_table(
                args.ranking_out,
                report,
            )
        )
    if args.pairwise_out is not None:
        outputs.append(
            write_comparative_regression_pairwise_table(
                args.pairwise_out,
                report,
            )
        )
    if args.excluded_taxa_out is not None:
        outputs.append(
            write_comparative_regression_excluded_taxa_table(
                args.excluded_taxa_out,
                report,
            )
        )
    outputs = _finalize_outputs(
        args,
        command="comparative",
        inputs=[args.tree, args.table],
        outputs=outputs,
    )
    selected_row = next(row for row in report.rows if row.selected)
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            warnings=report.warnings,
            metrics={
                "response": report.response,
                "model_family": report.model_family,
                "model_count": len(report.rows),
                "analysis_taxon_count": len(report.analysis_taxa),
                "excluded_taxon_count": len(report.excluded_taxa),
                "pairwise_comparison_count": len(report.pairwise_rows),
                "best_formula": report.best_formula,
                "selected_criterion": report.selected_criterion,
                "selected_log_likelihood": selected_row.log_likelihood,
            },
            data=report,
        ),
        json_output=args.json,
    )
    return 0

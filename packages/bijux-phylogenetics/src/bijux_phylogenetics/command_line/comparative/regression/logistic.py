from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.regression import (
    summarize_phylogenetic_logistic,
    write_phylogenetic_logistic_coefficient_table,
    write_phylogenetic_logistic_excluded_taxa_table,
    write_phylogenetic_logistic_fitted_table,
)
from bijux_phylogenetics.evidence.provenance.method_tiers import (
    method_tier_metrics,
    method_tier_warnings,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_logistic_command(comparative_subparsers: Any) -> None:
    comparative_logistic = comparative_subparsers.add_parser(
        "logistic",
        help="Fit a binary phylogenetic logistic approximation with a phylogenetic working correlation.",
    )
    comparative_logistic.add_argument("tree", type=Path)
    comparative_logistic.add_argument("table", type=Path)
    comparative_logistic.add_argument("--response")
    comparative_logistic.add_argument("--predictors", nargs="+")
    comparative_logistic.add_argument(
        "--formula",
        help="Formula-style specification such as 'presence ~ body_mass + habitat'.",
    )
    comparative_logistic.add_argument("--taxon-column")
    comparative_logistic.add_argument(
        "--lambda-value",
        default="1.0",
        help="Use a numeric Pagel lambda value between 0 and 1 for the working correlation.",
    )
    comparative_logistic.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the fitted logistic coefficient ledger as TSV or CSV.",
    )
    comparative_logistic.add_argument(
        "--fitted-out",
        type=Path,
        help="Write the fitted taxon-level probability ledger as TSV or CSV.",
    )
    comparative_logistic.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the excluded-taxa ledger as TSV or CSV.",
    )
    comparative_logistic.add_argument(
        "--json",
        action="store_true",
        help="Emit the logistic result as JSON.",
    )
    _add_manifest_argument(comparative_logistic)


def run_comparative_logistic_command(args: Any) -> int | None:
    if args.comparative_command != "logistic":
        return None

    report = summarize_phylogenetic_logistic(
        args.tree,
        args.table,
        response=args.response,
        predictors=list(args.predictors or []),
        formula=args.formula,
        taxon_column=args.taxon_column,
        lambda_value=float(args.lambda_value),
    )
    outputs: list[Path | str] = []
    if args.coefficients_out is not None:
        outputs.append(
            write_phylogenetic_logistic_coefficient_table(
                args.coefficients_out,
                report,
            )
        )
    if args.fitted_out is not None:
        outputs.append(
            write_phylogenetic_logistic_fitted_table(
                args.fitted_out,
                report,
            )
        )
    if args.excluded_taxa_out is not None:
        outputs.append(
            write_phylogenetic_logistic_excluded_taxa_table(
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
    _print_result(
        build_command_result(
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=outputs,
            metrics={
                "taxon_count": report.taxon_count,
                "success_count": report.success_count,
                "failure_count": report.failure_count,
                "coefficient_count": len(report.coefficients),
                "fitted_row_count": len(report.fitted_rows),
                "lambda_value": report.lambda_value,
                "approximation_method": report.approximation_method,
                "converged": report.converged,
                "iteration_count": report.iteration_count,
                "binomial_log_likelihood": report.binomial_log_likelihood,
                "separation_detected": report.separation_detected,
                "warning_count": len(report.warnings)
                + len(method_tier_warnings(report.method_tier)),
                "coefficient_inference_distribution": (
                    report.coefficients[0].inference_distribution
                    if report.coefficients
                    else None
                ),
                **method_tier_metrics(report.method_tier),
            },
            warnings=method_tier_warnings(report.method_tier)
            + [warning.message for warning in report.warnings],
            data=report,
        ),
        json_output=args.json,
    )
    return 0

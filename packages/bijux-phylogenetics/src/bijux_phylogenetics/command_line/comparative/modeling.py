from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.regression import (
    run_multivariate_comparative_regression,
    write_multivariate_excluded_taxa_table,
    write_multivariate_residual_association_table,
    write_multivariate_residual_correlation_table,
    write_multivariate_residual_covariance_table,
    write_multivariate_response_coefficient_table,
    write_multivariate_response_model_table,
)
from bijux_phylogenetics.comparative.regression_model_selection import (
    compare_comparative_regression_models,
    write_comparative_regression_excluded_taxa_table,
    write_comparative_regression_model_ranking_table,
    write_comparative_regression_pairwise_table,
)
from bijux_phylogenetics.comparative.reporting.analysis_package import (
    build_comparative_report_package,
)
from bijux_phylogenetics.comparative.reporting import (
    build_comparative_method_report,
    build_trait_influence_report,
    compare_comparative_results_across_pruning,
    compare_comparative_results_across_trees,
    write_comparative_method_report,
    write_comparative_methods_summary_text,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_modeling_commands(comparative_subparsers: Any) -> None:
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

    comparative_multivariate = comparative_subparsers.add_parser(
        "multivariate",
        help="Fit shared-taxon comparative regressions across multiple response traits.",
    )
    comparative_multivariate.add_argument("tree", type=Path)
    comparative_multivariate.add_argument("table", type=Path)
    comparative_multivariate.add_argument("--responses", nargs="+", required=True)
    comparative_multivariate.add_argument("--predictors", nargs="+", required=True)
    comparative_multivariate.add_argument("--taxon-column")
    comparative_multivariate.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multivariate.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the residual covariance ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--correlation-out",
        type=Path,
        help="Write the residual correlation ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--associations-out",
        type=Path,
        help="Write the residual trait-association ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--coefficients-out",
        type=Path,
        help="Write the per-response coefficient ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--response-models-out",
        type=Path,
        help="Write the per-response model summary ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write the explicit excluded-taxa ledger as TSV or CSV.",
    )
    comparative_multivariate.add_argument(
        "--json",
        action="store_true",
        help="Emit the multivariate regression report as JSON.",
    )
    _add_manifest_argument(comparative_multivariate)

    comparative_report = comparative_subparsers.add_parser(
        "report",
        help="Build an integrated comparative-method report.",
    )
    comparative_report.add_argument("tree", type=Path)
    comparative_report.add_argument("table", type=Path)
    comparative_report.add_argument("--response")
    comparative_report.add_argument("--predictors", nargs="+")
    comparative_report.add_argument("--formula")
    comparative_report.add_argument("--taxon-column")
    comparative_report.add_argument("--out", type=Path)
    comparative_report.add_argument(
        "--methods-summary-out",
        type=Path,
        help="Write reviewer-facing Markdown methods text for the comparative analysis.",
    )
    comparative_report.add_argument(
        "--out-dir",
        type=Path,
        help="Write a full comparative analysis package directory with HTML and reviewer TSV ledgers.",
    )
    comparative_report.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_report.add_argument(
        "--json", action="store_true", help="Emit the comparative report as JSON."
    )
    _add_manifest_argument(comparative_report)

    comparative_influence = comparative_subparsers.add_parser(
        "influence",
        help="Identify predictor terms and taxa driving one comparative result.",
    )
    comparative_influence.add_argument("tree", type=Path)
    comparative_influence.add_argument("table", type=Path)
    comparative_influence.add_argument("--response")
    comparative_influence.add_argument("--predictors", nargs="+")
    comparative_influence.add_argument("--formula")
    comparative_influence.add_argument("--taxon-column")
    comparative_influence.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_influence.add_argument(
        "--json", action="store_true", help="Emit the influence report as JSON."
    )
    _add_manifest_argument(comparative_influence)

    comparative_compare_trees = comparative_subparsers.add_parser(
        "compare-trees",
        help="Compare comparative results across two alternative trees.",
    )
    comparative_compare_trees.add_argument("left_tree", type=Path)
    comparative_compare_trees.add_argument("right_tree", type=Path)
    comparative_compare_trees.add_argument("table", type=Path)
    comparative_compare_trees.add_argument("--response")
    comparative_compare_trees.add_argument("--predictors", nargs="+")
    comparative_compare_trees.add_argument("--formula")
    comparative_compare_trees.add_argument("--taxon-column")
    comparative_compare_trees.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_trees.add_argument(
        "--json", action="store_true", help="Emit the comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_trees)

    comparative_compare_pruning = comparative_subparsers.add_parser(
        "compare-pruning",
        help="Compare comparative results before and after explicit pruning.",
    )
    comparative_compare_pruning.add_argument("tree", type=Path)
    comparative_compare_pruning.add_argument("table", type=Path)
    comparative_compare_pruning.add_argument("--response")
    comparative_compare_pruning.add_argument("--predictors", nargs="+")
    comparative_compare_pruning.add_argument("--formula")
    comparative_compare_pruning.add_argument("--drop-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--keep-taxa", nargs="+")
    comparative_compare_pruning.add_argument("--taxon-column")
    comparative_compare_pruning.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_compare_pruning.add_argument(
        "--json", action="store_true", help="Emit the pruning comparison as JSON."
    )
    _add_manifest_argument(comparative_compare_pruning)


def run_comparative_modeling_command(
    args: Any,
    *,
    parser: Any,
) -> int | None:
    del parser
    lambda_value: float | str
    if hasattr(args, "lambda_value"):
        if args.lambda_value == "estimate":
            lambda_value = "estimate"
        else:
            lambda_value = float(args.lambda_value)
    else:
        lambda_value = "estimate"
    if args.comparative_command == "model-selection":
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
    if args.comparative_command == "multivariate":
        report = run_multivariate_comparative_regression(
            args.tree,
            args.table,
            responses=list(args.responses),
            predictors=list(args.predictors),
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs: list[Path | str] = []
        if args.covariance_out is not None:
            outputs.append(
                write_multivariate_residual_covariance_table(args.covariance_out, report)
            )
        if args.correlation_out is not None:
            outputs.append(
                write_multivariate_residual_correlation_table(
                    args.correlation_out,
                    report,
                )
            )
        if args.associations_out is not None:
            outputs.append(
                write_multivariate_residual_association_table(
                    args.associations_out,
                    report,
                )
            )
        if args.coefficients_out is not None:
            outputs.append(
                write_multivariate_response_coefficient_table(
                    args.coefficients_out,
                    report,
                )
            )
        if args.response_models_out is not None:
            outputs.append(
                write_multivariate_response_model_table(
                    args.response_models_out,
                    report,
                )
            )
        if args.excluded_taxa_out is not None:
            outputs.append(
                write_multivariate_excluded_taxa_table(
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
                    "response_count": len(report.responses),
                    "predictor_count": len(report.predictors),
                    "analysis_taxa": len(report.analysis_taxa),
                    "excluded_taxa": len(report.excluded_taxa),
                    "residual_covariance_response_count": (
                        report.covariance_diagnostics.response_count
                    ),
                    "residual_covariance_matrix_rank": (
                        report.covariance_diagnostics.matrix_rank
                    ),
                    "residual_covariance_condition_number": (
                        report.covariance_diagnostics.condition_number
                    ),
                    "residual_covariance_singular": (
                        report.covariance_diagnostics.is_singular
                    ),
                    "residual_covariance_near_singular": (
                        report.covariance_diagnostics.is_near_singular
                    ),
                    "residual_covariance_row_count": len(report.covariance_rows),
                    "residual_correlation_row_count": len(report.correlation_rows),
                    "residual_association_count": len(report.association_rows),
                    "response_model_count": len(report.response_model_rows),
                    "coefficient_row_count": len(report.coefficient_rows),
                    "warning_count": len(report.warnings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "report":
        package_result = None
        if args.out_dir is not None:
            package_result = build_comparative_report_package(
                args.tree,
                args.table,
                out_dir=args.out_dir,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
            report = package_result.report
        else:
            report = build_comparative_method_report(
                args.tree,
                args.table,
                response=args.response,
                predictors=list(args.predictors or []),
                formula=args.formula,
                taxon_column=args.taxon_column,
                lambda_value=lambda_value,
            )
        output_paths: list[Path | str] = []
        if args.out is not None:
            write_comparative_method_report(args.out, report)
            output_paths.append(args.out)
        if args.methods_summary_out is not None:
            write_comparative_methods_summary_text(args.methods_summary_out, report)
            output_paths.append(args.methods_summary_out)
        if package_result is not None:
            output_paths.extend(
                [
                    package_result.report_path,
                    package_result.methods_summary_path,
                    package_result.reviewer_audit_checklist_path,
                    package_result.summary_table_path,
                    package_result.coefficient_table_path,
                    package_result.residual_table_path,
                    package_result.signal_table_path,
                    package_result.model_comparison_table_path,
                    package_result.interpretation_table_path,
                    package_result.audit_table_path,
                    package_result.contrast_table_path,
                    package_result.manifest_path,
                ]
            )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.tree, args.table],
            outputs=output_paths,
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.tree, args.table],
                outputs=outputs,
                metrics={
                    "taxon_count": report.snapshot.pgls_model.taxon_count,
                    "selected_model": report.snapshot.model_comparison.better_model,
                    "audit_row_count": len(report.snapshot.audit_rows),
                    "excluded_taxa": len(
                        report.snapshot.pgls_inputs.formula_audit.excluded_taxa
                    ),
                    "limitation_count": len(report.snapshot.limitations),
                    "coefficient_count": len(report.snapshot.pgls_model.coefficients),
                    "package_output_count": len(outputs),
                },
                data=report if package_result is None else package_result,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "influence":
        report = build_trait_influence_report(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
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
                    "predictor_count": len(report.predictor_rows),
                    "taxon_count": len(report.taxon_rows),
                    "top_predictor_terms": len(report.top_predictor_terms),
                    "top_taxa": len(report.top_taxa),
                    "selected_model": report.selected_model,
                },
                warnings=report.warnings,
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "compare-trees":
        report = compare_comparative_results_across_trees(
            args.left_tree,
            args.right_tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs = _finalize_outputs(
            args,
            command="comparative",
            inputs=[args.left_tree, args.right_tree, args.table],
        )
        _print_result(
            build_command_result(
                command="comparative",
                inputs=[args.left_tree, args.right_tree, args.table],
                outputs=outputs,
                metrics={
                    "coefficient_delta_count": len(report.coefficient_deltas),
                    "sign_changed_terms": len(report.sign_changed_terms),
                    "conclusion_changed": report.conclusion_changed,
                    "left_selected_model": report.left_selected_model,
                    "right_selected_model": report.right_selected_model,
                },
                warnings=report.warnings,
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "compare-pruning":
        report = compare_comparative_results_across_pruning(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            drop_taxa=list(args.drop_taxa or []),
            keep_taxa=list(args.keep_taxa or []),
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
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
                    "baseline_taxa": len(report.baseline_taxa),
                    "pruned_taxa": len(report.pruned_taxa),
                    "dropped_taxa": len(report.dropped_taxa),
                    "sign_changed_terms": len(report.sign_changed_terms),
                    "conclusion_changed": report.conclusion_changed,
                },
                warnings=report.warnings,
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    return None

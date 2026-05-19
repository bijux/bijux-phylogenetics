from __future__ import annotations

from pathlib import Path
from typing import Any

from bijux_phylogenetics.command_line.arguments import _add_manifest_argument
from bijux_phylogenetics.command_line.output import _print_result
from bijux_phylogenetics.command_line.routing import _finalize_outputs
from bijux_phylogenetics.comparative.covariance_audit import (
    summarize_comparative_covariance_audit,
    write_comparative_covariance_audit_candidate_table,
    write_comparative_covariance_audit_excluded_taxa_table,
    write_comparative_covariance_audit_summary_table,
)
from bijux_phylogenetics.comparative.pgls import (
    inspect_pgls_inputs,
    run_pgls,
    run_pgls_multiple_testing,
    write_pgls_model_matrix_table,
)
from bijux_phylogenetics.comparative.pgls_brownian_covariance import (
    summarize_brownian_covariance_pgls,
    write_brownian_covariance_table,
)
from bijux_phylogenetics.comparative.pgls_categorical_contrasts import (
    summarize_pgls_categorical_contrasts,
    write_pgls_categorical_contrast_table,
)
from bijux_phylogenetics.comparative.pgls_interaction_coefficients import (
    summarize_pgls_interaction_coefficients,
    write_pgls_interaction_coefficient_table,
)
from bijux_phylogenetics.comparative.pgls_lambda_fit import (
    write_pgls_lambda_profile_table,
)
from bijux_phylogenetics.comparative.pgls_ou_covariance import (
    summarize_ou_covariance_pgls,
    write_ou_alpha_profile_table,
    write_ou_covariance_table,
)
from bijux_phylogenetics.runtime.results import build_command_result


def add_comparative_pgls_commands(comparative_subparsers: Any) -> None:
    comparative_covariance_audit = comparative_subparsers.add_parser(
        "covariance-audit",
        help="Inspect the covariance matrix shape and stability for supported comparative analyses.",
    )
    comparative_covariance_audit.add_argument("tree", type=Path)
    comparative_covariance_audit.add_argument("table", type=Path)
    comparative_covariance_audit.add_argument(
        "--analysis",
        choices=("pgls", "brownian-trait", "ou-trait"),
        required=True,
    )
    comparative_covariance_audit.add_argument("--trait")
    comparative_covariance_audit.add_argument("--response")
    comparative_covariance_audit.add_argument("--predictors", nargs="+")
    comparative_covariance_audit.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_covariance_audit.add_argument("--taxon-column")
    comparative_covariance_audit.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_covariance_audit.add_argument(
        "--alpha",
        default="estimate",
        help="Use 'estimate' or a positive OU alpha value.",
    )
    comparative_covariance_audit.add_argument(
        "--summary-out",
        type=Path,
        help="Write one covariance-audit summary row as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--candidates-out",
        type=Path,
        help="Write candidate covariance audit rows as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--excluded-taxa-out",
        type=Path,
        help="Write excluded covariance-audit taxa as TSV or CSV.",
    )
    comparative_covariance_audit.add_argument(
        "--json", action="store_true", help="Emit the covariance audit as JSON."
    )
    _add_manifest_argument(comparative_covariance_audit)

    comparative_pgls = comparative_subparsers.add_parser(
        "pgls",
        help="Fit a phylogenetic generalized least-squares model.",
    )
    comparative_pgls.add_argument("tree", type=Path)
    comparative_pgls.add_argument("table", type=Path)
    comparative_pgls.add_argument("--response")
    comparative_pgls.add_argument("--predictors", nargs="+")
    comparative_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_pgls.add_argument("--taxon-column")
    comparative_pgls.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_pgls.add_argument(
        "--model-matrix-out",
        type=Path,
        help="Write the encoded comparative model matrix as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--categorical-contrasts-out",
        type=Path,
        help="Write categorical predictor contrast rows as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--interaction-coefficients-out",
        type=Path,
        help="Write interaction coefficient rows as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--lambda-profile-out",
        type=Path,
        help="Write the fitted Pagel lambda likelihood profile as TSV or CSV.",
    )
    comparative_pgls.add_argument(
        "--json", action="store_true", help="Emit the model result as JSON."
    )
    _add_manifest_argument(comparative_pgls)

    comparative_brownian_pgls = comparative_subparsers.add_parser(
        "brownian-pgls",
        help="Fit a PGLS model under fixed Brownian shared-path covariance.",
    )
    comparative_brownian_pgls.add_argument("tree", type=Path)
    comparative_brownian_pgls.add_argument("table", type=Path)
    comparative_brownian_pgls.add_argument("--response")
    comparative_brownian_pgls.add_argument("--predictors", nargs="+")
    comparative_brownian_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_brownian_pgls.add_argument("--taxon-column")
    comparative_brownian_pgls.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the pairwise Brownian covariance ledger as TSV or CSV.",
    )
    comparative_brownian_pgls.add_argument(
        "--json", action="store_true", help="Emit the Brownian PGLS result as JSON."
    )
    _add_manifest_argument(comparative_brownian_pgls)

    comparative_ou_pgls = comparative_subparsers.add_parser(
        "ou-pgls",
        help="Fit a PGLS model under stationary-root OU covariance.",
    )
    comparative_ou_pgls.add_argument("tree", type=Path)
    comparative_ou_pgls.add_argument("table", type=Path)
    comparative_ou_pgls.add_argument("--response")
    comparative_ou_pgls.add_argument("--predictors", nargs="+")
    comparative_ou_pgls.add_argument(
        "--formula",
        help="Formula-style specification such as 'response ~ body_mass * habitat'.",
    )
    comparative_ou_pgls.add_argument("--taxon-column")
    comparative_ou_pgls.add_argument(
        "--alpha",
        default="estimate",
        help="Use 'estimate' or a positive numeric OU alpha value.",
    )
    comparative_ou_pgls.add_argument(
        "--covariance-out",
        type=Path,
        help="Write the pairwise OU covariance ledger as TSV or CSV.",
    )
    comparative_ou_pgls.add_argument(
        "--alpha-profile-out",
        type=Path,
        help="Write the fitted OU alpha likelihood profile as TSV or CSV.",
    )
    comparative_ou_pgls.add_argument(
        "--json",
        action="store_true",
        help="Emit the OU covariance PGLS result as JSON.",
    )
    _add_manifest_argument(comparative_ou_pgls)

    comparative_multiple_testing = comparative_subparsers.add_parser(
        "multiple-testing",
        help="Adjust PGLS coefficient p-values across many response traits.",
    )
    comparative_multiple_testing.add_argument("tree", type=Path)
    comparative_multiple_testing.add_argument("table", type=Path)
    comparative_multiple_testing.add_argument("--responses", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--predictors", nargs="+", required=True)
    comparative_multiple_testing.add_argument("--taxon-column")
    comparative_multiple_testing.add_argument(
        "--lambda-value",
        default="estimate",
        help="Use 'estimate' or a numeric Pagel lambda value between 0 and 1.",
    )
    comparative_multiple_testing.add_argument(
        "--json", action="store_true", help="Emit the correction report as JSON."
    )
    _add_manifest_argument(comparative_multiple_testing)


def run_comparative_pgls_command(
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
    if args.comparative_command == "covariance-audit":
        resolved_lambda: float | str
        if args.lambda_value == "estimate":
            resolved_lambda = "estimate"
        else:
            resolved_lambda = float(args.lambda_value)
        resolved_alpha: float | str
        if args.alpha == "estimate":
            resolved_alpha = "estimate"
        else:
            resolved_alpha = float(args.alpha)
        report = summarize_comparative_covariance_audit(
            args.tree,
            args.table,
            analysis=args.analysis,
            trait=args.trait,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=resolved_lambda,
            alpha=resolved_alpha,
        )
        outputs: list[Path | str] = []
        if args.summary_out is not None:
            outputs.append(
                write_comparative_covariance_audit_summary_table(
                    args.summary_out,
                    report,
                )
            )
        if args.candidates_out is not None:
            outputs.append(
                write_comparative_covariance_audit_candidate_table(
                    args.candidates_out,
                    report,
                )
            )
        if args.excluded_taxa_out is not None:
            outputs.append(
                write_comparative_covariance_audit_excluded_taxa_table(
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
                warnings=report.warnings,
                metrics={
                    "analysis": report.analysis,
                    "covariance_model": report.covariance_model,
                    "matrix_dimension": report.matrix_dimension,
                    "matrix_rank": report.matrix_rank,
                    "condition_number": report.condition_number,
                    "fit_strategy": report.fit_strategy,
                    "singular": report.singular,
                    "near_singular": report.near_singular,
                    "matched_taxon_count": len(report.matched_taxa),
                    "analysis_taxon_count": len(report.analysis_taxa),
                    "duplicate_tree_taxon_count": len(report.duplicate_tree_taxa),
                    "duplicate_trait_taxon_count": len(report.duplicate_trait_taxa),
                    "candidate_row_count": len(report.candidate_rows),
                    "blocker_count": len(report.blockers),
                    "warning_count": len(report.warnings),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "multiple-testing":
        report = run_pgls_multiple_testing(
            args.tree,
            args.table,
            responses=list(args.responses),
            predictors=list(args.predictors),
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
                    "response_count": len(report.responses),
                    "test_count": len(report.rows),
                    "family_size": report.family_size,
                    "raw_significant_count": report.raw_significant_count,
                    "significant_count": sum(1 for row in report.rows if row.significant),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "brownian-pgls":
        report = summarize_brownian_covariance_pgls(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
        )
        outputs: list[Path | str] = []
        if args.covariance_out is not None:
            outputs.append(
                write_brownian_covariance_table(
                    args.covariance_out,
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
                    "predictor_count": len(report.model.predictors),
                    "coefficient_count": len(report.model.coefficients),
                    "covariance_row_count": len(report.rows),
                    "lambda_value": report.model.lambda_value,
                    "covariance_model": "brownian-shared-path",
                    "tree_is_ultrametric": report.tree_is_ultrametric,
                    "minimum_root_to_tip_depth": report.minimum_root_to_tip_depth,
                    "maximum_root_to_tip_depth": report.maximum_root_to_tip_depth,
                    "raw_log_determinant": report.raw_log_determinant,
                    "positive_definite_before_stabilization": (
                        report.positive_definite_before_stabilization
                    ),
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "ou-pgls":
        report = summarize_ou_covariance_pgls(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            alpha=args.alpha,
        )
        outputs: list[Path | str] = []
        if args.covariance_out is not None:
            outputs.append(write_ou_covariance_table(args.covariance_out, report))
        if args.alpha_profile_out is not None:
            outputs.append(
                write_ou_alpha_profile_table(
                    args.alpha_profile_out,
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
                    "predictor_count": len(report.model.predictors),
                    "coefficient_count": len(report.model.coefficients),
                    "covariance_row_count": len(report.rows),
                    "alpha": report.alpha,
                    "alpha_estimation_mode": report.alpha_estimation_mode,
                    "alpha_profile_point_count": len(report.alpha_profile_rows),
                    "alpha_lower_95_confidence_interval": (
                        report.lower_95_confidence_interval
                    ),
                    "alpha_upper_95_confidence_interval": (
                        report.upper_95_confidence_interval
                    ),
                    "covariance_model": "ou-stationary-root",
                    "tree_is_ultrametric": report.tree_is_ultrametric,
                    "raw_log_determinant": report.raw_log_determinant,
                    "positive_definite_before_stabilization": (
                        report.positive_definite_before_stabilization
                    ),
                    "log_likelihood": report.model.log_likelihood,
                    "aic": report.model.aic,
                },
                data=report,
            ),
            json_output=args.json,
        )
        return 0
    if args.comparative_command == "pgls":
        input_report = inspect_pgls_inputs(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
        )
        report = run_pgls(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        categorical_contrasts = summarize_pgls_categorical_contrasts(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        interaction_coefficients = summarize_pgls_interaction_coefficients(
            args.tree,
            args.table,
            response=args.response,
            predictors=list(args.predictors or []),
            formula=args.formula,
            taxon_column=args.taxon_column,
            lambda_value=lambda_value,
        )
        outputs: list[Path | str] = []
        if args.model_matrix_out is not None:
            outputs.append(
                write_pgls_model_matrix_table(
                    args.model_matrix_out,
                    input_report.model_matrix,
                )
            )
        if args.categorical_contrasts_out is not None:
            outputs.append(
                write_pgls_categorical_contrast_table(
                    args.categorical_contrasts_out,
                    categorical_contrasts,
                )
            )
        if args.interaction_coefficients_out is not None:
            outputs.append(
                write_pgls_interaction_coefficient_table(
                    args.interaction_coefficients_out,
                    interaction_coefficients,
                )
            )
        if args.lambda_profile_out is not None:
            outputs.append(
                write_pgls_lambda_profile_table(
                    args.lambda_profile_out,
                    report.lambda_fit,
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
                warnings=input_report.warnings,
                metrics={
                    "taxon_count": report.taxon_count,
                    "predictor_count": len(report.predictors),
                    "coefficient_count": len(report.coefficients),
                    "confidence_interval_count": len(report.coefficients),
                    "categorical_contrast_predictor_count": (
                        categorical_contrasts.categorical_predictor_count
                    ),
                    "categorical_contrast_row_count": len(categorical_contrasts.rows),
                    "interaction_term_count": (
                        interaction_coefficients.interaction_term_count
                    ),
                    "interaction_coefficient_row_count": len(
                        interaction_coefficients.rows
                    ),
                    "intercept_included": input_report.formula.include_intercept,
                    "model_matrix_row_count": input_report.model_matrix.row_count,
                    "model_matrix_column_count": len(
                        input_report.model_matrix.encoded_columns
                    ),
                    "residual_degrees_of_freedom": (
                        report.coefficients[0].degrees_of_freedom
                        if report.coefficients
                        else 0
                    ),
                    "coefficient_inference_distribution": (
                        report.coefficients[0].inference_distribution
                        if report.coefficients
                        else None
                    ),
                    "encoded_predictor_count": len(
                        input_report.model_matrix.encoded_columns
                    )
                    - (1 if input_report.formula.include_intercept else 0),
                    "categorical_predictor_count": len(
                        input_report.categorical_predictors
                    ),
                    "transformed_term_count": len(
                        input_report.formula_audit.transformed_terms
                    ),
                    "lambda_value": report.lambda_value,
                    "lambda_estimation_mode": report.lambda_fit.mode,
                    "lambda_profile_point_count": len(report.lambda_fit.profile_rows),
                    "lambda_lower_95_confidence_interval": (
                        report.lambda_fit.lower_95_confidence_interval
                    ),
                    "lambda_upper_95_confidence_interval": (
                        report.lambda_fit.upper_95_confidence_interval
                    ),
                    "aic": report.aic,
                    "r_squared": report.r_squared,
                },
                data={
                    "inputs": input_report,
                    "model": report,
                    "categorical_contrasts": categorical_contrasts,
                    "interaction_coefficients": interaction_coefficients,
                },
            ),
            json_output=args.json,
        )
        return 0
    return None

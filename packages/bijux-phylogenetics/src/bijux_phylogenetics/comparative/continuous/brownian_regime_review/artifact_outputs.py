from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from .contracts import BrownianRegimeFitSummaryReport


def write_brownian_regime_summary_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "branch_id_column",
            "regime_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "regime_count",
            "root_state",
            "root_state_lower_95",
            "root_state_upper_95",
            "log_likelihood",
            "aic",
            "aicc",
            "better_model",
            "likelihood_ratio_statistic",
            "likelihood_ratio_degrees_of_freedom",
            "likelihood_ratio_p_value",
            "identifiability_warning_count",
            "residual_variance",
            "max_abs_standardized_residual",
            "phylogenetic_residual_lambda",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "branch_id_column": report.branch_id_column,
                "regime_column": report.regime_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "regime_count": len(report.regime_rows),
                "root_state": format(report.root_state, ".15g"),
                "root_state_lower_95": format(
                    report.root_state_interval.lower_95,
                    ".15g",
                ),
                "root_state_upper_95": format(
                    report.root_state_interval.upper_95,
                    ".15g",
                ),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "aic": format(report.aic, ".15g"),
                "aicc": format(report.aicc, ".15g"),
                "better_model": report.better_model,
                "likelihood_ratio_statistic": format(
                    report.likelihood_ratio_statistic,
                    ".15g",
                ),
                "likelihood_ratio_degrees_of_freedom": report.likelihood_ratio_degrees_of_freedom,
                "likelihood_ratio_p_value": format(
                    report.likelihood_ratio_p_value,
                    ".15g",
                ),
                "identifiability_warning_count": len(report.identifiability_warnings),
                "residual_variance": format(
                    report.residual_diagnostics.residual_variance,
                    ".15g",
                ),
                "max_abs_standardized_residual": format(
                    report.residual_diagnostics.max_abs_standardized_residual,
                    ".15g",
                ),
                "phylogenetic_residual_lambda": format(
                    report.residual_diagnostics.phylogenetic_residual_lambda,
                    ".15g",
                ),
            }
        ],
    )


def write_brownian_regime_comparison_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    best_aicc = min(row.aicc for row in report.comparison_rows)
    rows = [
        {
            "row_kind": "model_fit",
            "model": row.model,
            "comparison_id": "",
            "parameter_count": row.parameter_count,
            "log_likelihood": format(row.log_likelihood, ".15g"),
            "aic": format(row.aic, ".15g"),
            "aicc": format(row.aicc, ".15g"),
            "delta_aicc": format(row.aicc - best_aicc, ".15g"),
            "selected": str(row.selected).lower(),
            "left_model": "",
            "right_model": "",
            "statistic": "",
            "degrees_of_freedom": "",
            "p_value": "",
            "p_value_method": "",
        }
        for row in report.comparison_rows
    ]
    rows.append(
        {
            "row_kind": "likelihood_ratio_test",
            "model": "",
            "comparison_id": "brownian-vs-brownian-regimes",
            "parameter_count": "",
            "log_likelihood": "",
            "aic": "",
            "aicc": "",
            "delta_aicc": "",
            "selected": "",
            "left_model": "brownian",
            "right_model": "brownian-regimes",
            "statistic": format(report.likelihood_ratio_statistic, ".15g"),
            "degrees_of_freedom": report.likelihood_ratio_degrees_of_freedom,
            "p_value": format(report.likelihood_ratio_p_value, ".15g"),
            "p_value_method": report.likelihood_ratio_p_value_method,
        }
    )
    return write_taxon_rows(
        path,
        columns=[
            "row_kind",
            "model",
            "comparison_id",
            "parameter_count",
            "log_likelihood",
            "aic",
            "aicc",
            "delta_aicc",
            "selected",
            "left_model",
            "right_model",
            "statistic",
            "degrees_of_freedom",
            "p_value",
            "p_value_method",
        ],
        rows=rows,
    )


def write_brownian_regime_rate_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "regime",
            "branch_count",
            "contributing_branch_count",
            "total_branch_length",
            "contributing_branch_length",
            "sigma_squared",
            "sigma_squared_lower_95",
            "sigma_squared_upper_95",
            "interval_method",
        ],
        rows=[
            {
                "regime": row.regime,
                "branch_count": row.branch_count,
                "contributing_branch_count": row.contributing_branch_count,
                "total_branch_length": format(row.total_branch_length, ".15g"),
                "contributing_branch_length": format(
                    row.contributing_branch_length,
                    ".15g",
                ),
                "sigma_squared": format(row.sigma_squared, ".15g"),
                "sigma_squared_lower_95": format(row.lower_95, ".15g"),
                "sigma_squared_upper_95": format(row.upper_95, ".15g"),
                "interval_method": row.interval_method,
            }
            for row in report.regime_rows
        ],
    )


def write_brownian_regime_profile_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "regime",
            "sigma_squared",
            "log_likelihood",
            "delta_log_likelihood",
            "in_support_interval",
            "selected",
        ],
        rows=[
            {
                "regime": row.regime,
                "sigma_squared": format(row.sigma_squared, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "delta_log_likelihood": format(row.delta_log_likelihood, ".15g"),
                "in_support_interval": str(row.in_support_interval).lower(),
                "selected": str(row.selected).lower(),
            }
            for row in report.profile_rows
        ],
    )


def write_brownian_regime_branch_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=[
            "branch_id",
            "regime",
            "branch_length",
            "descendant_taxa",
            "analyzed_descendant_taxa",
            "contributes_to_analysis",
        ],
        rows=[
            {
                "branch_id": row.branch_id,
                "regime": row.regime,
                "branch_length": format(row.branch_length, ".15g"),
                "descendant_taxa": ",".join(row.descendant_taxa),
                "analyzed_descendant_taxa": ",".join(row.analyzed_descendant_taxa),
                "contributes_to_analysis": str(row.contributes_to_analysis).lower(),
            }
            for row in report.branch_rows
        ],
    )


def write_brownian_regime_exclusion_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )

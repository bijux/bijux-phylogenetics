from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    fit_brownian_motion_model,
)
from bijux_phylogenetics.comparative.continuous.brownian_regime_review.contracts import (
    BrownianRegimeBranchRow,
    BrownianRegimeExclusion,
    BrownianRegimeFitSummaryReport,
    BrownianRegimeIdentifiabilityWarning,
    BrownianRegimeProfileRow,
    BrownianRegimeRateRow,
)
from bijux_phylogenetics.comparative.continuous.brownian_regime_review.covariance_fitting import (
    fit_multirate_brownian_model as _fit_multirate_brownian_model,
)
from bijux_phylogenetics.comparative.continuous.brownian_regime_review.regime_mapping import (
    build_excluded_taxa as _build_excluded_taxa,
    build_regime_covariance_components as _build_regime_covariance_components,
    load_branch_regime_rows as _load_branch_regime_rows,
)
from bijux_phylogenetics.comparative.continuous.brownian_regime_review.review_metrics import (
    build_comparison_rows as _build_comparison_rows,
    build_identifiability_warnings as _build_identifiability_warnings,
    build_multirate_residual_diagnostics as _build_multirate_residual_diagnostics,
    build_profile_rows as _build_profile_rows,
    build_regime_rows as _build_regime_rows,
    chi_square_survival as _chi_square_survival,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

__all__ = [
    "BrownianRegimeBranchRow",
    "BrownianRegimeExclusion",
    "BrownianRegimeFitSummaryReport",
    "BrownianRegimeIdentifiabilityWarning",
    "BrownianRegimeProfileRow",
    "BrownianRegimeRateRow",
    "summarize_brownian_regime_rates",
    "write_brownian_regime_branch_table",
    "write_brownian_regime_comparison_table",
    "write_brownian_regime_exclusion_table",
    "write_brownian_regime_profile_table",
    "write_brownian_regime_rate_table",
    "write_brownian_regime_summary_table",
]
def summarize_brownian_regime_rates(
    tree_path: Path,
    traits_path: Path,
    regime_map_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    branch_id_column: str | None = None,
    regime_column: str = "regime",
) -> BrownianRegimeFitSummaryReport:
    """Summarize a multi-rate Brownian model from a user-provided branch regime map."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    source_tree = load_tree(tree_path)
    branch_rows, resolved_branch_id_column = _load_branch_regime_rows(
        regime_map_path,
        tree=source_tree,
        analyzed_taxa=dataset.taxa,
        branch_id_column=branch_id_column,
        regime_column=regime_column,
    )
    regime_matrices = _build_regime_covariance_components(dataset.taxa, branch_rows)
    if len(regime_matrices) < 2:
        raise ComparativeMethodError(
            "multi-rate Brownian fitting requires at least two contributing regimes"
        )
    baseline = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    fit = _fit_multirate_brownian_model(
        dataset,
        regime_matrices=regime_matrices,
        baseline_rate=baseline.rate,
    )
    comparison_rows = _build_comparison_rows(
        taxon_count=dataset.taxon_count
        if hasattr(dataset, "taxon_count")
        else len(dataset.taxa),
        baseline_log_likelihood=baseline.log_likelihood,
        multirate_log_likelihood=fit.log_likelihood,
        regime_count=len(regime_matrices),
    )
    better_model = next(row.model for row in comparison_rows if row.selected)
    likelihood_ratio_degrees_of_freedom = len(regime_matrices) - 1
    likelihood_ratio_statistic = max(
        0.0,
        -2.0 * (baseline.log_likelihood - fit.log_likelihood),
    )
    likelihood_ratio_p_value = _chi_square_survival(
        likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom,
    )
    profile_rows = _build_profile_rows(
        dataset,
        regime_matrices=regime_matrices,
        fitted_rates=fit.regime_rates,
        baseline_rate=baseline.rate,
    )
    regime_rows = _build_regime_rows(branch_rows, fit.regime_rates, profile_rows)
    identifiability_warnings = _build_identifiability_warnings(
        regime_rows,
        profile_rows,
        better_model=better_model,
    )
    residual_diagnostics = _build_multirate_residual_diagnostics(
        dataset,
        fit.covariance,
        fit.residuals,
    )
    warnings = list(
        dict.fromkeys(
            [
                *readiness.warnings,
                *residual_diagnostics.warnings,
                *[warning.message for warning in identifiability_warnings],
            ]
        )
    )
    return BrownianRegimeFitSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        regime_map_path=regime_map_path,
        taxon_column=dataset.taxon_column,
        branch_id_column=resolved_branch_id_column,
        regime_column=regime_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(dataset.taxa),
        analyzed_taxon_count=len(dataset.taxa),
        excluded_taxa=_build_excluded_taxa(readiness),
        branch_rows=branch_rows,
        regime_rows=regime_rows,
        profile_rows=profile_rows,
        root_state=fit.root_state,
        root_state_interval=fit.root_state_interval,
        log_likelihood=fit.log_likelihood,
        aic=next(row.aic for row in comparison_rows if row.model == "brownian-regimes"),
        aicc=next(
            row.aicc for row in comparison_rows if row.model == "brownian-regimes"
        ),
        comparison_rows=comparison_rows,
        better_model=better_model,
        likelihood_ratio_statistic=likelihood_ratio_statistic,
        likelihood_ratio_degrees_of_freedom=likelihood_ratio_degrees_of_freedom,
        likelihood_ratio_p_value=likelihood_ratio_p_value,
        likelihood_ratio_p_value_method="wilson-hilferty-approximation",
        identifiability_warnings=identifiability_warnings,
        residual_diagnostics=residual_diagnostics,
        assumptions=[
            "Every non-root branch must be assigned to one user-provided regime.",
            "Each regime contributes its own Brownian sigma-squared rate to shared-path covariance.",
            "Regime-specific uncertainty is reported from conditional likelihood profiles with other regimes fixed at their best-supported values.",
        ],
        warnings=warnings,
        readiness=readiness,
    )


def write_brownian_regime_summary_table(
    path: Path,
    report: BrownianRegimeFitSummaryReport,
) -> Path:
    """Write one overall summary ledger for a multi-rate Brownian fit."""
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
                    report.root_state_interval.lower_95, ".15g"
                ),
                "root_state_upper_95": format(
                    report.root_state_interval.upper_95, ".15g"
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
    """Write one single-rate versus multi-rate Brownian comparison ledger."""
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
    """Write one per-regime rate ledger for a multi-rate Brownian fit."""
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
    """Write one conditional rate-profile ledger for a multi-rate Brownian fit."""
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
    """Write one normalized branch-regime assignment ledger."""
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
    """Write one excluded-taxon ledger for a multi-rate Brownian fit."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )

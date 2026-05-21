from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    matrix_multiply,
    transpose,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    load_comparative_dataset,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeResidualOutlier,
    ComparativeResidualSummary,
    _comparison_row,
    _estimate_lambda_for_values,
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
    fit_multirate_covariance as _fit_multirate_covariance,
    logspace as _logspace,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)
from bijux_phylogenetics.comparative.continuous.brownian_regime_review.regime_mapping import (
    build_excluded_taxa as _build_excluded_taxa,
    build_regime_covariance_components as _build_regime_covariance_components,
    load_branch_regime_rows as _load_branch_regime_rows,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

_PROFILE_CONFIDENCE_DELTA = 1.920729410347062

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


def _build_comparison_rows(
    *,
    taxon_count: int,
    baseline_log_likelihood: float,
    multirate_log_likelihood: float,
    regime_count: int,
) -> list[ComparativeModelComparisonRow]:
    rows = [
        _comparison_row("brownian", 2, baseline_log_likelihood, taxon_count),
        _comparison_row(
            "brownian-regimes",
            regime_count + 1,
            multirate_log_likelihood,
            taxon_count,
        ),
    ]
    best_aicc = min(row.aicc for row in rows)
    for row in rows:
        row.selected = math.isclose(row.aicc, best_aicc, rel_tol=0.0, abs_tol=1e-12)
    return rows


def _build_profile_rows(
    dataset: ComparativeDataset,
    *,
    regime_matrices: dict[str, list[list[float]]],
    fitted_rates: dict[str, float],
    baseline_rate: float,
) -> list[BrownianRegimeProfileRow]:
    lower = max(baseline_rate * 0.02, 1e-6)
    upper = max(baseline_rate * 50.0, lower * 10.0)
    rows: list[BrownianRegimeProfileRow] = []
    for regime in sorted(regime_matrices):
        candidates = sorted(
            {
                *_logspace(lower, upper, 81),
                fitted_rates[regime],
            }
        )
        profile = []
        best_log_likelihood = -math.inf
        for candidate in candidates:
            trial_rates = dict(fitted_rates)
            trial_rates[regime] = candidate
            fit = _fit_multirate_covariance(dataset, regime_matrices, trial_rates)
            profile.append((candidate, fit.log_likelihood))
            best_log_likelihood = max(best_log_likelihood, fit.log_likelihood)
        for candidate, log_likelihood in profile:
            rows.append(
                BrownianRegimeProfileRow(
                    regime=regime,
                    sigma_squared=candidate,
                    log_likelihood=log_likelihood,
                    delta_log_likelihood=best_log_likelihood - log_likelihood,
                    in_support_interval=(
                        (best_log_likelihood - log_likelihood)
                        <= _PROFILE_CONFIDENCE_DELTA
                    ),
                    selected=math.isclose(
                        candidate,
                        fitted_rates[regime],
                        rel_tol=0.0,
                        abs_tol=1e-12,
                    ),
                )
            )
    return rows


def _build_regime_rows(
    branch_rows: list[BrownianRegimeBranchRow],
    regime_rates: dict[str, float],
    profile_rows: list[BrownianRegimeProfileRow],
) -> list[BrownianRegimeRateRow]:
    by_regime: dict[str, list[BrownianRegimeBranchRow]] = {}
    for row in branch_rows:
        by_regime.setdefault(row.regime, []).append(row)
    profile_by_regime: dict[str, list[BrownianRegimeProfileRow]] = {}
    for row in profile_rows:
        profile_by_regime.setdefault(row.regime, []).append(row)
    rows: list[BrownianRegimeRateRow] = []
    for regime in sorted(regime_rates):
        branches = by_regime[regime]
        support = [
            row.sigma_squared
            for row in profile_by_regime.get(regime, [])
            if row.in_support_interval
        ]
        rows.append(
            BrownianRegimeRateRow(
                regime=regime,
                branch_count=len(branches),
                contributing_branch_count=sum(
                    1 for branch in branches if branch.contributes_to_analysis
                ),
                total_branch_length=sum(branch.branch_length for branch in branches),
                contributing_branch_length=sum(
                    branch.branch_length
                    for branch in branches
                    if branch.contributes_to_analysis
                ),
                sigma_squared=regime_rates[regime],
                lower_95=min(support) if support else regime_rates[regime],
                upper_95=max(support) if support else regime_rates[regime],
                interval_method="conditional_profile_likelihood_95",
            )
        )
    return rows


def _build_identifiability_warnings(
    regime_rows: list[BrownianRegimeRateRow],
    profile_rows: list[BrownianRegimeProfileRow],
    *,
    better_model: str,
) -> list[BrownianRegimeIdentifiabilityWarning]:
    warnings: list[BrownianRegimeIdentifiabilityWarning] = []
    profile_by_regime: dict[str, list[BrownianRegimeProfileRow]] = {}
    for row in profile_rows:
        profile_by_regime.setdefault(row.regime, []).append(row)
    for row in regime_rows:
        profile = profile_by_regime.get(row.regime, [])
        if not profile:
            continue
        supported_count = sum(
            1 for candidate in profile if candidate.in_support_interval
        )
        if supported_count >= math.ceil(len(profile) / 2.0):
            warnings.append(
                BrownianRegimeIdentifiabilityWarning(
                    regime=row.regime,
                    kind="flat_profile",
                    message=(
                        f"regime '{row.regime}' has a broad conditional likelihood "
                        "profile, so its Brownian rate is weakly identified"
                    ),
                )
            )
        if math.isclose(
            row.lower_95,
            min(candidate.sigma_squared for candidate in profile),
            rel_tol=0.0,
            abs_tol=1e-12,
        ) or math.isclose(
            row.upper_95,
            max(candidate.sigma_squared for candidate in profile),
            rel_tol=0.0,
            abs_tol=1e-12,
        ):
            warnings.append(
                BrownianRegimeIdentifiabilityWarning(
                    regime=row.regime,
                    kind="boundary_rate",
                    message=(
                        f"regime '{row.regime}' reaches a profile search boundary, "
                        "so its Brownian rate is weakly identified"
                    ),
                )
            )
        if row.contributing_branch_length <= 0.5:
            warnings.append(
                BrownianRegimeIdentifiabilityWarning(
                    regime=row.regime,
                    kind="low_branch_exposure",
                    message=(
                        f"regime '{row.regime}' contributes very little analyzed "
                        "branch length, so its Brownian rate is weakly identified"
                    ),
                )
            )
    if better_model != "brownian-regimes":
        warnings.append(
            BrownianRegimeIdentifiabilityWarning(
                regime="*",
                kind="comparison_not_preferred",
                message=(
                    "single-rate Brownian remains preferred over the multi-rate "
                    "regime model by AICc on the analyzed taxon set"
                ),
            )
        )
    return warnings


def _build_multirate_residual_diagnostics(
    dataset: ComparativeDataset,
    covariance: list[list[float]],
    residuals: list[float],
) -> ComparativeResidualSummary:
    residual_mean = sum(residuals) / len(residuals)
    residual_variance = sum((value - residual_mean) ** 2 for value in residuals) / max(
        1,
        len(residuals) - 1,
    )
    residual_skewness = 0.0
    if residual_variance > 0.0:
        residual_sd = math.sqrt(residual_variance)
        residual_skewness = sum(
            ((value - residual_mean) / residual_sd) ** 3 for value in residuals
        ) / len(residuals)
    inverse_covariance = invert_matrix(covariance)
    hat = _hat_matrix(len(dataset.taxa), covariance, inverse_covariance)
    outliers: list[ComparativeResidualOutlier] = []
    max_abs_standardized = 0.0
    for index, taxon in enumerate(dataset.taxa):
        leverage = min(max(hat[index][index], 0.0), 0.999999)
        denominator = math.sqrt(max(covariance[index][index] * (1.0 - leverage), 1e-12))
        standardized = residuals[index] / denominator
        max_abs_standardized = max(max_abs_standardized, abs(standardized))
        if abs(standardized) >= 2.0:
            outliers.append(
                ComparativeResidualOutlier(
                    taxon=taxon,
                    residual=residuals[index],
                    standardized_residual=standardized,
                    note="absolute standardized residual exceeds 2.0",
                )
            )
    residual_lambda = _estimate_lambda_for_values(dataset, residuals)
    warnings: list[str] = []
    if abs(residual_skewness) > 1.0:
        warnings.append("residual distribution is noticeably skewed")
    if residual_lambda > 0.5:
        warnings.append("residuals retain moderate phylogenetic structure")
    if outliers:
        warnings.append("one or more taxa have unusually large residuals")
    return ComparativeResidualSummary(
        residual_mean=residual_mean,
        residual_variance=residual_variance,
        residual_skewness=residual_skewness,
        max_abs_standardized_residual=max_abs_standardized,
        phylogenetic_residual_lambda=residual_lambda,
        outlier_taxa=outliers,
        warnings=warnings,
    )


def _hat_matrix(
    sample_size: int,
    covariance: list[list[float]],
    inverse_covariance: list[list[float]],
) -> list[list[float]]:
    design = [[1.0] for _ in range(sample_size)]
    xt_vinv = matrix_multiply(transpose(design), inverse_covariance)
    xt_vinv_x_inverse = invert_matrix(matrix_multiply(xt_vinv, design))
    return matrix_multiply(design, matrix_multiply(xt_vinv_x_inverse, xt_vinv))


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


def _chi_square_survival(statistic: float, degrees_of_freedom: int) -> float:
    if statistic <= 0.0:
        return 1.0
    if degrees_of_freedom <= 0:
        return 1.0
    if degrees_of_freedom == 1:
        return math.erfc(math.sqrt(statistic / 2.0))
    z_score = (
        ((statistic / degrees_of_freedom) ** (1.0 / 3.0))
        - (1.0 - (2.0 / (9.0 * degrees_of_freedom)))
    ) / math.sqrt(2.0 / (9.0 * degrees_of_freedom))
    return 0.5 * math.erfc(z_score / math.sqrt(2.0))

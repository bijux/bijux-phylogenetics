from __future__ import annotations

import math

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    matrix_multiply,
    transpose,
)
from bijux_phylogenetics.comparative.common import ComparativeDataset
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeResidualOutlier,
    ComparativeResidualSummary,
    _comparison_row,
    _estimate_lambda_for_values,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)

from .contracts import (
    BrownianRegimeBranchRow,
    BrownianRegimeIdentifiabilityWarning,
    BrownianRegimeProfileRow,
    BrownianRegimeRateRow,
)
from .covariance_fitting import fit_multirate_covariance, logspace

_PROFILE_CONFIDENCE_DELTA = 1.920729410347062


def build_comparison_rows(
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


def build_profile_rows(
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
        candidates = sorted({*logspace(lower, upper, 81), fitted_rates[regime]})
        profile = []
        best_log_likelihood = -math.inf
        for candidate in candidates:
            trial_rates = dict(fitted_rates)
            trial_rates[regime] = candidate
            fit = fit_multirate_covariance(dataset, regime_matrices, trial_rates)
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


def build_regime_rows(
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


def build_identifiability_warnings(
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


def build_multirate_residual_diagnostics(
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


def chi_square_survival(statistic: float, degrees_of_freedom: int) -> float:
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


def _hat_matrix(
    sample_size: int,
    covariance: list[list[float]],
    inverse_covariance: list[list[float]],
) -> list[list[float]]:
    design = [[1.0] for _ in range(sample_size)]
    xt_vinv = matrix_multiply(transpose(design), inverse_covariance)
    xt_vinv_x_inverse = invert_matrix(matrix_multiply(xt_vinv, design))
    return matrix_multiply(design, matrix_multiply(xt_vinv_x_inverse, xt_vinv))

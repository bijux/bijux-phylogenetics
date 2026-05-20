from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    matrix_multiply,
    quadratic_form,
    stable_covariance,
    transpose,
)
from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    build_brownian_covariance_matrix,
    build_ou_covariance_matrix,
    load_comparative_dataset,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonReport,
    ComparativeModelComparisonRow,
)

_Z_95 = 1.959963984540054


@dataclass(slots=True)
class ComparativeParameterInterval:
    """Approximate 95% interval for one fitted comparative-model parameter."""

    name: str
    estimate: float
    lower_95: float
    upper_95: float
    method: str


@dataclass(slots=True)
class ComparativeResidualOutlier:
    """One taxon with an unusually large model residual."""

    taxon: str
    residual: float
    standardized_residual: float
    note: str


@dataclass(slots=True)
class ComparativeResidualSummary:
    """Residual diagnostics for one fitted comparative model."""

    residual_mean: float
    residual_variance: float
    residual_skewness: float
    max_abs_standardized_residual: float
    phylogenetic_residual_lambda: float
    outlier_taxa: list[ComparativeResidualOutlier]
    warnings: list[str]


@dataclass(slots=True)
class BrownianMotionFitReport:
    """Standalone Brownian-motion continuous-trait fit."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    taxon_count: int
    taxa: list[str]
    root_state: float
    rate: float
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]
    assumptions: list[str]
    confidence_intervals: list[ComparativeParameterInterval]
    residual_diagnostics: ComparativeResidualSummary


@dataclass(slots=True)
class OUIdentifiabilityWarning:
    """Warning that an OU fit may not be statistically identifiable."""

    kind: str
    message: str


@dataclass(slots=True)
class OUTraitModelReport:
    """Standalone Ornstein-Uhlenbeck continuous-trait fit."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    taxon_count: int
    taxa: list[str]
    alpha: float
    theta: float
    sigma_squared: float
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]
    convergence_status: str
    assumptions: list[str]
    confidence_intervals: list[ComparativeParameterInterval]
    identifiability_warnings: list[OUIdentifiabilityWarning]
    residual_diagnostics: ComparativeResidualSummary


def fit_brownian_motion_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> BrownianMotionFitReport:
    """Fit a standalone Brownian-motion model for one continuous trait."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
    )
    fit = _fit_intercept_only_model(dataset, covariance)
    intervals = _brownian_parameter_intervals(fit.theta, fit.sigma_squared, covariance)
    residual_diagnostics = _build_residual_diagnostics(
        dataset, covariance, fit.residuals, fit.sigma_squared
    )
    return BrownianMotionFitReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        taxon_count=len(dataset.taxa),
        taxa=dataset.taxa,
        root_state=fit.theta,
        rate=fit.sigma_squared,
        log_likelihood=fit.log_likelihood,
        fitted_values=fit.fitted_values,
        residuals=fit.residuals,
        assumptions=[
            "Brownian motion assumes trait variance accumulates proportionally with shared branch length",
            "Brownian motion assumes no directional optimum and a constant diffusion rate across the tree",
        ],
        confidence_intervals=intervals,
        residual_diagnostics=residual_diagnostics,
    )


def fit_ornstein_uhlenbeck_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> OUTraitModelReport:
    """Fit a stationary-root Ornstein-Uhlenbeck model for one continuous trait."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    alpha_grid = _alpha_grid(dataset)
    best_alpha = alpha_grid[0]
    best_fit: _InterceptOnlyFit | None = None
    profile: list[tuple[float, float]] = []
    for alpha in alpha_grid:
        covariance = stable_covariance(
            build_ou_covariance_matrix(dataset.tree, dataset.taxa, alpha=alpha)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        profile.append((alpha, fit.log_likelihood))
        if best_fit is None or fit.log_likelihood > best_fit.log_likelihood:
            best_alpha = alpha
            best_fit = fit
    if best_fit is None:
        raise ValueError(
            "failed to fit OU intercept-only model for the provided alpha grid"
        )
    best_covariance = stable_covariance(
        build_ou_covariance_matrix(dataset.tree, dataset.taxa, alpha=best_alpha)
    )
    intervals = _ou_parameter_intervals(
        best_alpha, best_fit.theta, best_fit.sigma_squared, best_covariance, profile
    )
    identifiability_warnings = _ou_identifiability_warnings(
        dataset,
        best_alpha,
        profile,
    )
    residual_diagnostics = _build_residual_diagnostics(
        dataset,
        best_covariance,
        best_fit.residuals,
        best_fit.sigma_squared,
    )
    return OUTraitModelReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=dataset.taxon_column,
        trait=trait,
        taxon_count=len(dataset.taxa),
        taxa=dataset.taxa,
        alpha=best_alpha,
        theta=best_fit.theta,
        sigma_squared=best_fit.sigma_squared,
        log_likelihood=best_fit.log_likelihood,
        fitted_values=best_fit.fitted_values,
        residuals=best_fit.residuals,
        convergence_status="grid-search-converged",
        assumptions=[
            "OU assumes trait evolution is pulled toward a stationary optimum theta",
            "OU fit here uses a rooted stationary-process covariance and grid-search over alpha",
        ],
        confidence_intervals=intervals,
        identifiability_warnings=identifiability_warnings,
        residual_diagnostics=residual_diagnostics,
    )


def compare_brownian_and_ou_models(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
) -> ComparativeModelComparisonReport:
    """Compare standalone Brownian-motion and OU models by likelihood and information criteria."""
    brownian = fit_brownian_motion_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    ou = fit_ornstein_uhlenbeck_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    n = brownian.taxon_count
    rows = [
        _comparison_row("brownian", 2, brownian.log_likelihood, n),
        _comparison_row("ou", 3, ou.log_likelihood, n),
    ]
    best_aicc = min(row.aicc for row in rows)
    for row in rows:
        row.selected = math.isclose(row.aicc, best_aicc, abs_tol=1e-12)
    better_model = next(row.model for row in rows if row.selected)
    return ComparativeModelComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=n,
        rows=rows,
        better_model=better_model,
        warnings=[],
    )


@dataclass(slots=True)
class _InterceptOnlyFit:
    theta: float
    sigma_squared: float
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]


def _fit_intercept_only_model(
    dataset: ComparativeDataset,
    covariance: list[list[float]],
) -> _InterceptOnlyFit:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(dataset.trait_values)
    denom = quadratic_form(ones, inverse_covariance)
    theta = (
        sum(
            ones[row_index]
            * sum(
                inverse_covariance[row_index][column_index]
                * dataset.trait_values[column_index]
                for column_index in range(len(dataset.trait_values))
            )
            for row_index in range(len(dataset.trait_values))
        )
        / denom
    )
    fitted_values = [theta] * len(dataset.trait_values)
    residuals = [value - theta for value in dataset.trait_values]
    sigma_squared = quadratic_form(residuals, inverse_covariance) / len(
        dataset.trait_values
    )
    log_likelihood = -0.5 * (
        len(dataset.trait_values) * math.log(2.0 * math.pi * sigma_squared)
        + log_determinant(covariance)
        + len(dataset.trait_values)
    )
    return _InterceptOnlyFit(
        theta=theta,
        sigma_squared=sigma_squared,
        log_likelihood=log_likelihood,
        fitted_values=fitted_values,
        residuals=residuals,
    )


def _brownian_parameter_intervals(
    root_state: float,
    sigma_squared: float,
    covariance: list[list[float]],
) -> list[ComparativeParameterInterval]:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(covariance)
    denom = quadratic_form(ones, inverse_covariance)
    root_se = math.sqrt(sigma_squared / denom)
    rate_lower, rate_upper = _variance_interval(
        sigma_squared,
        degrees_of_freedom=max(1, len(covariance) - 1),
    )
    return [
        ComparativeParameterInterval(
            name="root_state",
            estimate=root_state,
            lower_95=root_state - (_Z_95 * root_se),
            upper_95=root_state + (_Z_95 * root_se),
            method="wald",
        ),
        ComparativeParameterInterval(
            name="rate",
            estimate=sigma_squared,
            lower_95=rate_lower,
            upper_95=rate_upper,
            method="chi-square-approximation",
        ),
    ]


def _ou_parameter_intervals(
    alpha: float,
    theta: float,
    sigma_squared: float,
    covariance: list[list[float]],
    profile: list[tuple[float, float]],
) -> list[ComparativeParameterInterval]:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(covariance)
    denom = quadratic_form(ones, inverse_covariance)
    theta_se = math.sqrt(sigma_squared / denom)
    sigma_lower, sigma_upper = _variance_interval(
        sigma_squared,
        degrees_of_freedom=max(1, len(covariance) - 1),
    )
    best_log_likelihood = max(log_likelihood for _, log_likelihood in profile)
    supported = [
        candidate
        for candidate, log_likelihood in profile
        if log_likelihood >= best_log_likelihood - 1.92
    ]
    return [
        ComparativeParameterInterval(
            name="alpha",
            estimate=alpha,
            lower_95=min(supported),
            upper_95=max(supported),
            method="profile-likelihood-grid",
        ),
        ComparativeParameterInterval(
            name="theta",
            estimate=theta,
            lower_95=theta - (_Z_95 * theta_se),
            upper_95=theta + (_Z_95 * theta_se),
            method="wald",
        ),
        ComparativeParameterInterval(
            name="sigma_squared",
            estimate=sigma_squared,
            lower_95=sigma_lower,
            upper_95=sigma_upper,
            method="chi-square-approximation",
        ),
    ]


def _variance_interval(
    sigma_squared: float, *, degrees_of_freedom: int
) -> tuple[float, float]:
    chi_upper = _chi_square_quantile(0.975, degrees_of_freedom)
    chi_lower = _chi_square_quantile(0.025, degrees_of_freedom)
    lower = (degrees_of_freedom * sigma_squared) / chi_upper
    upper = (degrees_of_freedom * sigma_squared) / chi_lower
    return lower, upper


def _chi_square_quantile(probability: float, degrees_of_freedom: int) -> float:
    if probability not in {0.025, 0.975}:
        raise ValueError(
            "supported chi-square approximation probabilities are 0.025 and 0.975"
        )
    z = -_Z_95 if probability < 0.5 else _Z_95
    factor = (
        1.0
        - (2.0 / (9.0 * degrees_of_freedom))
        + z * math.sqrt(2.0 / (9.0 * degrees_of_freedom))
    )
    return degrees_of_freedom * (factor**3)


def _alpha_grid(dataset: ComparativeDataset) -> list[float]:
    max_depth = max(
        max(row) for row in build_brownian_covariance_matrix(dataset.tree, dataset.taxa)
    )
    scale = 1.0 / max(max_depth, 1e-6)
    coarse = [
        round(scale * multiplier, 6)
        for multiplier in (0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0)
    ]
    return sorted({value for value in coarse if value > 0.0})


def _ou_identifiability_warnings(
    dataset: ComparativeDataset,
    alpha: float,
    profile: list[tuple[float, float]],
) -> list[OUIdentifiabilityWarning]:
    warnings: list[OUIdentifiabilityWarning] = []
    if len(dataset.taxa) < 5:
        warnings.append(
            OUIdentifiabilityWarning(
                kind="small_sample_size",
                message="OU alpha is hard to identify with fewer than five taxa",
            )
        )
    ordered_alphas = sorted(candidate for candidate, _ in profile)
    if math.isclose(
        alpha, ordered_alphas[0], rel_tol=0.0, abs_tol=1e-9
    ) or math.isclose(alpha, ordered_alphas[-1], rel_tol=0.0, abs_tol=1e-9):
        warnings.append(
            OUIdentifiabilityWarning(
                kind="boundary_alpha",
                message="best-supported OU alpha falls on the search boundary and may not be well identified",
            )
        )
    best_log_likelihood = max(log_likelihood for _, log_likelihood in profile)
    second_best = sorted(
        (log_likelihood for _, log_likelihood in profile), reverse=True
    )[1]
    if best_log_likelihood - second_best < 0.5:
        warnings.append(
            OUIdentifiabilityWarning(
                kind="flat_likelihood",
                message="OU likelihood surface is shallow across alpha values, so model choice may be unstable",
            )
        )
    if alpha < ordered_alphas[len(ordered_alphas) // 3]:
        warnings.append(
            OUIdentifiabilityWarning(
                kind="weak_pull_to_optimum",
                message="best-supported OU alpha is weak and may be difficult to distinguish from Brownian motion",
            )
        )
    return warnings


def _build_residual_diagnostics(
    dataset: ComparativeDataset,
    covariance: list[list[float]],
    residuals: list[float],
    sigma_squared: float,
) -> ComparativeResidualSummary:
    residual_mean = sum(residuals) / len(residuals)
    residual_variance = sum((value - residual_mean) ** 2 for value in residuals) / max(
        1, len(residuals) - 1
    )
    residual_skewness = 0.0
    if residual_variance > 0.0:
        residual_sd = math.sqrt(residual_variance)
        residual_skewness = sum(
            ((value - residual_mean) / residual_sd) ** 3 for value in residuals
        ) / len(residuals)
    inverse_covariance = invert_matrix(covariance)
    residual_lambda = _estimate_lambda_for_values(dataset, residuals)
    standardized_rows = _standardized_residual_rows(
        dataset.taxa,
        covariance,
        inverse_covariance,
        residuals,
        sigma_squared,
    )
    outliers = [
        ComparativeResidualOutlier(
            taxon=taxon,
            residual=residual,
            standardized_residual=standardized,
            note="absolute standardized residual exceeds 2.0",
        )
        for taxon, residual, standardized in standardized_rows
        if abs(standardized) >= 2.0
    ]
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
        max_abs_standardized_residual=max(abs(row[2]) for row in standardized_rows),
        phylogenetic_residual_lambda=residual_lambda,
        outlier_taxa=outliers,
        warnings=warnings,
    )


def _standardized_residual_rows(
    taxa: list[str],
    covariance: list[list[float]],
    inverse_covariance: list[list[float]],
    residuals: list[float],
    sigma_squared: float,
) -> list[tuple[str, float, float]]:
    hat = _hat_matrix(len(taxa), covariance, inverse_covariance)
    rows: list[tuple[str, float, float]] = []
    for index, taxon in enumerate(taxa):
        leverage = min(max(hat[index][index], 0.0), 0.999999)
        denominator = math.sqrt(max(sigma_squared * (1.0 - leverage), 1e-12))
        rows.append((taxon, residuals[index], residuals[index] / denominator))
    return rows


def _hat_matrix(
    sample_size: int,
    covariance: list[list[float]],
    inverse_covariance: list[list[float]],
) -> list[list[float]]:
    design = [[1.0] for _ in range(sample_size)]
    xt_vinv = matrix_multiply(transpose(design), inverse_covariance)
    xt_vinv_x_inverse = invert_matrix(matrix_multiply(xt_vinv, design))
    return matrix_multiply(design, matrix_multiply(xt_vinv_x_inverse, xt_vinv))


def _estimate_lambda_for_values(
    dataset: ComparativeDataset, values: list[float]
) -> float:
    candidates = [index / 20.0 for index in range(21)]
    best_lambda = 0.0
    best_score = -math.inf
    for lambda_value in candidates:
        covariance = [
            [
                value if row_index == column_index else value * lambda_value
                for column_index, value in enumerate(row)
            ]
            for row_index, row in enumerate(dataset.covariance_matrix)
        ]
        covariance = stable_covariance(covariance)
        inverse_covariance = invert_matrix(covariance)
        ones = [1.0] * len(values)
        denom = quadratic_form(ones, inverse_covariance)
        mean_value = (
            sum(
                ones[row_index]
                * sum(
                    inverse_covariance[row_index][column_index] * values[column_index]
                    for column_index in range(len(values))
                )
                for row_index in range(len(values))
            )
            / denom
        )
        residuals = [value - mean_value for value in values]
        sigma_squared = quadratic_form(residuals, inverse_covariance) / len(values)
        score = -0.5 * (
            len(values) * math.log(2.0 * math.pi * sigma_squared)
            + log_determinant(covariance)
            + len(values)
        )
        if score > best_score:
            best_score = score
            best_lambda = lambda_value
    return best_lambda


def _comparison_row(
    model: str,
    parameter_count: int,
    log_likelihood: float,
    sample_size: int,
    *,
    likelihood_constant_policy: str | None = None,
) -> ComparativeModelComparisonRow:
    aic = (2.0 * parameter_count) - (2.0 * log_likelihood)
    if sample_size <= parameter_count + 1:
        aicc = math.inf
    else:
        aicc = aic + (
            (2.0 * parameter_count * (parameter_count + 1))
            / (sample_size - parameter_count - 1)
        )
    return ComparativeModelComparisonRow(
        model=model,
        parameter_count=parameter_count,
        log_likelihood=log_likelihood,
        aic=aic,
        aicc=aicc,
        delta_aic=None,
        delta_aicc=None,
        rank=None,
        comparable=True,
        comparability_note=None,
        selected=False,
        likelihood_constant_policy=likelihood_constant_policy,
    )

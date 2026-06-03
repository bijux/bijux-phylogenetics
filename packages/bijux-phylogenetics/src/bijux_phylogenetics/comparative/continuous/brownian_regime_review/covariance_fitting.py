from __future__ import annotations

from dataclasses import dataclass
import math

from bijux_phylogenetics.comparative._math import (
    invert_matrix,
    log_determinant,
    quadratic_form,
    stable_covariance,
)
from bijux_phylogenetics.comparative.common import ComparativeDataset
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeParameterInterval,
)

_Z_95 = 1.959963984540054


@dataclass(slots=True)
class MultirateFit:
    regime_rates: dict[str, float]
    root_state: float
    root_state_interval: ComparativeParameterInterval
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]
    covariance: list[list[float]]


@dataclass(slots=True)
class CovarianceFit:
    root_state: float
    log_likelihood: float
    fitted_values: list[float]
    residuals: list[float]
    covariance: list[list[float]]


def fit_multirate_brownian_model(
    dataset: ComparativeDataset,
    *,
    regime_matrices: dict[str, list[list[float]]],
    baseline_rate: float,
) -> MultirateFit:
    lower = max(baseline_rate * 0.02, 1e-6)
    upper = max(baseline_rate * 50.0, lower * 10.0)
    regimes = sorted(regime_matrices)
    bounds = dict.fromkeys(regimes, (lower, upper))
    regime_rates = dict.fromkeys(regimes, baseline_rate)
    for _ in range(6):
        for regime in regimes:
            search = logspace(bounds[regime][0], bounds[regime][1], 41)
            best_index = 0
            best_log_likelihood = -math.inf
            best_rate = search[0]
            for index, candidate in enumerate(search):
                trial_rates = dict(regime_rates)
                trial_rates[regime] = candidate
                fit = fit_multirate_covariance(
                    dataset,
                    regime_matrices,
                    trial_rates,
                )
                if fit.log_likelihood > best_log_likelihood:
                    best_log_likelihood = fit.log_likelihood
                    best_index = index
                    best_rate = candidate
            regime_rates[regime] = best_rate
            bounds[regime] = (
                search[max(0, best_index - 3)],
                search[min(len(search) - 1, best_index + 3)],
            )
    fit = fit_multirate_covariance(dataset, regime_matrices, regime_rates)
    return MultirateFit(
        regime_rates=regime_rates,
        root_state=fit.root_state,
        root_state_interval=root_state_interval(fit.root_state, fit.covariance),
        log_likelihood=fit.log_likelihood,
        fitted_values=fit.fitted_values,
        residuals=fit.residuals,
        covariance=fit.covariance,
    )


def fit_multirate_covariance(
    dataset: ComparativeDataset,
    regime_matrices: dict[str, list[list[float]]],
    regime_rates: dict[str, float],
) -> CovarianceFit:
    covariance = stable_covariance(
        combine_regime_covariance(regime_matrices, regime_rates)
    )
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(dataset.trait_values)
    denom = quadratic_form(ones, inverse_covariance)
    root_state = (
        sum(
            sum(
                inverse_covariance[row_index][column_index]
                * dataset.trait_values[column_index]
                for column_index in range(len(dataset.trait_values))
            )
            for row_index in range(len(dataset.trait_values))
        )
        / denom
    )
    fitted_values = [root_state] * len(dataset.trait_values)
    residuals = [value - root_state for value in dataset.trait_values]
    log_likelihood = -0.5 * (
        len(dataset.trait_values) * math.log(2.0 * math.pi)
        + log_determinant(covariance)
        + quadratic_form(residuals, inverse_covariance)
    )
    return CovarianceFit(
        root_state=root_state,
        log_likelihood=log_likelihood,
        fitted_values=fitted_values,
        residuals=residuals,
        covariance=covariance,
    )


def combine_regime_covariance(
    regime_matrices: dict[str, list[list[float]]],
    regime_rates: dict[str, float],
) -> list[list[float]]:
    size = len(next(iter(regime_matrices.values())))
    combined = [[0.0] * size for _ in range(size)]
    for regime, matrix in regime_matrices.items():
        rate = regime_rates[regime]
        for row_index in range(size):
            for column_index in range(size):
                combined[row_index][column_index] += (
                    rate * matrix[row_index][column_index]
                )
    return combined


def root_state_interval(
    root_state: float,
    covariance: list[list[float]],
) -> ComparativeParameterInterval:
    inverse_covariance = invert_matrix(covariance)
    ones = [1.0] * len(covariance)
    denom = quadratic_form(ones, inverse_covariance)
    root_se = math.sqrt(max(1.0 / denom, 1e-12))
    return ComparativeParameterInterval(
        name="root_state",
        estimate=root_state,
        lower_95=root_state - (_Z_95 * root_se),
        upper_95=root_state + (_Z_95 * root_se),
        method="wald",
    )


def logspace(lower: float, upper: float, count: int) -> list[float]:
    if count < 2:
        return [lower]
    low_log = math.log(lower)
    high_log = math.log(upper)
    return [
        math.exp(low_log + ((high_log - low_log) * (index / float(count - 1))))
        for index in range(count)
    ]

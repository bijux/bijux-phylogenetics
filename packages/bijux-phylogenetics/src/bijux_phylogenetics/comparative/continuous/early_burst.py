from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    ComparativeReadinessReport,
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    stable_covariance,
    summarize_numeric_trait_readiness,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    ComparativeParameterInterval,
    ComparativeResidualSummary,
    _comparison_row,
    _fit_intercept_only_model,
)
from bijux_phylogenetics.comparative.evolutionary_modes import (
    ContinuousEvolutionaryModeComparisonReport,
    LikelihoodRatioTestResult,
    compare_continuous_evolutionary_modes,
    fit_continuous_evolutionary_mode,
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
)
from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

_PROFILE_CONFIDENCE_DELTA = 1.920729410347062


@dataclass(slots=True)
class EarlyBurstTraitEvolutionExclusion:
    """One taxon excluded before early-burst trait-evolution fitting."""

    taxon: str
    reason: str


@dataclass(slots=True)
class EarlyBurstIdentifiabilityWarning:
    """Warning that the early-burst rate-change surface is weakly identified."""

    kind: str
    message: str


@dataclass(slots=True)
class EarlyBurstRateChangeProfileRow:
    """One fixed rate-change fit on the early-burst likelihood profile."""

    rate_change: float
    log_likelihood: float
    aic: float
    aicc: float
    delta_log_likelihood: float
    in_support_interval: bool
    selected: bool


@dataclass(slots=True)
class EarlyBurstTraitEvolutionSummaryReport:
    """Reviewer-facing early-burst trait-evolution fit with comparison context."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    tree_taxon_count: int
    analyzed_taxa: list[str]
    analyzed_taxon_count: int
    excluded_taxa: list[EarlyBurstTraitEvolutionExclusion]
    rate_change_bounds: tuple[float, float]
    rate_change: float
    root_state: float
    sigma_squared: float
    log_likelihood: float
    aic: float
    aicc: float
    confidence_intervals: list[ComparativeParameterInterval]
    residual_diagnostics: ComparativeResidualSummary
    comparison_rows: list[ComparativeModelComparisonRow]
    likelihood_ratio_tests: list[LikelihoodRatioTestResult]
    better_model: str
    profile_rows: list[EarlyBurstRateChangeProfileRow]
    identifiability_warnings: list[EarlyBurstIdentifiabilityWarning]
    assumptions: list[str]
    warnings: list[str]
    readiness: ComparativeReadinessReport


def summarize_early_burst_trait_evolution(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    rate_change_bounds: tuple[float, float] = (0.0, 50.0),
) -> EarlyBurstTraitEvolutionSummaryReport:
    """Summarize one early-burst continuous-trait fit with explicit comparison context."""
    readiness = summarize_numeric_trait_readiness(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    fit = fit_continuous_evolutionary_mode(
        tree_path,
        traits_path,
        trait=trait,
        mode="early-burst",
        taxon_column=taxon_column,
        early_burst_bounds=rate_change_bounds,
    )
    comparison = compare_continuous_evolutionary_modes(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        early_burst_bounds=rate_change_bounds,
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
    profile_rows = _build_rate_change_profile(dataset, rate_change_bounds)
    rate_change_interval = _build_rate_change_interval(
        profile_rows,
        fit.parameter_value or 0.0,
    )
    identifiability_warnings = _detect_identifiability_warnings(
        rate_change=fit.parameter_value or 0.0,
        rate_change_bounds=rate_change_bounds,
        profile_rows=profile_rows,
        comparison=comparison,
    )
    warnings = list(
        dict.fromkeys(
            [
                *readiness.warnings,
                *fit.residual_diagnostics.warnings,
                *[warning.message for warning in identifiability_warnings],
            ]
        )
    )
    assumptions = [
        *fit.assumptions,
        "The early-burst fit is reviewed against Brownian and OU fits on the same pruned taxon set.",
        "Weak identifiability is flagged from the bounded rate-change profile and model-selection context.",
    ]
    fit_intervals = [
        interval
        for interval in fit.confidence_intervals
        if isinstance(interval, ComparativeParameterInterval)
    ]
    return EarlyBurstTraitEvolutionSummaryReport(
        tree_path=tree_path,
        traits_path=traits_path,
        taxon_column=fit.taxon_column,
        trait=trait,
        tree_taxon_count=readiness.tree_taxa,
        analyzed_taxa=list(fit.taxa),
        analyzed_taxon_count=fit.taxon_count,
        excluded_taxa=_build_excluded_taxa(readiness),
        rate_change_bounds=rate_change_bounds,
        rate_change=fit.parameter_value or 0.0,
        root_state=fit.root_state,
        sigma_squared=fit.rate,
        log_likelihood=fit.log_likelihood,
        aic=fit.aic,
        aicc=fit.aicc,
        confidence_intervals=[rate_change_interval, *fit_intervals],
        residual_diagnostics=fit.residual_diagnostics,
        comparison_rows=list(comparison.rows),
        likelihood_ratio_tests=list(comparison.likelihood_ratio_tests),
        better_model=comparison.better_model,
        profile_rows=profile_rows,
        identifiability_warnings=identifiability_warnings,
        assumptions=assumptions,
        warnings=warnings,
        readiness=readiness,
    )


def write_early_burst_trait_evolution_summary_table(
    path: Path,
    report: EarlyBurstTraitEvolutionSummaryReport,
) -> Path:
    """Write one summary ledger for an early-burst trait-evolution fit."""
    interval_by_name = {
        interval.name: interval for interval in report.confidence_intervals
    }
    rate_change_interval = interval_by_name.get("rate_change")
    sigma_interval = interval_by_name.get("rate")
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "taxon_column",
            "tree_taxon_count",
            "analyzed_taxon_count",
            "excluded_taxon_count",
            "rate_change",
            "rate_change_lower_95",
            "rate_change_upper_95",
            "root_state",
            "sigma_squared",
            "sigma_squared_lower_95",
            "sigma_squared_upper_95",
            "log_likelihood",
            "aic",
            "aicc",
            "better_model",
            "identifiability_warning_count",
            "residual_variance",
            "max_abs_standardized_residual",
            "phylogenetic_residual_lambda",
        ],
        rows=[
            {
                "trait": report.trait,
                "taxon_column": report.taxon_column,
                "tree_taxon_count": report.tree_taxon_count,
                "analyzed_taxon_count": report.analyzed_taxon_count,
                "excluded_taxon_count": len(report.excluded_taxa),
                "rate_change": format(report.rate_change, ".15g"),
                "rate_change_lower_95": _format_interval_bound(
                    rate_change_interval, "lower_95"
                ),
                "rate_change_upper_95": _format_interval_bound(
                    rate_change_interval, "upper_95"
                ),
                "root_state": format(report.root_state, ".15g"),
                "sigma_squared": format(report.sigma_squared, ".15g"),
                "sigma_squared_lower_95": _format_interval_bound(
                    sigma_interval, "lower_95"
                ),
                "sigma_squared_upper_95": _format_interval_bound(
                    sigma_interval, "upper_95"
                ),
                "log_likelihood": format(report.log_likelihood, ".15g"),
                "aic": format(report.aic, ".15g"),
                "aicc": format(report.aicc, ".15g"),
                "better_model": report.better_model,
                "identifiability_warning_count": len(report.identifiability_warnings),
                "residual_variance": format(
                    report.residual_diagnostics.residual_variance, ".15g"
                ),
                "max_abs_standardized_residual": format(
                    report.residual_diagnostics.max_abs_standardized_residual, ".15g"
                ),
                "phylogenetic_residual_lambda": format(
                    report.residual_diagnostics.phylogenetic_residual_lambda, ".15g"
                ),
            }
        ],
    )


def write_early_burst_trait_evolution_exclusion_table(
    path: Path,
    report: EarlyBurstTraitEvolutionSummaryReport,
) -> Path:
    """Write one excluded-taxon ledger for an early-burst trait-evolution fit."""
    return write_taxon_rows(
        path,
        columns=["taxon", "reason"],
        rows=[
            {"taxon": row.taxon, "reason": row.reason} for row in report.excluded_taxa
        ],
    )


def write_early_burst_trait_evolution_comparison_table(
    path: Path,
    report: EarlyBurstTraitEvolutionSummaryReport,
) -> Path:
    """Write one combined comparison ledger for early-burst, BM, and OU fits."""
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
            "left_mode": "",
            "right_mode": "",
            "statistic": "",
            "degrees_of_freedom": "",
            "p_value": "",
        }
        for row in report.comparison_rows
    ]
    rows.extend(
        {
            "row_kind": "likelihood_ratio_test",
            "model": "",
            "comparison_id": row.comparison_id,
            "parameter_count": "",
            "log_likelihood": "",
            "aic": "",
            "aicc": "",
            "delta_aicc": "",
            "selected": "",
            "left_mode": row.left_mode,
            "right_mode": row.right_mode,
            "statistic": format(row.statistic, ".15g"),
            "degrees_of_freedom": row.degrees_of_freedom,
            "p_value": format(row.p_value, ".15g"),
        }
        for row in report.likelihood_ratio_tests
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
            "left_mode",
            "right_mode",
            "statistic",
            "degrees_of_freedom",
            "p_value",
        ],
        rows=rows,
    )


def write_early_burst_rate_change_profile_table(
    path: Path,
    report: EarlyBurstTraitEvolutionSummaryReport,
) -> Path:
    """Write one fixed-parameter likelihood profile for the early-burst fit."""
    return write_taxon_rows(
        path,
        columns=[
            "trait",
            "rate_change",
            "log_likelihood",
            "aic",
            "aicc",
            "delta_log_likelihood",
            "in_support_interval",
            "selected",
        ],
        rows=[
            {
                "trait": report.trait,
                "rate_change": format(row.rate_change, ".15g"),
                "log_likelihood": format(row.log_likelihood, ".15g"),
                "aic": format(row.aic, ".15g"),
                "aicc": format(row.aicc, ".15g"),
                "delta_log_likelihood": format(row.delta_log_likelihood, ".15g"),
                "in_support_interval": str(row.in_support_interval).lower(),
                "selected": str(row.selected).lower(),
            }
            for row in report.profile_rows
        ],
    )


def _build_excluded_taxa(
    readiness: ComparativeReadinessReport,
) -> list[EarlyBurstTraitEvolutionExclusion]:
    rows: list[EarlyBurstTraitEvolutionExclusion] = []
    rows.extend(
        EarlyBurstTraitEvolutionExclusion(
            taxon=taxon,
            reason="missing_from_trait_table",
        )
        for taxon in readiness.missing_from_traits
    )
    rows.extend(
        EarlyBurstTraitEvolutionExclusion(
            taxon=taxon,
            reason="missing_trait_value",
        )
        for taxon in readiness.pruned_missing_value_taxa
    )
    rows.extend(
        EarlyBurstTraitEvolutionExclusion(
            taxon=taxon,
            reason="non_numeric_trait_value",
        )
        for taxon in readiness.pruned_non_numeric_taxa
    )
    rows.extend(
        EarlyBurstTraitEvolutionExclusion(
            taxon=taxon,
            reason="absent_from_tree",
        )
        for taxon in readiness.extra_trait_taxa
    )
    return rows


def _build_rate_change_profile(
    dataset: ComparativeDataset,
    rate_change_bounds: tuple[float, float],
) -> list[EarlyBurstRateChangeProfileRow]:
    lower, upper = rate_change_bounds
    grid = {lower + ((upper - lower) * (index / 160.0)) for index in range(161)}
    if math.isclose(lower, upper, rel_tol=0.0, abs_tol=1e-12):
        grid = {lower}
    fits = [
        _fit_early_burst_candidate(
            dataset,
            rate_change=rate_change,
        )
        for rate_change in sorted(grid)
    ]
    best_log_likelihood = max(row.log_likelihood for row in fits)
    best_rate_change = max(
        (
            row.rate_change
            for row in fits
            if math.isclose(
                row.log_likelihood,
                best_log_likelihood,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ),
        key=float,
    )
    rows: list[EarlyBurstRateChangeProfileRow] = []
    for row in fits:
        delta = best_log_likelihood - row.log_likelihood
        rows.append(
            EarlyBurstRateChangeProfileRow(
                rate_change=row.rate_change,
                log_likelihood=row.log_likelihood,
                aic=row.aic,
                aicc=row.aicc,
                delta_log_likelihood=delta,
                in_support_interval=delta <= _PROFILE_CONFIDENCE_DELTA,
                selected=math.isclose(
                    row.rate_change,
                    best_rate_change,
                    rel_tol=0.0,
                    abs_tol=1e-12,
                ),
            )
        )
    return rows


@dataclass(slots=True)
class _ProfileFitRow:
    rate_change: float
    log_likelihood: float
    aic: float
    aicc: float


def _fit_early_burst_candidate(
    dataset: ComparativeDataset,
    *,
    rate_change: float,
) -> _ProfileFitRow:
    transformed_tree = transform_tree_for_evolutionary_mode(
        dataset.tree,
        mode="early-burst",
        parameter_value=rate_change,
    )
    covariance = stable_covariance(
        build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
    )
    fit = _fit_intercept_only_model(dataset, covariance)
    comparison_row = _comparison_row(
        "early-burst",
        3,
        fit.log_likelihood,
        len(dataset.taxa),
    )
    return _ProfileFitRow(
        rate_change=rate_change,
        log_likelihood=fit.log_likelihood,
        aic=comparison_row.aic,
        aicc=comparison_row.aicc,
    )


def _build_rate_change_interval(
    profile_rows: list[EarlyBurstRateChangeProfileRow],
    estimate: float,
) -> ComparativeParameterInterval:
    supported = [row.rate_change for row in profile_rows if row.in_support_interval]
    lower = min(supported) if supported else estimate
    upper = max(supported) if supported else estimate
    return ComparativeParameterInterval(
        name="rate_change",
        estimate=estimate,
        lower_95=lower,
        upper_95=upper,
        method="profile_likelihood_95",
    )


def _detect_identifiability_warnings(
    *,
    rate_change: float,
    rate_change_bounds: tuple[float, float],
    profile_rows: list[EarlyBurstRateChangeProfileRow],
    comparison: ContinuousEvolutionaryModeComparisonReport,
) -> list[EarlyBurstIdentifiabilityWarning]:
    lower, upper = rate_change_bounds
    span = upper - lower
    interval_row = _build_rate_change_interval(profile_rows, rate_change)
    support_rows = [row for row in profile_rows if row.in_support_interval]
    warnings: list[EarlyBurstIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if (
        interval_row.lower_95 <= lower + boundary_tolerance
        or interval_row.upper_95 >= upper - boundary_tolerance
    ):
        warnings.append(
            EarlyBurstIdentifiabilityWarning(
                kind="boundary_rate_change",
                message=(
                    "the supported early-burst rate-change interval reaches a search "
                    "boundary, so the fitted rate-change parameter is weakly identified"
                ),
            )
        )
    if support_rows and (
        (interval_row.upper_95 - interval_row.lower_95) >= (0.75 * span)
        or len(support_rows) >= math.ceil(len(profile_rows) / 2.0)
    ):
        warnings.append(
            EarlyBurstIdentifiabilityWarning(
                kind="flat_likelihood_profile",
                message=(
                    "the early-burst likelihood profile stays broad across much of the "
                    "bounded search range, so rate-change estimates are weakly identified"
                ),
            )
        )
    best_aicc = min(row.aicc for row in comparison.rows)
    early_burst_row = next(row for row in comparison.rows if row.model == "early-burst")
    brownian_row = next(row for row in comparison.rows if row.model == "brownian")
    if (
        comparison.better_model != "early-burst"
        or (early_burst_row.aicc - best_aicc) <= 2.0
    ) and interval_row.lower_95 <= lower + boundary_tolerance:
        warnings.append(
            EarlyBurstIdentifiabilityWarning(
                kind="brownian_like_rate_change",
                message=(
                    "the early-burst fit remains close to the zero-change boundary and "
                    "Brownian AICc remains competitive, so time-varying rate claims are weak"
                ),
            )
        )
    if brownian_row.aicc < early_burst_row.aicc:
        warnings.append(
            EarlyBurstIdentifiabilityWarning(
                kind="comparison_not_preferred",
                message=(
                    "Brownian or OU remains preferred over early-burst by AICc on the "
                    "analyzed taxon set"
                ),
            )
        )
    return warnings


def _format_interval_bound(
    interval: ComparativeParameterInterval | None,
    field_name: str,
) -> str:
    if interval is None:
        return ""
    return format(getattr(interval, field_name), ".15g")

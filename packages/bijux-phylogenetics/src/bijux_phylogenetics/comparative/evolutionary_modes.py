from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    build_brownian_covariance_matrix,
    lambda_transform_covariance,
    load_comparative_dataset,
    node_signature,
    stable_covariance,
)
from bijux_phylogenetics.comparative.models import (
    ComparativeModelComparisonRow,
    ComparativeResidualSummary,
    _brownian_parameter_intervals,
    _build_residual_diagnostics,
    _comparison_row,
    _fit_intercept_only_model,
)
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree

ALLOWED_EVOLUTIONARY_MODES = {
    "brownian",
    "pagel-lambda",
    "pagel-kappa",
    "ornstein-uhlenbeck",
    "early-burst",
}


@dataclass(slots=True)
class EvolutionaryModeBranchLengthRow:
    """One deterministic branch-length change under a governed tree rescaling."""

    node: str
    descendant_taxa: list[str]
    original_branch_length: float
    transformed_branch_length: float
    parent_depth: float
    child_depth: float


@dataclass(slots=True)
class ComparativeTreeRescalingReport:
    """Canonical summary of one OU or early-burst tree-rescaling surface."""

    tree_path: Path
    mode: str
    parameter_name: str
    parameter_value: float
    tip_count: int
    original_total_branch_length: float
    transformed_total_branch_length: float
    transformed_tree_newick: str
    branch_rows: list[EvolutionaryModeBranchLengthRow]


@dataclass(slots=True)
class ContinuousEvolutionaryModeFitReport:
    """Intercept-only continuous-trait fit under one governed evolutionary mode."""

    tree_path: Path
    traits_path: Path
    taxon_column: str
    trait: str
    taxon_count: int
    taxa: list[str]
    mode: str
    parameter_name: str | None
    parameter_value: float | None
    root_state: float
    rate: float
    log_likelihood: float
    aic: float
    aicc: float
    fitted_values: list[float]
    residuals: list[float]
    transformed_tree_newick: str
    confidence_intervals: list[object]
    residual_diagnostics: ComparativeResidualSummary
    optimizer_diagnostics: ContinuousModeOptimizerDiagnostics | None
    identifiability_warnings: list[EvolutionaryModeIdentifiabilityWarning]
    assumptions: list[str]


@dataclass(slots=True)
class EvolutionaryModeIdentifiabilityWarning:
    """One identifiability or boundary warning on an evolutionary-mode fit."""

    kind: str
    message: str


@dataclass(slots=True)
class ContinuousModeOptimizerDiagnostics:
    """Optimizer diagnostics for one governed evolutionary-mode parameter search."""

    optimizer_name: str
    lower_bound: float
    upper_bound: float
    coarse_grid_point_count: int
    fine_grid_point_count: int
    function_evaluation_count: int
    coarse_best_parameter: float
    coarse_best_log_likelihood: float
    fine_search_start: float
    fine_search_stop: float
    converged: bool
    hit_lower_boundary: bool
    hit_upper_boundary: bool


@dataclass(slots=True)
class _TransformedModeSearchResult:
    parameter_value: float
    transformed_tree: PhyloTree
    covariance: list[list[float]]
    optimizer_diagnostics: ContinuousModeOptimizerDiagnostics
    profile: list[tuple[float, float]]


@dataclass(slots=True)
class LikelihoodRatioTestResult:
    """Likelihood-ratio comparison between two governed evolutionary-mode fits."""

    comparison_id: str
    left_mode: str
    right_mode: str
    statistic: float
    degrees_of_freedom: int
    p_value: float


@dataclass(slots=True)
class ContinuousEvolutionaryModeComparisonReport:
    """Model-comparison summary over Brownian, OU, and early-burst fits."""

    tree_path: Path
    traits_path: Path
    trait: str
    taxon_count: int
    rows: list[ComparativeModelComparisonRow]
    better_model: str
    likelihood_ratio_tests: list[LikelihoodRatioTestResult]


def rescale_tree_ornstein_uhlenbeck(
    tree_path: Path,
    *,
    alpha: float,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style OU branch rescaling to a rooted tree."""
    return _build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="ornstein-uhlenbeck",
        parameter_name="alpha",
        parameter_value=alpha,
        sigsq=sigsq,
    )


def rescale_tree_early_burst(
    tree_path: Path,
    *,
    rate_change: float,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style early-burst branch rescaling to a rooted tree."""
    return _build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="early-burst",
        parameter_name="rate_change",
        parameter_value=rate_change,
        sigsq=sigsq,
    )


def rescale_tree_pagel_kappa(
    tree_path: Path,
    *,
    kappa: float,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style Pagel-kappa branch-length rescaling to a rooted tree."""
    return _build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="pagel-kappa",
        parameter_name="kappa",
        parameter_value=kappa,
        sigsq=1.0,
    )


def transform_tree_for_evolutionary_mode(
    tree: PhyloTree,
    *,
    mode: str,
    parameter_value: float,
    sigsq: float = 1.0,
) -> PhyloTree:
    """Transform an in-memory tree under a governed continuous-mode branch rule."""
    return _transform_tree(
        tree,
        mode=mode,
        parameter_value=parameter_value,
        sigsq=sigsq,
    )


def fit_continuous_evolutionary_mode(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    mode: str,
    taxon_column: str | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeFitReport:
    """Fit a Brownian, Pagel-lambda, Pagel-kappa, OU, or early-burst intercept-only trait model."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    return _fit_evolutionary_mode_from_dataset(
        dataset,
        mode=mode,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def compare_continuous_evolutionary_modes(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare Brownian, OU-rescaled, and early-burst-rescaled intercept-only fits."""
    dataset = load_comparative_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    fits = [
        _fit_evolutionary_mode_from_dataset(
            dataset,
            mode=mode,
            lambda_bounds=lambda_bounds,
            kappa_bounds=kappa_bounds,
            ou_bounds=ou_bounds,
            early_burst_bounds=early_burst_bounds,
        )
        for mode in ("brownian", "ornstein-uhlenbeck", "early-burst")
    ]
    rows = [
        _comparison_row(
            fit.mode,
            2 if fit.parameter_value is None else 3,
            fit.log_likelihood,
            fit.taxon_count,
        )
        for fit in fits
    ]
    best_aicc = min(row.aicc for row in rows)
    for row in rows:
        row.selected = math.isclose(row.aicc, best_aicc, rel_tol=0.0, abs_tol=1e-12)
    likelihood_ratio_tests = [
        _likelihood_ratio_test(
            comparison_id="brownian-vs-ornstein-uhlenbeck",
            left_fit=fits[0],
            right_fit=fits[1],
        ),
        _likelihood_ratio_test(
            comparison_id="brownian-vs-early-burst",
            left_fit=fits[0],
            right_fit=fits[2],
        ),
        _likelihood_ratio_test(
            comparison_id="ornstein-uhlenbeck-vs-early-burst",
            left_fit=fits[1],
            right_fit=fits[2],
        ),
    ]
    better_model = next(row.model for row in rows if row.selected)
    return ContinuousEvolutionaryModeComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=dataset.readiness.tree_taxa,
        rows=rows,
        better_model=better_model,
        likelihood_ratio_tests=likelihood_ratio_tests,
    )


def _fit_evolutionary_mode_from_dataset(
    dataset: ComparativeDataset,
    *,
    mode: str,
    lambda_bounds: tuple[float, float],
    kappa_bounds: tuple[float, float],
    ou_bounds: tuple[float, float],
    early_burst_bounds: tuple[float, float],
) -> ContinuousEvolutionaryModeFitReport:
    if mode not in ALLOWED_EVOLUTIONARY_MODES:
        raise ComparativeMethodError(
            "unsupported evolutionary mode; expected one of: "
            + ", ".join(sorted(ALLOWED_EVOLUTIONARY_MODES))
        )

    if mode == "brownian":
        transformed_tree = _clone_tree(dataset.tree)
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        parameter_name = None
        parameter_value = None
        optimizer_diagnostics = None
        identifiability_warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
        assumptions = [
            "Brownian mode retains the original rooted branch lengths.",
            "Trait variance accumulates proportionally with shared branch length.",
        ]
    elif mode == "pagel-lambda":
        parameter_name = "lambda"
        search_result = _best_pagel_lambda_fit(
            dataset,
            bounds=lambda_bounds,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        identifiability_warnings = _lambda_identifiability_warnings_from_profile(
            parameter_value,
            search_result.profile,
            lambda_bounds,
        )
        assumptions = [
            "Pagel-lambda mode follows the geiger-style covariance transformation, scaling shared covariance while keeping each tip variance fixed.",
            "Unlike phytools::phylosig-style signal summaries, this mode keeps the full fitContinuous surface visible with root state, sigma-squared-backed rate, AIC, and AICc.",
        ]
    elif mode == "pagel-kappa":
        parameter_name = "kappa"
        search_result = _best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=kappa_bounds,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        identifiability_warnings = _kappa_identifiability_warnings_from_profile(
            parameter_value,
            search_result.profile,
            kappa_bounds,
        )
        assumptions = [
            "Pagel-kappa mode follows the geiger-style branch-length power transformation, raising each branch length to kappa before Brownian intercept fitting.",
            "The fitted covariance comes from the transformed branch-length tree rather than from a topology-only change count surface.",
        ]
    elif mode == "ornstein-uhlenbeck":
        parameter_name = "alpha"
        search_result = _best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=ou_bounds,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        identifiability_warnings = _ou_identifiability_warnings_from_profile(
            dataset,
            parameter_value,
            search_result.profile,
        )
        assumptions = [
            "OU mode follows the lecture-style tree rescaling before Brownian intercept fitting.",
            "Alpha is selected by maximizing log likelihood over a governed bounded search grid.",
        ]
    else:
        parameter_name = "rate_change"
        search_result = _best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=early_burst_bounds,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        identifiability_warnings = _early_burst_identifiability_warnings_from_profile(
            parameter_value,
            search_result.profile,
            early_burst_bounds,
        )
        assumptions = [
            "Early-burst mode follows the geiger-style branch-rescaling rule before Brownian intercept fitting.",
            "The rate-change parameter is selected by maximizing log likelihood over a governed bounded search grid.",
        ]

    row = _comparison_row(
        mode,
        2 if parameter_value is None else 3,
        fit.log_likelihood,
        len(dataset.taxa),
    )
    intervals = _brownian_parameter_intervals(
        fit.theta,
        fit.sigma_squared,
        covariance,
    )
    if parameter_name is not None and parameter_value is not None:
        intervals = [
            *[
                interval
                for interval in intervals
                if getattr(interval, "name", None) != "root_state"
            ],
        ]
    residual_diagnostics = _build_residual_diagnostics(
        dataset,
        covariance,
        fit.residuals,
        fit.sigma_squared,
    )
    return ContinuousEvolutionaryModeFitReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=_stable_value(parameter_value)
        if parameter_value is not None
        else None,
        root_state=_stable_value(fit.theta),
        rate=_stable_value(fit.sigma_squared),
        log_likelihood=_stable_value(fit.log_likelihood),
        aic=_stable_value(row.aic),
        aicc=_stable_value(row.aicc),
        fitted_values=[_stable_value(value) for value in fit.fitted_values],
        residuals=[_stable_value(value) for value in fit.residuals],
        transformed_tree_newick=dumps_newick(transformed_tree),
        confidence_intervals=intervals,
        residual_diagnostics=residual_diagnostics,
        optimizer_diagnostics=optimizer_diagnostics,
        identifiability_warnings=identifiability_warnings,
        assumptions=assumptions,
    )


def _best_transformed_mode_fit(
    dataset: ComparativeDataset,
    *,
    mode: str,
    bounds: tuple[float, float],
) -> _TransformedModeSearchResult:
    lower, upper = bounds
    if upper <= lower:
        raise ComparativeMethodError("parameter bounds must be strictly increasing")
    if mode == "ornstein-uhlenbeck":
        lower = max(lower, 1e-6)

    coarse = _linspace(lower, upper, 81)
    best_parameter = coarse[0]
    best_tree = _transform_tree(dataset.tree, mode=mode, parameter_value=best_parameter)
    best_covariance = stable_covariance(
        build_brownian_covariance_matrix(best_tree, dataset.taxa)
    )
    best_fit = _fit_intercept_only_model(dataset, best_covariance)
    best_index = 0
    profile: list[tuple[float, float]] = [(best_parameter, best_fit.log_likelihood)]
    coarse_best_parameter = best_parameter
    coarse_best_log_likelihood = best_fit.log_likelihood
    for index, candidate in enumerate(coarse[1:], start=1):
        transformed_tree = _transform_tree(
            dataset.tree,
            mode=mode,
            parameter_value=candidate,
        )
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        profile.append((candidate, fit.log_likelihood))
        if fit.log_likelihood > best_fit.log_likelihood:
            best_parameter = candidate
            best_tree = transformed_tree
            best_covariance = covariance
            best_fit = fit
            best_index = index
            coarse_best_parameter = candidate
            coarse_best_log_likelihood = fit.log_likelihood

    left = coarse[max(0, best_index - 1)]
    right = coarse[min(len(coarse) - 1, best_index + 1)]
    fine_candidates = _linspace(left, right, 81)
    for candidate in _linspace(left, right, 81):
        if math.isclose(candidate, best_parameter, rel_tol=0.0, abs_tol=1e-12):
            continue
        transformed_tree = _transform_tree(
            dataset.tree,
            mode=mode,
            parameter_value=candidate,
        )
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        profile.append((candidate, fit.log_likelihood))
        if fit.log_likelihood > best_fit.log_likelihood:
            best_parameter = candidate
            best_tree = transformed_tree
            best_covariance = covariance
            best_fit = fit
    fine_step = 0.0
    if len(fine_candidates) > 1:
        fine_step = fine_candidates[1] - fine_candidates[0]
    tolerance = max(abs(fine_step) / 2.0, 1e-9)
    diagnostics = ContinuousModeOptimizerDiagnostics(
        optimizer_name="governed-two-stage-grid-search",
        lower_bound=_stable_value(lower),
        upper_bound=_stable_value(upper),
        coarse_grid_point_count=len(coarse),
        fine_grid_point_count=len(fine_candidates),
        function_evaluation_count=len(profile),
        coarse_best_parameter=_stable_value(coarse_best_parameter),
        coarse_best_log_likelihood=_stable_value(coarse_best_log_likelihood),
        fine_search_start=_stable_value(left),
        fine_search_stop=_stable_value(right),
        converged=True,
        hit_lower_boundary=math.isclose(
            best_parameter,
            lower,
            rel_tol=0.0,
            abs_tol=tolerance,
        ),
        hit_upper_boundary=math.isclose(
            best_parameter,
            upper,
            rel_tol=0.0,
            abs_tol=tolerance,
        ),
    )
    normalized_profile = sorted(
        {round(parameter, 12): (parameter, log_likelihood) for parameter, log_likelihood in profile}.values(),
        key=lambda item: item[0],
    )
    return _TransformedModeSearchResult(
        parameter_value=best_parameter,
        transformed_tree=best_tree,
        covariance=best_covariance,
        optimizer_diagnostics=diagnostics,
        profile=normalized_profile,
    )


def _best_pagel_lambda_fit(
    dataset: ComparativeDataset,
    *,
    bounds: tuple[float, float],
) -> _TransformedModeSearchResult:
    lower, upper = bounds
    if lower < 0.0 or upper > 1.0 or upper <= lower:
        raise ComparativeMethodError(
            "Pagel-lambda bounds must be strictly increasing within [0, 1]"
        )

    coarse = _linspace(lower, upper, 81)
    best_parameter = coarse[0]
    best_tree = _transform_tree(
        dataset.tree,
        mode="pagel-lambda",
        parameter_value=best_parameter,
    )
    best_covariance = lambda_transform_covariance(
        dataset.covariance_matrix,
        best_parameter,
    )
    best_fit = _fit_intercept_only_model(dataset, best_covariance)
    best_index = 0
    profile: list[tuple[float, float]] = [(best_parameter, best_fit.log_likelihood)]
    coarse_best_parameter = best_parameter
    coarse_best_log_likelihood = best_fit.log_likelihood
    for index, candidate in enumerate(coarse[1:], start=1):
        transformed_tree = _transform_tree(
            dataset.tree,
            mode="pagel-lambda",
            parameter_value=candidate,
        )
        covariance = lambda_transform_covariance(
            dataset.covariance_matrix,
            candidate,
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        profile.append((candidate, fit.log_likelihood))
        if fit.log_likelihood > best_fit.log_likelihood:
            best_parameter = candidate
            best_tree = transformed_tree
            best_covariance = covariance
            best_fit = fit
            best_index = index
            coarse_best_parameter = candidate
            coarse_best_log_likelihood = fit.log_likelihood

    left = coarse[max(0, best_index - 1)]
    right = coarse[min(len(coarse) - 1, best_index + 1)]
    fine_candidates = _linspace(left, right, 81)
    for candidate in fine_candidates:
        if math.isclose(candidate, best_parameter, rel_tol=0.0, abs_tol=1e-12):
            continue
        transformed_tree = _transform_tree(
            dataset.tree,
            mode="pagel-lambda",
            parameter_value=candidate,
        )
        covariance = lambda_transform_covariance(
            dataset.covariance_matrix,
            candidate,
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        profile.append((candidate, fit.log_likelihood))
        if fit.log_likelihood > best_fit.log_likelihood:
            best_parameter = candidate
            best_tree = transformed_tree
            best_covariance = covariance
            best_fit = fit

    fine_step = 0.0
    if len(fine_candidates) > 1:
        fine_step = fine_candidates[1] - fine_candidates[0]
    tolerance = max(abs(fine_step) / 2.0, 1e-9)
    diagnostics = ContinuousModeOptimizerDiagnostics(
        optimizer_name="governed-two-stage-grid-search",
        lower_bound=_stable_value(lower),
        upper_bound=_stable_value(upper),
        coarse_grid_point_count=len(coarse),
        fine_grid_point_count=len(fine_candidates),
        function_evaluation_count=len(profile),
        coarse_best_parameter=_stable_value(coarse_best_parameter),
        coarse_best_log_likelihood=_stable_value(coarse_best_log_likelihood),
        fine_search_start=_stable_value(left),
        fine_search_stop=_stable_value(right),
        converged=True,
        hit_lower_boundary=math.isclose(
            best_parameter,
            lower,
            rel_tol=0.0,
            abs_tol=tolerance,
        ),
        hit_upper_boundary=math.isclose(
            best_parameter,
            upper,
            rel_tol=0.0,
            abs_tol=tolerance,
        ),
    )
    normalized_profile = sorted(
        {
            round(parameter, 12): (parameter, log_likelihood)
            for parameter, log_likelihood in profile
        }.values(),
        key=lambda item: item[0],
    )
    return _TransformedModeSearchResult(
        parameter_value=best_parameter,
        transformed_tree=best_tree,
        covariance=best_covariance,
        optimizer_diagnostics=diagnostics,
        profile=normalized_profile,
    )


def _linspace(start: float, stop: float, count: int) -> list[float]:
    if count < 2:
        return [start]
    step = (stop - start) / float(count - 1)
    return [start + (step * index) for index in range(count)]


def _ou_identifiability_warnings_from_profile(
    dataset: ComparativeDataset,
    alpha: float,
    profile: list[tuple[float, float]],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    if len(dataset.taxa) < 5:
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="small_sample_size",
                message="OU alpha is hard to identify with fewer than five taxa",
            )
        )
    ordered_alphas = sorted(candidate for candidate, _ in profile)
    if math.isclose(
        alpha,
        ordered_alphas[0],
        rel_tol=0.0,
        abs_tol=1e-9,
    ) or math.isclose(
        alpha,
        ordered_alphas[-1],
        rel_tol=0.0,
        abs_tol=1e-9,
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_alpha",
                message="best-supported OU alpha falls on the search boundary and may not be well identified",
            )
        )
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="OU likelihood surface is shallow across alpha values, so model choice may be unstable",
            )
        )
    if alpha < ordered_alphas[len(ordered_alphas) // 3]:
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="weak_pull_to_optimum",
                message="best-supported OU alpha is weak and may be difficult to distinguish from Brownian motion",
            )
        )
    return warnings


def _early_burst_identifiability_warnings_from_profile(
    rate_change: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        rate_change,
        lower,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    ) or math.isclose(
        rate_change,
        upper,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_rate_change",
                message="best-supported early-burst rate change falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood_profile",
                message="early-burst likelihood stays shallow across rate-change values, so model choice may be unstable",
            )
        )
    if rate_change <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="brownian_like_rate_change",
                message="best-supported early-burst rate change remains close to the zero-change boundary and may be difficult to distinguish from Brownian motion",
            )
        )
    return warnings


def _lambda_identifiability_warnings_from_profile(
    lambda_value: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        lambda_value,
        lower,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    ) or math.isclose(
        lambda_value,
        upper,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_lambda",
                message="best-supported Pagel lambda falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="Pagel-lambda likelihood stays shallow across the bounded search, so the covariance scaling may be unstable",
            )
        )
    if lambda_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="weak_phylogenetic_signal",
                message="best-supported Pagel lambda remains close to the zero-signal boundary and may be difficult to distinguish from a star-like covariance surface",
            )
        )
    if lambda_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="brownian_limit",
                message="best-supported Pagel lambda remains close to the Brownian boundary and may be difficult to distinguish from an untransformed covariance surface",
            )
        )
    return warnings


def _kappa_identifiability_warnings_from_profile(
    kappa_value: float,
    profile: list[tuple[float, float]],
    bounds: tuple[float, float],
) -> list[EvolutionaryModeIdentifiabilityWarning]:
    lower, upper = bounds
    span = upper - lower
    ordered_log_likelihoods = sorted(
        (log_likelihood for _, log_likelihood in profile),
        reverse=True,
    )
    warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
    boundary_tolerance = max(span / 160.0, 1e-9)
    if math.isclose(
        kappa_value,
        lower,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    ) or math.isclose(
        kappa_value,
        upper,
        rel_tol=0.0,
        abs_tol=boundary_tolerance,
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="boundary_kappa",
                message="best-supported Pagel kappa falls on the search boundary and may not be well identified",
            )
        )
    if len(ordered_log_likelihoods) > 1 and (
        ordered_log_likelihoods[0] - ordered_log_likelihoods[1] < 0.5
    ):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="flat_likelihood",
                message="Pagel-kappa likelihood stays shallow across the bounded search, so branch-length transformation support may be unstable",
            )
        )
    if kappa_value <= lower + max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="punctuational_limit",
                message="best-supported Pagel kappa remains close to the equal-length punctuational boundary and may be difficult to distinguish from a branch-count surface",
            )
        )
    if kappa_value >= upper - max(boundary_tolerance, 1e-6):
        warnings.append(
            EvolutionaryModeIdentifiabilityWarning(
                kind="upper_search_limit",
                message="best-supported Pagel kappa remains close to the upper search boundary, so a wider branch-length review may be needed",
            )
        )
    return warnings


def _build_tree_rescaling_report(
    tree: PhyloTree,
    tree_path: Path,
    *,
    mode: str,
    parameter_name: str,
    parameter_value: float,
    sigsq: float,
) -> ComparativeTreeRescalingReport:
    transformed_tree = _transform_tree(
        tree,
        mode=mode,
        parameter_value=parameter_value,
        sigsq=sigsq,
    )
    branch_rows = _branch_length_rows(
        original_tree=tree,
        transformed_tree=transformed_tree,
    )
    return ComparativeTreeRescalingReport(
        tree_path=tree_path,
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=_stable_value(parameter_value),
        tip_count=tree.tip_count,
        original_total_branch_length=_stable_value(tree.total_branch_length()),
        transformed_total_branch_length=_stable_value(
            transformed_tree.total_branch_length()
        ),
        transformed_tree_newick=dumps_newick(transformed_tree),
        branch_rows=branch_rows,
    )


def _transform_tree(
    tree: PhyloTree,
    *,
    mode: str,
    parameter_value: float,
    sigsq: float = 1.0,
) -> PhyloTree:
    if mode not in {
        "ornstein-uhlenbeck",
        "early-burst",
        "pagel-lambda",
        "pagel-kappa",
    }:
        raise ComparativeMethodError(
            "tree transformation mode must be 'ornstein-uhlenbeck', 'early-burst', 'pagel-lambda', or 'pagel-kappa'"
        )
    if mode == "ornstein-uhlenbeck" and parameter_value < 0.0:
        raise ComparativeMethodError("OU alpha must be non-negative")
    if mode == "pagel-lambda" and not 0.0 <= parameter_value <= 1.0:
        raise ComparativeMethodError("Pagel lambda must lie within [0, 1]")
    if mode == "pagel-kappa" and parameter_value < 0.0:
        raise ComparativeMethodError("Pagel kappa must be non-negative")
    cloned_root = _clone_node(tree.root)
    if mode == "pagel-lambda":
        def visit_pagel_lambda(node: TreeNode, depth: float) -> None:
            for child in node.children:
                original_length = float(child.branch_length or 0.0)
                if child.is_leaf():
                    child.branch_length = original_length + (
                        (1.0 - parameter_value) * depth
                    )
                else:
                    child.branch_length = original_length * parameter_value
                visit_pagel_lambda(child, depth + original_length)

        visit_pagel_lambda(cloned_root, 0.0)
        return PhyloTree(
            root=cloned_root,
            source_format=tree.source_format,
            rooted=tree.rooted,
        )
    if mode == "pagel-kappa":
        def visit_pagel_kappa(node: TreeNode) -> None:
            for child in node.children:
                original_length = float(child.branch_length or 0.0)
                child.branch_length = _kappa_branch_length(
                    original_length,
                    kappa=parameter_value,
                )
                visit_pagel_kappa(child)

        visit_pagel_kappa(cloned_root)
        return PhyloTree(
            root=cloned_root,
            source_format=tree.source_format,
            rooted=tree.rooted,
        )

    total_depth = _max_tip_depth(tree.root, depth=0.0)

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            original_length = float(child.branch_length or 0.0)
            child_depth = depth + original_length
            if mode == "ornstein-uhlenbeck":
                child.branch_length = _ou_branch_length(
                    parent_depth=depth,
                    child_depth=child_depth,
                    total_depth=total_depth,
                    alpha=parameter_value,
                    sigsq=sigsq,
                )
            else:
                child.branch_length = _early_burst_branch_length(
                    parent_depth=depth,
                    child_depth=child_depth,
                    rate_change=parameter_value,
                    sigsq=sigsq,
                )
            visit(child, child_depth)

    visit(cloned_root, 0.0)
    return PhyloTree(
        root=cloned_root,
        source_format=tree.source_format,
        rooted=tree.rooted,
    )


def _branch_length_rows(
    *,
    original_tree: PhyloTree,
    transformed_tree: PhyloTree,
) -> list[EvolutionaryModeBranchLengthRow]:
    original_rows = _tree_branch_lookup(original_tree)
    transformed_rows = _tree_branch_lookup(transformed_tree)
    rows: list[EvolutionaryModeBranchLengthRow] = []
    for node_id in sorted(original_rows):
        original = original_rows[node_id]
        transformed = transformed_rows[node_id]
        rows.append(
            EvolutionaryModeBranchLengthRow(
                node=node_id,
                descendant_taxa=list(original["descendant_taxa"]),
                original_branch_length=_stable_value(original["branch_length"]),
                transformed_branch_length=_stable_value(transformed["branch_length"]),
                parent_depth=_stable_value(original["parent_depth"]),
                child_depth=_stable_value(original["child_depth"]),
            )
        )
    return rows


def _tree_branch_lookup(tree: PhyloTree) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            branch_length = float(child.branch_length or 0.0)
            child_depth = depth + branch_length
            branch_id = node_signature(child)
            rows[branch_id] = {
                "branch_length": branch_length,
                "parent_depth": depth,
                "child_depth": child_depth,
                "descendant_taxa": _descendant_taxa(child),
            }
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return rows


def _ou_branch_length(
    *,
    parent_depth: float,
    child_depth: float,
    total_depth: float,
    alpha: float,
    sigsq: float,
) -> float:
    if alpha <= 0.0:
        raise ComparativeMethodError("OU alpha must be positive")

    def _term(depth: float) -> float:
        return (
            (1.0 / (2.0 * alpha))
            * math.exp(-2.0 * alpha * (total_depth - depth))
            * (1.0 - math.exp(-2.0 * alpha * depth))
        )

    return max(0.0, (_term(child_depth) - _term(parent_depth)) * sigsq)


def _early_burst_branch_length(
    *,
    parent_depth: float,
    child_depth: float,
    rate_change: float,
    sigsq: float,
) -> float:
    if math.isclose(rate_change, 0.0, rel_tol=0.0, abs_tol=1e-12):
        return max(0.0, (child_depth - parent_depth) * sigsq)
    transformed = (
        math.exp(-rate_change * parent_depth)
        - math.exp(-rate_change * child_depth)
    ) / rate_change
    return max(0.0, transformed * sigsq)


def _kappa_branch_length(branch_length: float, *, kappa: float) -> float:
    if kappa < 0.0:
        raise ComparativeMethodError("Pagel kappa must be non-negative")
    if branch_length < 0.0:
        raise ComparativeMethodError(
            "Pagel kappa cannot transform negative branch lengths"
        )
    transformed = math.pow(branch_length, kappa)
    if not math.isfinite(transformed) or transformed < 0.0:
        raise ComparativeMethodError(
            "Pagel kappa produced an invalid transformed branch length"
        )
    return transformed


def _likelihood_ratio_test(
    *,
    comparison_id: str,
    left_fit: ContinuousEvolutionaryModeFitReport,
    right_fit: ContinuousEvolutionaryModeFitReport,
) -> LikelihoodRatioTestResult:
    statistic = max(0.0, -2.0 * (left_fit.log_likelihood - right_fit.log_likelihood))
    return LikelihoodRatioTestResult(
        comparison_id=comparison_id,
        left_mode=left_fit.mode,
        right_mode=right_fit.mode,
        statistic=_stable_value(statistic),
        degrees_of_freedom=1,
        p_value=_stable_value(math.erfc(math.sqrt(statistic / 2.0))),
    )


def _clone_node(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=[_clone_node(child) for child in node.children],
    )


def _clone_tree(tree: PhyloTree) -> PhyloTree:
    return PhyloTree(
        root=_clone_node(tree.root),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )


def _max_tip_depth(node: TreeNode, *, depth: float) -> float:
    if node.is_leaf():
        return depth
    return max(
        _max_tip_depth(
            child,
            depth=depth + float(child.branch_length or 0.0),
        )
        for child in node.children
    )


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _stable_value(value: float | None) -> float:
    if value is None:
        raise ValueError("expected a float value")
    return float(format(round(float(value), 15), ".15g"))

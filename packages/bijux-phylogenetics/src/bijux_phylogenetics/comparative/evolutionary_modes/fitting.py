from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    build_brownian_covariance_matrix,
    load_comparative_dataset,
    stable_covariance,
)
from bijux_phylogenetics.comparative.continuous.model_fitting import (
    _brownian_parameter_intervals,
    _build_residual_diagnostics,
    _comparison_row,
    _fit_intercept_only_model,
)
from bijux_phylogenetics.comparative.evolutionary_modes.models import (
    ALLOWED_EVOLUTIONARY_MODES,
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY,
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    EXCLUDED_GEIGER_TREND_MODE_ALIASES,
    ContinuousEvolutionaryModeFitReport,
    ContinuousModeOptimizerProfileRow,
    ContinuousModeSearchControls,
    EvolutionaryModeIdentifiabilityWarning,
)
from bijux_phylogenetics.comparative.evolutionary_modes.numeric import stable_float
from bijux_phylogenetics.comparative.evolutionary_modes.search import (
    best_pagel_lambda_fit,
    best_transformed_mode_fit,
    continuous_boundary_assessment,
    delta_identifiability_warnings_from_profile,
    early_burst_identifiability_warnings_from_profile,
    kappa_identifiability_warnings_from_profile,
    lambda_identifiability_warnings_from_profile,
    normalized_search_controls,
    ou_identifiability_warnings_from_profile,
    reject_nonparameterized_search_controls,
)
from bijux_phylogenetics.comparative.evolutionary_modes.tree_transforms import (
    clone_tree,
    identity_covariance_matrix,
    transform_tree,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


def fit_continuous_evolutionary_mode(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    mode: str,
    taxon_column: str | None = None,
    standard_error_trait: str | None = None,
    search_controls: ContinuousModeSearchControls | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    delta_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeFitReport:
    """Fit a governed intercept-only continuous-trait evolutionary mode."""
    if standard_error_trait is not None:
        raise ComparativeMethodError(
            "geiger::fitContinuous standard-error parity is explicitly excluded in this round; "
            f"received standard_error_trait='{standard_error_trait}' while Bijux still lacks "
            "measurement-error variance handling on the owned continuous-mode fit surface"
        )
    if mode in EXCLUDED_GEIGER_TREND_MODE_ALIASES:
        raise ComparativeMethodError(
            "trend-mode parity is explicitly excluded in this round because geiger exposes distinct `rate_trend` and `mean_trend` likelihoods rather than one unambiguous `trend` contract"
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
    return fit_evolutionary_mode_from_dataset(
        dataset,
        mode=mode,
        search_controls=search_controls,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def fit_evolutionary_mode_from_dataset(
    dataset: ComparativeDataset,
    *,
    mode: str,
    search_controls: ContinuousModeSearchControls | None = None,
    lambda_bounds: tuple[float, float],
    kappa_bounds: tuple[float, float],
    delta_bounds: tuple[float, float],
    ou_bounds: tuple[float, float],
    early_burst_bounds: tuple[float, float],
) -> ContinuousEvolutionaryModeFitReport:
    if mode not in ALLOWED_EVOLUTIONARY_MODES:
        raise ComparativeMethodError(
            "unsupported evolutionary mode; expected one of: "
            + ", ".join(sorted(ALLOWED_EVOLUTIONARY_MODES))
        )

    if mode in {"brownian", "white-noise"}:
        reject_nonparameterized_search_controls(mode, search_controls)
    else:
        search_controls = normalized_search_controls(search_controls)

    if mode == "brownian":
        transformed_tree = clone_tree(dataset.tree)
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
        fit = _fit_intercept_only_model(dataset, covariance)
        parameter_name = None
        parameter_value = None
        optimizer_diagnostics = None
        optimizer_profile_rows = None
        identifiability_warnings: list[EvolutionaryModeIdentifiabilityWarning] = []
        assumptions = [
            "Brownian mode retains the original rooted branch lengths.",
            "Trait variance accumulates proportionally with shared branch length.",
        ]
    elif mode == "white-noise":
        transformed_tree = transform_tree(
            dataset.tree,
            mode="white-noise",
            parameter_value=1.0,
            sigsq=1.0,
        )
        covariance = stable_covariance(identity_covariance_matrix(len(dataset.taxa)))
        fit = _fit_intercept_only_model(dataset, covariance)
        parameter_name = None
        parameter_value = None
        optimizer_diagnostics = None
        optimizer_profile_rows = None
        identifiability_warnings = [
            EvolutionaryModeIdentifiabilityWarning(
                kind="no_phylogenetic_correlation",
                message="white-noise mode assumes independent residual variance across taxa and ignores shared phylogenetic covariance",
            )
        ]
        assumptions = [
            "White-noise mode treats taxa as independent and uses an identity covariance surface instead of branch-length-derived shared variance.",
            "This mode is the no-phylogenetic-correlation baseline for comparing whether Brownian or transformed phylogenetic models explain trait covariance better.",
        ]
    elif mode == "pagel-lambda":
        parameter_name = "lambda"
        search_result = best_pagel_lambda_fit(
            dataset,
            bounds=lambda_bounds,
            search_controls=search_controls,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        optimizer_profile_rows = [
            ContinuousModeOptimizerProfileRow(
                parameter_value=stable_float(parameter),
                log_likelihood=stable_float(log_likelihood),
            )
            for parameter, log_likelihood in search_result.profile
        ]
        identifiability_warnings = lambda_identifiability_warnings_from_profile(
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
        search_result = best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=kappa_bounds,
            search_controls=search_controls,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        optimizer_profile_rows = [
            ContinuousModeOptimizerProfileRow(
                parameter_value=stable_float(parameter),
                log_likelihood=stable_float(log_likelihood),
            )
            for parameter, log_likelihood in search_result.profile
        ]
        identifiability_warnings = kappa_identifiability_warnings_from_profile(
            parameter_value,
            search_result.profile,
            kappa_bounds,
        )
        assumptions = [
            "Pagel-kappa mode follows the geiger-style branch-length power transformation, raising each branch length to kappa before Brownian intercept fitting.",
            "The fitted covariance comes from the transformed branch-length tree rather than from a topology-only change count surface.",
        ]
    elif mode == "pagel-delta":
        parameter_name = "delta"
        search_result = best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=delta_bounds,
            search_controls=search_controls,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        optimizer_profile_rows = [
            ContinuousModeOptimizerProfileRow(
                parameter_value=stable_float(parameter),
                log_likelihood=stable_float(log_likelihood),
            )
            for parameter, log_likelihood in search_result.profile
        ]
        identifiability_warnings = delta_identifiability_warnings_from_profile(
            parameter_value,
            search_result.profile,
            delta_bounds,
        )
        assumptions = [
            "Pagel-delta mode follows the geiger-style depth transformation, raising each node depth proportion to delta before recomputing branch lengths.",
            "The fitted covariance comes from the transformed depth-scaled tree rather than from an edge-wise power transform.",
        ]
    elif mode == "ornstein-uhlenbeck":
        parameter_name = "alpha"
        search_result = best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=ou_bounds,
            search_controls=search_controls,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        optimizer_profile_rows = [
            ContinuousModeOptimizerProfileRow(
                parameter_value=stable_float(parameter),
                log_likelihood=stable_float(log_likelihood),
            )
            for parameter, log_likelihood in search_result.profile
        ]
        identifiability_warnings = ou_identifiability_warnings_from_profile(
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
        search_result = best_transformed_mode_fit(
            dataset,
            mode=mode,
            bounds=early_burst_bounds,
            search_controls=search_controls,
        )
        parameter_value = search_result.parameter_value
        transformed_tree = search_result.transformed_tree
        covariance = search_result.covariance
        fit = _fit_intercept_only_model(dataset, covariance)
        optimizer_diagnostics = search_result.optimizer_diagnostics
        optimizer_profile_rows = [
            ContinuousModeOptimizerProfileRow(
                parameter_value=stable_float(parameter),
                log_likelihood=stable_float(log_likelihood),
            )
            for parameter, log_likelihood in search_result.profile
        ]
        identifiability_warnings = early_burst_identifiability_warnings_from_profile(
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
        likelihood_constant_policy=CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
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
    boundary_assessment = continuous_boundary_assessment(
        parameter_name=parameter_name,
        parameter_value=parameter_value,
        optimizer_diagnostics=optimizer_diagnostics,
        identifiability_warnings=identifiability_warnings,
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
        parameter_value=stable_float(parameter_value)
        if parameter_value is not None
        else None,
        root_state=stable_float(fit.theta),
        rate=stable_float(fit.sigma_squared),
        log_likelihood=stable_float(fit.log_likelihood),
        aic=stable_float(row.aic),
        aicc=stable_float(row.aicc),
        likelihood_constant_policy=CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
        likelihood_comparison_policy=CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY,
        fitted_values=[stable_float(value) for value in fit.fitted_values],
        residuals=[stable_float(value) for value in fit.residuals],
        transformed_tree_newick=dumps_newick(transformed_tree),
        confidence_intervals=intervals,
        residual_diagnostics=residual_diagnostics,
        optimizer_diagnostics=optimizer_diagnostics,
        optimizer_profile_rows=optimizer_profile_rows,
        identifiability_warnings=identifiability_warnings,
        assumptions=assumptions,
        boundary_assessment=boundary_assessment,
    )

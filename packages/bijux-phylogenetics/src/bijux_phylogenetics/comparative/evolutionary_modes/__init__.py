from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import (
    ComparativeDataset,
    load_comparative_dataset,
)
from bijux_phylogenetics.comparative.information_criteria import (
    rank_model_comparison_rows,
)
from bijux_phylogenetics.comparative.evolutionary_modes.models import (
    ALLOWED_EVOLUTIONARY_MODES,
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_COMPARISON_POLICY,
    CONTINUOUS_GAUSSIAN_LIKELIHOOD_CONSTANT_POLICY,
    EXCLUDED_GEIGER_STANDARD_ERROR_POLICY,
    EXCLUDED_GEIGER_TREND_MODE_ALIASES,
    FITCONTINUOUS_MODEL_COMPARISON_ORDER,
    FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
    FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS,
    FITCONTINUOUS_MODEL_RANKING_POLICY,
    ComparativeTreeRescalingReport,
    ContinuousEvolutionaryModeComparisonReport,
    ContinuousEvolutionaryModeFitReport,
    ContinuousModeBoundaryAssessment,
    ContinuousModeOptimizerDiagnostics,
    ContinuousModeOptimizerProfileRow,
    ContinuousModeSearchControls,
    EvolutionaryModeBranchLengthRow,
    EvolutionaryModeIdentifiabilityWarning,
    LikelihoodRatioTestResult,
)
from bijux_phylogenetics.comparative.evolutionary_modes.numeric import stable_float
from bijux_phylogenetics.comparative.evolutionary_modes.search import (
    TransformedModeSearchResult as _TransformedModeSearchResult,
    best_pagel_lambda_fit as _best_pagel_lambda_fit,
    best_transformed_mode_fit as _best_transformed_mode_fit,
    continuous_boundary_assessment as _continuous_boundary_assessment,
    delta_identifiability_warnings_from_profile as _delta_identifiability_warnings_from_profile,
    early_burst_identifiability_warnings_from_profile as _early_burst_identifiability_warnings_from_profile,
    kappa_identifiability_warnings_from_profile as _kappa_identifiability_warnings_from_profile,
    lambda_identifiability_warnings_from_profile as _lambda_identifiability_warnings_from_profile,
    normalized_search_controls as _normalized_search_controls,
    ou_identifiability_warnings_from_profile as _ou_identifiability_warnings_from_profile,
    reject_nonparameterized_search_controls as _reject_nonparameterized_search_controls,
)
from bijux_phylogenetics.comparative.evolutionary_modes.tree_transforms import (
    clone_tree as _clone_tree,
    identity_covariance_matrix as _identity_covariance_matrix,
    rescale_tree_early_burst,
    rescale_tree_ornstein_uhlenbeck,
    rescale_tree_pagel_delta,
    rescale_tree_pagel_kappa,
    rescale_tree_pagel_lambda,
    rescale_tree_white_noise,
    transform_tree as _transform_tree,
    transform_tree_for_evolutionary_mode,
)
from bijux_phylogenetics.comparative.models import (
    ComparativeModelComparisonRow,
    _brownian_parameter_intervals,
    _build_residual_diagnostics,
    _comparison_row,
    _fit_intercept_only_model,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError
from bijux_phylogenetics.io.newick import dumps_newick


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
    """Fit a Brownian, white-noise, Pagel-lambda, Pagel-kappa, Pagel-delta, OU, or early-burst intercept-only trait model."""
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
    return _fit_evolutionary_mode_from_dataset(
        dataset,
        mode=mode,
        search_controls=search_controls,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def compare_continuous_evolutionary_modes(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    standard_error_trait: str | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    delta_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare the legacy Brownian, OU, and early-burst mode trio."""
    return _compare_selected_continuous_modes(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        modes=("brownian", "ornstein-uhlenbeck", "early-burst"),
        standard_error_trait=standard_error_trait,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def compare_fitcontinuous_model_ranking(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    modes: tuple[str, ...] | None = None,
    standard_error_trait: str | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 3.0),
    delta_bounds: tuple[float, float] = (0.0, 3.0),
    ou_bounds: tuple[float, float] = (0.0, 10.0),
    early_burst_bounds: tuple[float, float] = (0.0, 50.0),
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare the governed `fitContinuous` model set by AIC and AICc."""
    return _compare_selected_continuous_modes(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        modes=FITCONTINUOUS_MODEL_COMPARISON_ORDER if modes is None else modes,
        standard_error_trait=standard_error_trait,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        ou_bounds=ou_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def _compare_selected_continuous_modes(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None,
    modes: tuple[str, ...],
    standard_error_trait: str | None,
    lambda_bounds: tuple[float, float],
    kappa_bounds: tuple[float, float],
    delta_bounds: tuple[float, float],
    ou_bounds: tuple[float, float],
    early_burst_bounds: tuple[float, float],
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare a selected governed continuous-trait model set by information criteria."""
    if standard_error_trait is not None:
        raise ComparativeMethodError(
            "geiger::fitContinuous standard-error parity is explicitly excluded in this round; "
            f"received standard_error_trait='{standard_error_trait}' while Bijux still lacks "
            "measurement-error variance handling on the owned continuous-mode comparison surface"
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
    selected_modes = _comparison_modes(modes)
    fits: dict[str, ContinuousEvolutionaryModeFitReport] = {}
    rows: list[ComparativeModelComparisonRow] = []
    comparison_warnings: list[str] = []
    for mode in selected_modes:
        try:
            fit = _fit_evolutionary_mode_from_dataset(
                dataset,
                mode=mode,
                lambda_bounds=lambda_bounds,
                kappa_bounds=kappa_bounds,
                delta_bounds=delta_bounds,
                ou_bounds=ou_bounds,
                early_burst_bounds=early_burst_bounds,
            )
        except ComparativeMethodError as error:
            rows.append(
                ComparativeModelComparisonRow(
                    model=mode,
                    parameter_count=_mode_parameter_count(mode),
                    log_likelihood=math.nan,
                    aic=math.inf,
                    aicc=math.inf,
                    delta_aic=math.inf,
                    delta_aicc=math.inf,
                    rank=None,
                    comparable=False,
                    comparability_note=str(error),
                    selected=False,
                )
            )
            comparison_warnings.append(
                f"{mode} is not comparable on this dataset because the owned fit failed: {error}"
            )
            continue
        fits[mode] = fit
        row = _comparison_row(
            fit.mode,
            _mode_parameter_count(fit.mode),
            fit.log_likelihood,
            fit.taxon_count,
            likelihood_constant_policy=fit.likelihood_constant_policy,
        )
        if not math.isfinite(row.aicc):
            row.comparable = False
            row.comparability_note = (
                "sample size is too small to compute finite AICc for this parameter count"
            )
            comparison_warnings.append(
                f"{mode} is not comparable on AICc because the retained taxon count is too small for a {row.parameter_count}-parameter fit"
            )
        rows.append(row)
    if not fits:
        raise ComparativeMethodError(
            "no continuous evolutionary mode remained comparable for the requested dataset"
        )
    likelihood_constant_policy, noncomparable_likelihood_models = (
        rank_model_comparison_rows(
            rows,
            delta_threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        )
    )
    if noncomparable_likelihood_models:
        blocked_models = ", ".join(noncomparable_likelihood_models)
        comparison_warnings.append(
            "continuous-mode ranking excluded models with incompatible likelihood "
            f"constant policies: {blocked_models}"
        )
    selected_rows = [row for row in rows if row.selected]
    if not selected_rows:
        if noncomparable_likelihood_models:
            raise ComparativeMethodError(
                "mixed likelihood constant policies prevent ranking incompatible continuous-mode models"
            )
        raise ComparativeMethodError(
            "no finite AICc model remained available for continuous-mode comparison"
        )
    better_model = selected_rows[0].model
    likelihood_ratio_tests = []
    nested_modes = ("brownian", "ornstein-uhlenbeck", "early-burst")
    if all(mode in fits for mode in nested_modes):
        likelihood_ratio_tests = [
            _likelihood_ratio_test(
                comparison_id="brownian-vs-ornstein-uhlenbeck",
                left_fit=fits["brownian"],
                right_fit=fits["ornstein-uhlenbeck"],
            ),
            _likelihood_ratio_test(
                comparison_id="brownian-vs-early-burst",
                left_fit=fits["brownian"],
                right_fit=fits["early-burst"],
            ),
            _likelihood_ratio_test(
                comparison_id="ornstein-uhlenbeck-vs-early-burst",
                left_fit=fits["ornstein-uhlenbeck"],
                right_fit=fits["early-burst"],
            ),
        ]
    if any(mode not in nested_modes for mode in selected_modes):
        comparison_warnings.append(
            "likelihood-ratio tests remain reported only for the nested brownian/ornstein-uhlenbeck/early-burst subset; full fitcontinuous model ranking across white, lambda, kappa, and delta is governed by AIC and AICc"
        )
    if len(selected_rows) > 1:
        tied_models = ", ".join(row.model for row in selected_rows)
        comparison_warnings.append(
            f"multiple models remain tied at the selected AICc boundary: {tied_models}"
        )
    selected_model_boundary_assessment = (
        None if better_model not in fits else fits[better_model].boundary_assessment
    )
    stable_conclusion_supported = True
    if (
        selected_model_boundary_assessment is not None
        and selected_model_boundary_assessment.boundary_dominates_interpretation
    ):
        stable_conclusion_supported = False
        comparison_warnings.append(
            "selected continuous-mode fit remains boundary dominated on "
            f"{selected_model_boundary_assessment.affected_parameter}; stable "
            "conclusion support is withheld until the bounded parameter surface is "
            "reviewed directly"
        )
    return ContinuousEvolutionaryModeComparisonReport(
        tree_path=tree_path,
        traits_path=traits_path,
        trait=trait,
        taxon_count=dataset.readiness.tree_taxa,
        rows=rows,
        better_model=better_model,
        likelihood_constant_policy=likelihood_constant_policy,
        likelihood_comparison_policy=FITCONTINUOUS_MODEL_RANKING_POLICY,
        noncomparable_likelihood_models=noncomparable_likelihood_models,
        likelihood_ratio_tests=likelihood_ratio_tests,
        model_confidence_weight_basis=FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS,
        model_confidence_delta_threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        selected_model_akaike_weight=_selected_model_akaike_weight(rows),
        models_within_delta_aic_threshold=_models_within_delta_threshold(
            rows,
            criterion="aic",
            threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        models_within_delta_aicc_threshold=_models_within_delta_threshold(
            rows,
            criterion="aicc",
            threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        uncertainty_language=_model_confidence_uncertainty_language(
            rows,
            better_model=better_model,
            threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        warnings=comparison_warnings,
        selected_model_boundary_assessment=selected_model_boundary_assessment,
        stable_conclusion_supported=stable_conclusion_supported,
    )


def _comparison_modes(modes: tuple[str, ...] | None) -> tuple[str, ...]:
    if modes is None:
        return FITCONTINUOUS_MODEL_COMPARISON_ORDER
    unknown = [mode for mode in modes if mode not in ALLOWED_EVOLUTIONARY_MODES]
    if unknown:
        raise ComparativeMethodError(
            "unsupported comparison mode(s): " + ", ".join(sorted(unknown))
        )
    deduplicated: list[str] = []
    for mode in modes:
        if mode not in deduplicated:
            deduplicated.append(mode)
    if not deduplicated:
        raise ComparativeMethodError(
            "at least one supported continuous evolutionary mode is required for comparison"
        )
    return tuple(deduplicated)


def _mode_parameter_count(mode: str) -> int:
    if mode in {"brownian", "white-noise"}:
        return 2
    return 3


def _selected_model_akaike_weight(
    rows: list[ComparativeModelComparisonRow],
) -> float | None:
    selected_row = next((row for row in rows if row.selected), None)
    if selected_row is None:
        return None
    return selected_row.akaike_weight


def _models_within_delta_threshold(
    rows: list[ComparativeModelComparisonRow],
    *,
    criterion: str,
    threshold: float,
) -> list[str]:
    selected_rows: list[str] = []
    for row in rows:
        if not row.comparable:
            continue
        within_threshold = (
            row.within_delta_aic_threshold
            if criterion == "aic"
            else row.within_delta_aicc_threshold
        )
        if within_threshold:
            selected_rows.append(row.model)
    return selected_rows


def _model_confidence_uncertainty_language(
    rows: list[ComparativeModelComparisonRow],
    *,
    better_model: str,
    threshold: float,
) -> str:
    selected_row = next((row for row in rows if row.model == better_model), None)
    if selected_row is None or selected_row.akaike_weight is None:
        return (
            "model confidence is unresolved because no comparable candidate retained "
            "one finite AICc surface for Akaike-weight review"
        )
    nearby_models = [
        row.model
        for row in rows
        if row.comparable
        and row.model != better_model
        and row.within_delta_aicc_threshold is True
    ]
    if nearby_models:
        tied_models = ", ".join([better_model, *nearby_models])
        return (
            "model confidence is limited because "
            f"{tied_models} remain within {threshold:.1f} AICc units of the selected "
            f"surface; {better_model} carries Akaike weight "
            f"{selected_row.akaike_weight:.3f}"
        )
    if selected_row.akaike_weight < 0.9:
        return (
            f"model confidence is moderate because {better_model} carries Akaike "
            f"weight {selected_row.akaike_weight:.3f} even though no runner-up remains "
            f"within {threshold:.1f} AICc units"
        )
    return (
        f"model confidence is strong because {better_model} carries Akaike weight "
        f"{selected_row.akaike_weight:.3f} and no runner-up remains within "
        f"{threshold:.1f} AICc units"
    )


def _fit_evolutionary_mode_from_dataset(
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
        _reject_nonparameterized_search_controls(mode, search_controls)
    else:
        search_controls = _normalized_search_controls(search_controls)

    if mode == "brownian":
        transformed_tree = _clone_tree(dataset.tree)
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
        transformed_tree = _transform_tree(
            dataset.tree,
            mode="white-noise",
            parameter_value=1.0,
            sigsq=1.0,
        )
        covariance = stable_covariance(
            build_brownian_covariance_matrix(transformed_tree, dataset.taxa)
        )
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
        search_result = _best_pagel_lambda_fit(
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
        identifiability_warnings = _kappa_identifiability_warnings_from_profile(
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
        search_result = _best_transformed_mode_fit(
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
        identifiability_warnings = _delta_identifiability_warnings_from_profile(
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
        search_result = _best_transformed_mode_fit(
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
    boundary_assessment = _continuous_boundary_assessment(
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
        statistic=stable_float(statistic),
        degrees_of_freedom=1,
        p_value=stable_float(math.erfc(math.sqrt(statistic / 2.0))),
    )

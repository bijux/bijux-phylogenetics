from __future__ import annotations

from collections.abc import Callable
import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import ComparativeDataset
from bijux_phylogenetics.comparative.continuous.model_fitting import _comparison_row
from bijux_phylogenetics.comparative.evolutionary_modes.models import (
    ALLOWED_EVOLUTIONARY_MODES,
    FITCONTINUOUS_MODEL_COMPARISON_ORDER,
    FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
    FITCONTINUOUS_MODEL_CONFIDENCE_WEIGHT_BASIS,
    FITCONTINUOUS_MODEL_RANKING_POLICY,
    ContinuousEvolutionaryModeComparisonReport,
    ContinuousEvolutionaryModeFitReport,
    LikelihoodRatioTestResult,
)
from bijux_phylogenetics.comparative.evolutionary_modes.numeric import stable_float
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
    rank_model_comparison_rows,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

DatasetLoader = Callable[..., ComparativeDataset]
DatasetFitFunction = Callable[..., ContinuousEvolutionaryModeFitReport]


def compare_selected_continuous_modes(
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
    dataset_loader: DatasetLoader,
    fit_mode_from_dataset: DatasetFitFunction,
) -> ContinuousEvolutionaryModeComparisonReport:
    """Compare a selected governed continuous-trait model set by information criteria."""
    if standard_error_trait is not None:
        raise ComparativeMethodError(
            "geiger::fitContinuous standard-error parity is explicitly excluded in this round; "
            f"received standard_error_trait='{standard_error_trait}' while Bijux still lacks "
            "measurement-error variance handling on the owned continuous-mode comparison surface"
        )
    dataset = dataset_loader(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        minimum_taxa=3,
        require_rooted=True,
        require_binary=False,
    )
    selected_modes = comparison_modes(modes)
    fits: dict[str, ContinuousEvolutionaryModeFitReport] = {}
    rows: list[ComparativeModelComparisonRow] = []
    comparison_warnings: list[str] = []
    for mode in selected_modes:
        try:
            fit = fit_mode_from_dataset(
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
                    parameter_count=mode_parameter_count(mode),
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
            mode_parameter_count(fit.mode),
            fit.log_likelihood,
            fit.taxon_count,
            likelihood_constant_policy=fit.likelihood_constant_policy,
        )
        if not math.isfinite(row.aicc):
            row.comparable = False
            row.comparability_note = "sample size is too small to compute finite AICc for this parameter count"
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
    likelihood_ratio_tests: list[LikelihoodRatioTestResult] = []
    nested_modes = ("brownian", "ornstein-uhlenbeck", "early-burst")
    if all(mode in fits for mode in nested_modes):
        likelihood_ratio_tests = [
            likelihood_ratio_test(
                comparison_id="brownian-vs-ornstein-uhlenbeck",
                left_fit=fits["brownian"],
                right_fit=fits["ornstein-uhlenbeck"],
            ),
            likelihood_ratio_test(
                comparison_id="brownian-vs-early-burst",
                left_fit=fits["brownian"],
                right_fit=fits["early-burst"],
            ),
            likelihood_ratio_test(
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
        selected_model_akaike_weight=selected_model_akaike_weight(rows),
        models_within_delta_aic_threshold=models_within_delta_threshold(
            rows,
            criterion="aic",
            threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        models_within_delta_aicc_threshold=models_within_delta_threshold(
            rows,
            criterion="aicc",
            threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        uncertainty_language=model_confidence_uncertainty_language(
            rows,
            better_model=better_model,
            threshold=FITCONTINUOUS_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        warnings=comparison_warnings,
        selected_model_boundary_assessment=selected_model_boundary_assessment,
        stable_conclusion_supported=stable_conclusion_supported,
    )


def comparison_modes(modes: tuple[str, ...] | None) -> tuple[str, ...]:
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


def mode_parameter_count(mode: str) -> int:
    if mode in {"brownian", "white-noise"}:
        return 2
    return 3


def selected_model_akaike_weight(
    rows: list[ComparativeModelComparisonRow],
) -> float | None:
    selected_row = next((row for row in rows if row.selected), None)
    if selected_row is None:
        return None
    return selected_row.akaike_weight


def models_within_delta_threshold(
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


def model_confidence_uncertainty_language(
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


def likelihood_ratio_test(
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

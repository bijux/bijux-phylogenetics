from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralDiscreteDataset,
    load_discrete_dataset,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_state_order as _resolve_state_order,
)
from bijux_phylogenetics.comparative.model_selection import (
    ComparativeModelComparisonRow,
    rank_model_comparison_rows,
)
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

from .fitting import fit_discrete_mk_model_from_dataset
from .models import (
    DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY,
    DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
    DISCRETE_MK_MODEL_CONFIDENCE_WEIGHT_BASIS,
    DISCRETE_MK_MODEL_RANKING_POLICY,
    DiscreteMkFitReport,
    DiscreteMkModelComparisonReport,
)
from .transforms import (
    DISCRETE_DELTA_LOWER_BOUND,
    DISCRETE_DELTA_UPPER_BOUND,
    DISCRETE_EARLY_BURST_LOWER_BOUND,
    DISCRETE_EARLY_BURST_UPPER_BOUND,
    comparison_models,
    comparison_parameter_count,
)


def compare_discrete_mk_model_ranking(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    models: tuple[str, ...] | None = None,
    ascertainment_policy: str = "none",
    transform: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 1.0),
    delta_bounds: tuple[float, float] = (
        DISCRETE_DELTA_LOWER_BOUND,
        DISCRETE_DELTA_UPPER_BOUND,
    ),
    early_burst_bounds: tuple[float, float] = (
        DISCRETE_EARLY_BURST_LOWER_BOUND,
        DISCRETE_EARLY_BURST_UPPER_BOUND,
    ),
) -> DiscreteMkModelComparisonReport:
    """Compare governed discrete Mk models by AIC and AICc."""
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return compare_discrete_mk_model_ranking_from_dataset(
        dataset,
        models=models,
        ascertainment_policy=ascertainment_policy,
        transform=transform,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
        allowed_transition_pairs=allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )


def compare_discrete_mk_model_ranking_from_dataset(
    dataset: AncestralDiscreteDataset,
    *,
    models: tuple[str, ...] | None = None,
    ascertainment_policy: str = "none",
    transform: str | None = None,
    state_ordering: str = "unordered",
    ordered_states: list[str] | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    lambda_bounds: tuple[float, float] = (0.0, 1.0),
    kappa_bounds: tuple[float, float] = (0.0, 1.0),
    delta_bounds: tuple[float, float] = (
        DISCRETE_DELTA_LOWER_BOUND,
        DISCRETE_DELTA_UPPER_BOUND,
    ),
    early_burst_bounds: tuple[float, float] = (
        DISCRETE_EARLY_BURST_LOWER_BOUND,
        DISCRETE_EARLY_BURST_UPPER_BOUND,
    ),
) -> DiscreteMkModelComparisonReport:
    """Compare a selected governed discrete Mk model set by information criteria."""
    selected_models = comparison_models(models)
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    rows: list[ComparativeModelComparisonRow] = []
    fits: dict[str, DiscreteMkFitReport] = {}
    comparison_warnings: list[str] = []
    for model in selected_models:
        resolved_model = _resolve_discrete_model_name(model)
        try:
            fit = fit_discrete_mk_model_from_dataset(
                dataset,
                model=resolved_model,
                ascertainment_policy=ascertainment_policy,
                transform=transform,
                state_ordering=state_ordering,
                ordered_states=state_order,
                allowed_transition_pairs=allowed_transition_pairs,
                lambda_bounds=lambda_bounds,
                kappa_bounds=kappa_bounds,
                delta_bounds=delta_bounds,
                early_burst_bounds=early_burst_bounds,
            )
        except ComparativeMethodError as error:
            rows.append(
                ComparativeModelComparisonRow(
                    model=resolved_model,
                    parameter_count=comparison_parameter_count(
                        state_order=state_order,
                        model=resolved_model,
                        transform=transform,
                        state_ordering=state_ordering,
                        allowed_transition_pairs=allowed_transition_pairs,
                    ),
                    log_likelihood=math.nan,
                    aic=math.inf,
                    aicc=math.inf,
                    comparable=False,
                    comparability_note=str(error),
                    selected=False,
                    likelihood_constant_policy=DISCRETE_MK_LIKELIHOOD_CONSTANT_POLICY,
                )
            )
            comparison_warnings.append(
                f"{resolved_model} is not comparable on this dataset because the owned fit failed: {error}"
            )
            continue
        fits[resolved_model] = fit
        row = ComparativeModelComparisonRow(
            model=fit.model,
            parameter_count=fit.parameter_count,
            log_likelihood=fit.log_likelihood,
            aic=fit.aic,
            aicc=fit.aicc,
            likelihood_constant_policy=fit.likelihood_constant_policy,
        )
        if not math.isfinite(row.aicc):
            row.comparable = False
            row.comparability_note = "sample size is too small to compute finite AICc for this parameter count"
            comparison_warnings.append(
                f"{resolved_model} is not comparable on AICc because the retained taxon count is too small for a {row.parameter_count}-parameter fit"
            )
        rows.append(row)
    if not fits:
        raise ComparativeMethodError(
            "no discrete Mk model remained comparable for the requested dataset"
        )
    likelihood_constant_policy, noncomparable_likelihood_models = (
        rank_model_comparison_rows(
            rows,
            delta_threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        )
    )
    if noncomparable_likelihood_models:
        blocked_models = ", ".join(noncomparable_likelihood_models)
        comparison_warnings.append(
            "discrete Mk ranking excluded models with incompatible likelihood "
            f"constant policies: {blocked_models}"
        )
    selected_rows = [row for row in rows if row.selected]
    if not selected_rows:
        if noncomparable_likelihood_models:
            raise ComparativeMethodError(
                "mixed likelihood constant policies prevent ranking incompatible discrete Mk models"
            )
        raise ComparativeMethodError(
            "no finite AICc model remained available for discrete Mk comparison"
        )
    better_model = selected_rows[0].model
    if len(selected_rows) > 1:
        tied_models = ", ".join(row.model for row in selected_rows)
        comparison_warnings.append(
            f"multiple discrete Mk models remain tied at the selected AICc boundary: {tied_models}"
        )
    selected_fit = fits[better_model]
    if selected_fit.overparameterized:
        comparison_warnings.append(
            "selected discrete Mk fit remains overparameterized relative to the analyzed taxon count and should be interpreted cautiously"
        )
    if not selected_fit.optimizer_diagnostics.converged:
        comparison_warnings.append(
            "selected discrete Mk fit did not converge cleanly, so model-ranking confidence is reduced"
        )
    if (
        selected_fit.optimizer_diagnostics.hit_lower_parameter_bound
        or selected_fit.optimizer_diagnostics.hit_upper_parameter_bound
    ):
        comparison_warnings.append(
            "selected discrete Mk fit hits one or more optimizer bounds, so the winning rate surface should be treated as weakly identified"
        )
    return DiscreteMkModelComparisonReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        ascertainment_policy=ascertainment_policy,
        taxon_count=len(dataset.taxa),
        rows=rows,
        better_model=better_model,
        likelihood_constant_policy=likelihood_constant_policy,
        likelihood_comparison_policy=DISCRETE_MK_MODEL_RANKING_POLICY,
        noncomparable_likelihood_models=noncomparable_likelihood_models,
        model_confidence_weight_basis=DISCRETE_MK_MODEL_CONFIDENCE_WEIGHT_BASIS,
        model_confidence_delta_threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        selected_model_akaike_weight=_selected_model_akaike_weight(rows),
        models_within_delta_aic_threshold=_models_within_delta_threshold(
            rows,
            criterion="aic",
            threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        models_within_delta_aicc_threshold=_models_within_delta_threshold(
            rows,
            criterion="aicc",
            threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        uncertainty_language=_model_confidence_uncertainty_language(
            rows,
            better_model=better_model,
            threshold=DISCRETE_MK_MODEL_CONFIDENCE_DELTA_THRESHOLD,
        ),
        warnings=comparison_warnings,
    )


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
            "model confidence is unresolved because no comparable discrete Mk "
            "candidate retained one finite AICc surface for Akaike-weight review"
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

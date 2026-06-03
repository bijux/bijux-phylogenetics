from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.ancestral.common import (
    AncestralDiscreteDataset,
    load_discrete_dataset,
)
from bijux_phylogenetics.ancestral.discrete import (
    DiscreteModelBaselineComparison,
    DiscreteTransitionRateRow,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    invariant_pattern_log_probability as _invariant_pattern_log_probability,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    tree_log_likelihood as _tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    variable_pattern_log_probability as _variable_pattern_log_probability,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    build_transition_rate_rows as _build_transition_rate_rows,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_allowed_transition_pairs as _resolve_allowed_transition_pairs,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_discrete_model_name as _resolve_discrete_model_name,
)
from bijux_phylogenetics.ancestral.discrete.policy import (
    resolve_state_order as _resolve_state_order,
)
from bijux_phylogenetics.ancestral.discrete.reconstruction import (
    _detect_discrete_overparameterization,
)
from bijux_phylogenetics.comparative.common import tip_root_depths
from bijux_phylogenetics.comparative.model_selection import (
    compute_aic,
    compute_aicc,
)
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    summarize_ultrametric_tip_depths,
)

from .models import (
    DISCRETE_MK_LIKELIHOOD_COMPARISON_POLICY,
    DiscreteMkFitReport,
    DiscreteMkInputAudit,
    DiscreteMkPatternLikelihoodRow,
    DiscreteMkTransformBaselineComparison,
    resolve_discrete_mk_likelihood_constant_policy,
    validate_discrete_mk_ascertainment_policy,
)
from .transforms import (
    DISCRETE_DELTA_LOWER_BOUND,
    DISCRETE_DELTA_UPPER_BOUND,
    DISCRETE_EARLY_BURST_LOWER_BOUND,
    DISCRETE_EARLY_BURST_UPPER_BOUND,
    discrete_parameter_count,
    fit_discrete_mk_surface,
    resolve_discrete_transform_name,
    validate_discrete_transform_request,
)


def _normalize_transition_rate_rows(
    rows: list[DiscreteTransitionRateRow],
    *,
    state_ordering: str,
) -> list[DiscreteTransitionRateRow]:
    if state_ordering != "unordered":
        return rows
    return [
        DiscreteTransitionRateRow(
            source_state=row.source_state,
            target_state=row.target_state,
            transition_allowed=row.transition_allowed,
            step_distance=(1 if row.transition_allowed else row.step_distance),
            rate=row.rate,
        )
        for row in rows
    ]


def _build_pattern_likelihood_rows(
    dataset: AncestralDiscreteDataset,
    *,
    fit_tree,
    state_order: list[str],
    rate_matrix,
    root_prior,
    ascertainment_policy: str,
) -> list[DiscreteMkPatternLikelihoodRow]:
    raw_log_likelihood = _tree_log_likelihood(
        fit_tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode="observed",
        ascertainment_policy="none",
    )
    conditioning_log_probability = (
        _variable_pattern_log_probability(
            fit_tree,
            taxa=list(dataset.taxa),
            state_order=state_order,
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            root_prior_mode="observed",
        )
        if ascertainment_policy == "lewis-variable-only"
        else None
    )
    row_log_likelihood = raw_log_likelihood
    if conditioning_log_probability is not None:
        if conditioning_log_probability == float("-inf"):
            row_log_likelihood = float("-inf")
        else:
            row_log_likelihood -= conditioning_log_probability
    return [
        DiscreteMkPatternLikelihoodRow(
            pattern_id="pattern-1",
            pattern_weight=1,
            tip_states=tuple(dataset.states_by_taxon[taxon] for taxon in dataset.taxa),
            raw_log_likelihood=raw_log_likelihood,
            ascertainment_conditioning_log_probability=conditioning_log_probability,
            log_likelihood=row_log_likelihood,
        )
    ]


def _validate_pattern_likelihood_reconstruction(
    rows: list[DiscreteMkPatternLikelihoodRow],
    *,
    expected_total_log_likelihood: float,
) -> None:
    if len(rows) != 1:
        raise ValueError(
            "discrete Mk pattern likelihood reconstruction currently expects exactly one observed trait pattern"
        )
    reconstructed_total = math.fsum(
        row.pattern_weight * row.log_likelihood for row in rows
    )
    if not math.isclose(
        reconstructed_total,
        expected_total_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "discrete Mk pattern likelihood rows did not reconstruct the declared total likelihood"
        )


def fit_discrete_mk_model(
    tree_path: Path,
    traits_path: Path,
    *,
    trait: str,
    taxon_column: str | None = None,
    model: str = "equal-rates",
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
) -> DiscreteMkFitReport:
    """Fit one Mk discrete-trait model on a rooted tree."""
    dataset = load_discrete_dataset(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
    )
    return fit_discrete_mk_model_from_dataset(
        dataset,
        model=model,
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


def fit_discrete_mk_model_from_dataset(
    dataset: AncestralDiscreteDataset,
    *,
    model: str = "equal-rates",
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
) -> DiscreteMkFitReport:
    """Fit one Mk discrete-trait model from a native discrete dataset."""
    resolved_model = _resolve_discrete_model_name(model)
    resolved_ascertainment_policy = validate_discrete_mk_ascertainment_policy(
        ascertainment_policy
    )
    resolved_transform = resolve_discrete_transform_name(transform)
    state_order = _resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    validate_discrete_transform_request(
        dataset,
        transform=resolved_transform,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )
    resolved_allowed_transition_pairs = _resolve_allowed_transition_pairs(
        state_order,
        model=resolved_model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    (
        fit_tree,
        rate_matrix,
        root_prior,
        optimizer_diagnostics,
        transform_fit,
        transform_warning_rows,
    ) = fit_discrete_mk_surface(
        dataset,
        model=resolved_model,
        ascertainment_policy=resolved_ascertainment_policy,
        transform=resolved_transform,
        state_ordering=state_ordering,
        state_order=state_order,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
        lambda_bounds=lambda_bounds,
        kappa_bounds=kappa_bounds,
        delta_bounds=delta_bounds,
        early_burst_bounds=early_burst_bounds,
    )
    log_likelihood = _tree_log_likelihood(
        fit_tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode="observed",
        ascertainment_policy=resolved_ascertainment_policy,
    )
    ascertainment_conditioning_log_probability: float | None = None
    invariant_pattern_log_probability: float | None = None
    if resolved_ascertainment_policy == "lewis-variable-only":
        ascertainment_conditioning_log_probability = _variable_pattern_log_probability(
            fit_tree,
            taxa=list(dataset.taxa),
            state_order=state_order,
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            root_prior_mode="observed",
        )
        invariant_pattern_log_probability = _invariant_pattern_log_probability(
            fit_tree,
            taxa=list(dataset.taxa),
            state_order=state_order,
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            root_prior_mode="observed",
        )
    parameter_count = discrete_parameter_count(
        state_count=len(state_order),
        model=resolved_model,
        transform=resolved_transform,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    pattern_likelihood_rows = _build_pattern_likelihood_rows(
        dataset,
        fit_tree=fit_tree,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        ascertainment_policy=resolved_ascertainment_policy,
    )
    _validate_pattern_likelihood_reconstruction(
        pattern_likelihood_rows,
        expected_total_log_likelihood=log_likelihood,
    )
    aic = compute_aic(log_likelihood, parameter_count=parameter_count)
    aicc = compute_aicc(
        aic,
        sample_size=len(dataset.taxa),
        parameter_count=parameter_count,
    )
    baseline_comparison: DiscreteModelBaselineComparison | None = None
    transform_baseline_comparison: DiscreteMkTransformBaselineComparison | None = None
    if (
        resolved_model != "equal-rates"
        and allowed_transition_pairs is None
        and state_ordering == "unordered"
    ):
        baseline_fit = fit_discrete_mk_model_from_dataset(
            dataset,
            model="equal-rates",
            ascertainment_policy=resolved_ascertainment_policy,
            transform=resolved_transform,
            state_ordering=state_ordering,
            ordered_states=state_order,
            allowed_transition_pairs=None,
            lambda_bounds=lambda_bounds,
            kappa_bounds=kappa_bounds,
            delta_bounds=delta_bounds,
            early_burst_bounds=early_burst_bounds,
        )
        baseline_comparison = DiscreteModelBaselineComparison(
            baseline_model="equal-rates",
            baseline_log_likelihood=baseline_fit.log_likelihood,
            baseline_parameter_count=baseline_fit.parameter_count,
            baseline_aic=baseline_fit.aic,
            delta_log_likelihood=log_likelihood - baseline_fit.log_likelihood,
            delta_aic=aic - baseline_fit.aic,
            preferred_model_by_aic=(
                resolved_model if aic <= baseline_fit.aic else "equal-rates"
            ),
        )
    if resolved_transform is not None:
        transform_baseline_fit = fit_discrete_mk_model_from_dataset(
            dataset,
            model=resolved_model,
            ascertainment_policy=resolved_ascertainment_policy,
            transform=None,
            state_ordering=state_ordering,
            ordered_states=state_order,
            allowed_transition_pairs=allowed_transition_pairs,
            lambda_bounds=lambda_bounds,
            kappa_bounds=kappa_bounds,
            delta_bounds=delta_bounds,
            early_burst_bounds=early_burst_bounds,
        )
        transform_baseline_comparison = DiscreteMkTransformBaselineComparison(
            baseline_transform="untransformed",
            baseline_log_likelihood=transform_baseline_fit.log_likelihood,
            baseline_parameter_count=transform_baseline_fit.parameter_count,
            baseline_aic=transform_baseline_fit.aic,
            delta_log_likelihood=log_likelihood - transform_baseline_fit.log_likelihood,
            delta_aic=aic - transform_baseline_fit.aic,
            preferred_transform_by_aic=(
                resolved_transform
                if aic <= transform_baseline_fit.aic
                else "untransformed"
            ),
        )
    overparameterized = _detect_discrete_overparameterization(
        taxon_count=len(dataset.taxa),
        parameter_count=parameter_count,
    )
    warnings = list(dataset.warnings)
    if overparameterized:
        warnings.append(
            "the discrete Mk likelihood fit is likely overparameterized relative to the analyzed taxon count"
        )
    for warning in transform_warning_rows:
        warnings.append(warning.message)
    if not optimizer_diagnostics.converged:
        warnings.append(
            "the discrete Mk optimizer did not converge and should be interpreted cautiously"
        )
    if (
        optimizer_diagnostics.hit_lower_parameter_bound
        or optimizer_diagnostics.hit_upper_parameter_bound
    ):
        warnings.append(
            "one or more discrete Mk rate parameters hit an optimizer bound and should be interpreted as weakly identified"
        )
    if (
        baseline_comparison is not None
        and baseline_comparison.preferred_model_by_aic == "equal-rates"
    ):
        warnings.append(
            "the equal-rates baseline remains preferred by AIC over the requested discrete Mk model"
        )
    if (
        transform_baseline_comparison is not None
        and transform_baseline_comparison.preferred_transform_by_aic == "untransformed"
    ):
        warnings.append(
            "the untransformed branch-length baseline remains preferred by AIC over the requested discrete Mk transform"
        )
    input_audit = build_discrete_mk_input_audit(dataset, warnings=warnings)
    return DiscreteMkFitReport(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        taxon_column=dataset.taxon_column,
        trait=dataset.trait,
        model=resolved_model,
        ascertainment_policy=resolved_ascertainment_policy,
        state_ordering=state_ordering,
        state_order=state_order,
        taxon_count=len(dataset.taxa),
        input_audit=input_audit,
        log_likelihood=log_likelihood,
        ascertainment_conditioning_log_probability=(
            ascertainment_conditioning_log_probability
        ),
        invariant_pattern_log_probability=invariant_pattern_log_probability,
        parameter_count=parameter_count,
        aic=aic,
        aicc=aicc,
        likelihood_constant_policy=resolve_discrete_mk_likelihood_constant_policy(
            resolved_ascertainment_policy
        ),
        likelihood_comparison_policy=DISCRETE_MK_LIKELIHOOD_COMPARISON_POLICY,
        pattern_likelihood_rows=pattern_likelihood_rows,
        transition_rate_rows=_normalize_transition_rate_rows(
            _build_transition_rate_rows(
                state_order=state_order,
                state_ordering=state_ordering,
                rate_matrix=rate_matrix,
                allowed_transition_pairs=resolved_allowed_transition_pairs,
            ),
            state_ordering=state_ordering,
        ),
        allowed_transition_pairs=[
            (state_order[left_index], state_order[right_index])
            for left_index, right_index in sorted(resolved_allowed_transition_pairs)
        ],
        optimizer_diagnostics=optimizer_diagnostics,
        overparameterized=overparameterized,
        transform_fit=transform_fit,
        transform_baseline_comparison=transform_baseline_comparison,
        baseline_comparison=baseline_comparison,
    )


def build_discrete_mk_input_audit(
    dataset: AncestralDiscreteDataset,
    *,
    warnings: list[str],
) -> DiscreteMkInputAudit:
    ultrametric_summary = summarize_ultrametric_tip_depths(
        tip_root_depths(dataset.tree, dataset.taxa),
        tolerance=1e-12,
    )
    return DiscreteMkInputAudit(
        tree_path=dataset.tree_path,
        traits_path=dataset.traits_path,
        trait=dataset.trait,
        taxon_count=len(dataset.taxa),
        taxa=list(dataset.taxa),
        observed_states=list(dataset.observed_states),
        state_counts=dict(dataset.state_counts),
        sparse_states=list(dataset.sparse_states),
        tree_is_ultrametric=ultrametric_summary.ultrametric,
        minimum_root_to_tip_depth=ultrametric_summary.minimum_tip_depth,
        maximum_root_to_tip_depth=ultrametric_summary.maximum_tip_depth,
        ultrametric_policy="accept-rooted-trees-and-report-ultrametricity",
        missing_value_policy="prune-overlapping-missing-values",
        missing_from_traits=list(dataset.alignment_report.dropped_tree_taxa),
        extra_trait_taxa=list(dataset.alignment_report.dropped_trait_taxa),
        pruned_missing_value_taxa=list(dataset.dropped_missing_taxa),
        warnings=warnings,
    )

from __future__ import annotations

import numpy

from bijux_phylogenetics.ancestral.common import node_signature

from ..models import (
    DiscreteAncestralEstimate,
    DiscreteLikelihoodFitResult,
    DiscreteOptimizerDiagnostics,
)
from ..policy import (
    parameter_count,
    reported_discrete_log_likelihood,
    rerooting_method_compatibility,
    resolve_allowed_transition_pairs,
    resolve_root_prior,
    resolve_state_order,
    uniform_root_prior,
)
from .likelihood_math import tree_log_likelihood
from .optimizer_search import (
    DISCRETE_LOG_PARAMETER_LOWER_BOUND,
    DISCRETE_LOG_PARAMETER_UPPER_BOUND,
    _DiscreteOptimizationRun,
    build_discrete_initial_candidates,
    optimize_log_parameters,
    optimize_single_log_parameter,
    regularize_plateau_log_parameters,
)
from .posterior_probabilities import estimate_marginal_state_probabilities
from .rate_matrix import build_transition_rate_rows, rate_matrix_from_log_parameters


def reconstruct_likelihood_estimates(
    dataset,
    *,
    model: str,
    state_ordering: str,
    ordered_states: list[str] | None,
    root_prior_mode: str = "equal",
    fixed_root_state: str | None = None,
    allowed_transition_pairs: list[tuple[str, str]] | None = None,
    include_equal_rates_baseline: bool = True,
    stable_probability_mapping,
    material_state_set,
    select_most_likely_state,
    build_discrete_estimate,
    detect_discrete_overparameterization,
) -> DiscreteLikelihoodFitResult:
    state_order = resolve_state_order(
        dataset.observed_states,
        state_ordering=state_ordering,
        ordered_states=ordered_states,
    )
    resolved_allowed_transition_pairs = resolve_allowed_transition_pairs(
        state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    rate_matrix, default_root_prior, optimizer_diagnostics = fit_discrete_mk_model(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    root_prior = resolve_root_prior(
        state_order,
        state_counts=dataset.state_counts,
        mode=root_prior_mode,
        fixed_root_state=fixed_root_state,
        default_root_prior=default_root_prior,
    )
    log_likelihood = tree_log_likelihood(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    reported_log_likelihood = reported_discrete_log_likelihood(
        log_likelihood,
        root_prior_mode=root_prior_mode,
        state_count=len(state_order),
    )
    posterior_by_node = estimate_marginal_state_probabilities(
        dataset.tree,
        dataset.states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    estimates: list[DiscreteAncestralEstimate] = []
    for node in dataset.tree.iter_nodes():
        signature = node_signature(node)
        raw_probabilities = {
            state: float(probability)
            for state, probability in posterior_by_node[signature].items()
        }
        probabilities = stable_probability_mapping(raw_probabilities)
        material_states = material_state_set(
            probabilities,
            preferred_order=state_order,
        )
        estimates.append(
            build_discrete_estimate(
                node=signature,
                node_name=node.name,
                is_tip=node.is_leaf(),
                descendant_taxa=[
                    descendant.name
                    for descendant in node.iter_leaves()
                    if descendant.name is not None
                ],
                state_set=material_states,
                most_likely_state=select_most_likely_state(
                    raw_probabilities,
                    preferred_order=state_order,
                ),
                state_probabilities=probabilities,
            )
        )
    count = parameter_count(
        len(state_order),
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=resolved_allowed_transition_pairs,
    )
    aic = (2.0 * count) - (2.0 * reported_log_likelihood)
    baseline_comparison = None
    if (
        include_equal_rates_baseline
        and model != "equal-rates"
        and allowed_transition_pairs is None
    ):
        try:
            baseline_fit = reconstruct_likelihood_estimates(
                dataset,
                model="equal-rates",
                state_ordering=state_ordering,
                ordered_states=state_order,
                root_prior_mode=root_prior_mode,
                fixed_root_state=fixed_root_state,
                allowed_transition_pairs=None,
                include_equal_rates_baseline=False,
                stable_probability_mapping=stable_probability_mapping,
                material_state_set=material_state_set,
                select_most_likely_state=select_most_likely_state,
                build_discrete_estimate=build_discrete_estimate,
                detect_discrete_overparameterization=detect_discrete_overparameterization,
            )
        except (RuntimeError, ValueError):
            baseline_fit = None
        if baseline_fit is not None:
            from ..models import DiscreteModelBaselineComparison

            baseline_comparison = DiscreteModelBaselineComparison(
                baseline_model="equal-rates",
                baseline_log_likelihood=baseline_fit.log_likelihood,
                baseline_parameter_count=baseline_fit.parameter_count,
                baseline_aic=baseline_fit.aic,
                delta_log_likelihood=(
                    reported_log_likelihood - baseline_fit.log_likelihood
                ),
                delta_aic=aic - baseline_fit.aic,
                preferred_model_by_aic=(
                    model if aic <= baseline_fit.aic else "equal-rates"
                ),
            )
    overparameterized = detect_discrete_overparameterization(
        taxon_count=len(dataset.taxa),
        parameter_count=count,
    )
    return DiscreteLikelihoodFitResult(
        estimates=estimates,
        ordered_states=(state_order if state_ordering == "ordered" else []),
        state_order=state_order,
        rerooting_method_compatibility=rerooting_method_compatibility(
            model=model,
            state_ordering=state_ordering,
            root_prior_mode=root_prior_mode,
        ),
        log_likelihood=reported_log_likelihood,
        parameter_count=count,
        aic=aic,
        transition_rate_rows=build_transition_rate_rows(
            state_order=state_order,
            state_ordering=state_ordering,
            rate_matrix=rate_matrix,
            allowed_transition_pairs=resolved_allowed_transition_pairs,
        ),
        allowed_transition_pairs=[
            (state_order[left_index], state_order[right_index])
            for left_index, right_index in sorted(resolved_allowed_transition_pairs)
        ],
        optimizer_diagnostics=optimizer_diagnostics,
        overparameterized=overparameterized,
        baseline_comparison=baseline_comparison,
    )


def fit_discrete_mk_model(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str = "equal",
    ascertainment_policy: str = "none",
) -> tuple[numpy.ndarray, numpy.ndarray, DiscreteOptimizerDiagnostics]:
    count = parameter_count(
        len(state_order),
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    if count < 1:
        raise ValueError(
            "discrete ancestral reconstruction requires at least one allowed transition"
        )
    if count == 1:
        best_run = optimize_single_log_parameter(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            ascertainment_policy=ascertainment_policy,
        )
        best_log_parameters = best_run.log_parameters
    else:
        initial_candidates = build_discrete_initial_candidates(
            parameter_count=count,
            model=model,
        )
        best_log_parameters: numpy.ndarray | None = None
        best_log_likelihood = float("-inf")
        best_run: _DiscreteOptimizationRun | None = None
        for initial_scale, initial in initial_candidates:
            run = optimize_log_parameters(
                tree,
                states_by_taxon,
                state_order=state_order,
                model=model,
                state_ordering=state_ordering,
                allowed_transition_pairs=allowed_transition_pairs,
                root_prior_mode=root_prior_mode,
                ascertainment_policy=ascertainment_policy,
                initial_log_parameters=initial,
                initial_scale=initial_scale,
            )
            if run.log_likelihood > best_log_likelihood:
                best_run = run
                best_log_parameters = run.log_parameters
                best_log_likelihood = run.log_likelihood
        if best_log_parameters is None or best_run is None:
            raise RuntimeError(
                "discrete ancestral optimization did not produce rate parameters"
            )
        best_log_parameters = regularize_plateau_log_parameters(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            ascertainment_policy=ascertainment_policy,
            log_parameters=best_log_parameters,
            reference_log_likelihood=best_log_likelihood,
        )
    rate_matrix = rate_matrix_from_log_parameters(
        best_log_parameters,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    root_prior = uniform_root_prior(len(state_order))
    return (
        rate_matrix,
        root_prior,
        DiscreteOptimizerDiagnostics(
            optimizer_name=best_run.optimizer_name,
            parameter_count=count,
            initial_candidate_count=best_run.initial_candidate_count,
            best_initial_scale=best_run.initial_scale,
            converged=best_run.converged,
            iteration_count=best_run.iteration_count,
            function_evaluation_count=best_run.function_evaluation_count,
            simplex_shrink_count=best_run.simplex_shrink_count,
            hit_lower_parameter_bound=bool(
                numpy.any(
                    best_log_parameters <= (DISCRETE_LOG_PARAMETER_LOWER_BOUND + 1e-9)
                )
            ),
            hit_upper_parameter_bound=bool(
                numpy.any(
                    best_log_parameters >= (DISCRETE_LOG_PARAMETER_UPPER_BOUND - 1e-9)
                )
            ),
        ),
    )

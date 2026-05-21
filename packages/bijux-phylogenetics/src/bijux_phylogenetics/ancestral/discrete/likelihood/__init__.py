from __future__ import annotations

from dataclasses import dataclass
import math
import random

import numpy

from bijux_phylogenetics.ancestral.common import node_signature

from ..models import (
    DiscreteAncestralEstimate,
    DiscreteLikelihoodFitResult,
    DiscreteOptimizerDiagnostics,
    DiscreteTransitionRateRow,
)
from ..policy import (
    normalize_array,
    parameter_count,
    reported_discrete_log_likelihood,
    rerooting_method_compatibility,
    resolve_allowed_transition_pairs,
    resolve_root_prior,
    resolve_state_order,
    transition_allowed,
    uniform_root_prior,
)
from .rate_matrix import build_transition_rate_rows
from .rate_matrix import rate_matrix_from_log_parameters

DISCRETE_LOG_PARAMETER_LOWER_BOUND = -10.0
DISCRETE_LOG_PARAMETER_UPPER_BOUND = 10.0
DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE = 1e-10


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
        )
        best_log_parameters = best_run.log_parameters
        best_log_likelihood = best_run.log_likelihood
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


def build_discrete_initial_candidates(
    *,
    parameter_count: int,
    model: str,
) -> list[tuple[float, numpy.ndarray]]:
    uniform_scales = (0.1, 1.0, 3.0)
    candidates = [
        (
            scale,
            numpy.full(parameter_count, math.log(scale), dtype=float),
        )
        for scale in uniform_scales
    ]
    if model != "all-rates-different" or parameter_count <= 2:
        return candidates
    exploratory_rng = random.Random(17)  # nosec B311
    for _ in range(32):
        log_parameters = numpy.array(
            [
                exploratory_rng.uniform(
                    DISCRETE_LOG_PARAMETER_LOWER_BOUND + 1.0,
                    DISCRETE_LOG_PARAMETER_UPPER_BOUND - 1.0,
                )
                for _ in range(parameter_count)
            ],
            dtype=float,
        )
        candidates.append(
            (float(math.exp(float(numpy.mean(log_parameters)))), log_parameters)
        )
    return candidates


@dataclass(slots=True)
class _DiscreteOptimizationRun:
    log_parameters: numpy.ndarray
    log_likelihood: float
    optimizer_name: str
    initial_candidate_count: int
    initial_scale: float
    converged: bool
    iteration_count: int
    function_evaluation_count: int
    simplex_shrink_count: int


def optimize_log_parameters(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
    initial_log_parameters: numpy.ndarray,
    initial_scale: float,
) -> _DiscreteOptimizationRun:
    simplex = [
        numpy.clip(
            initial_log_parameters.copy(),
            DISCRETE_LOG_PARAMETER_LOWER_BOUND,
            DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        )
    ]
    for index in range(initial_log_parameters.size):
        vertex = simplex[0].copy()
        vertex[index] += 0.75
        simplex.append(
            numpy.clip(
                vertex,
                DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
        )
    scores = [
        evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            log_parameters=vertex,
        )
        for vertex in simplex
    ]
    function_evaluation_count = len(scores)
    alpha = 1.0
    gamma = 2.0
    rho = 0.5
    sigma = 0.5
    converged = False
    iteration_count = 0
    simplex_shrink_count = 0
    for _iteration_count in range(1, 601):
        iteration_count = _iteration_count
        ordering = sorted(
            range(len(simplex)),
            key=lambda index: scores[index],
            reverse=True,
        )
        simplex = [simplex[index] for index in ordering]
        scores = [scores[index] for index in ordering]
        if (
            max(numpy.linalg.norm(vertex - simplex[0]) for vertex in simplex[1:]) < 1e-6
            and max(abs(score - scores[0]) for score in scores[1:]) < 1e-9
        ):
            converged = True
            break
        centroid = numpy.mean(simplex[:-1], axis=0)
        reflected = numpy.clip(
            centroid + alpha * (centroid - simplex[-1]),
            DISCRETE_LOG_PARAMETER_LOWER_BOUND,
            DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        )
        reflected_score = evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            log_parameters=reflected,
        )
        function_evaluation_count += 1
        if scores[0] >= reflected_score > scores[-2]:
            simplex[-1] = reflected
            scores[-1] = reflected_score
            continue
        if reflected_score > scores[0]:
            expanded = numpy.clip(
                centroid + gamma * (reflected - centroid),
                DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
            expanded_score = evaluate_log_likelihood(
                tree,
                states_by_taxon,
                state_order=state_order,
                model=model,
                state_ordering=state_ordering,
                allowed_transition_pairs=allowed_transition_pairs,
                root_prior_mode=root_prior_mode,
                log_parameters=expanded,
            )
            function_evaluation_count += 1
            if expanded_score > reflected_score:
                simplex[-1] = expanded
                scores[-1] = expanded_score
            else:
                simplex[-1] = reflected
                scores[-1] = reflected_score
            continue
        if reflected_score > scores[-1]:
            contracted = numpy.clip(
                centroid + rho * (reflected - centroid),
                DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
        else:
            contracted = numpy.clip(
                centroid + rho * (simplex[-1] - centroid),
                DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
        contracted_score = evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            log_parameters=contracted,
        )
        function_evaluation_count += 1
        if contracted_score > scores[-1]:
            simplex[-1] = contracted
            scores[-1] = contracted_score
            continue
        best_vertex = simplex[0]
        new_simplex = [best_vertex]
        new_scores = [scores[0]]
        simplex_shrink_count += 1
        for vertex in simplex[1:]:
            shrunk = numpy.clip(
                best_vertex + sigma * (vertex - best_vertex),
                DISCRETE_LOG_PARAMETER_LOWER_BOUND,
                DISCRETE_LOG_PARAMETER_UPPER_BOUND,
            )
            new_simplex.append(shrunk)
            new_scores.append(
                evaluate_log_likelihood(
                    tree,
                    states_by_taxon,
                    state_order=state_order,
                    model=model,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                    root_prior_mode=root_prior_mode,
                    log_parameters=shrunk,
                )
            )
            function_evaluation_count += 1
        simplex = new_simplex
        scores = new_scores
    best_index = max(range(len(scores)), key=lambda index: scores[index])
    return _DiscreteOptimizationRun(
        log_parameters=simplex[best_index],
        log_likelihood=scores[best_index],
        optimizer_name="nelder-mead",
        initial_candidate_count=len(simplex),
        initial_scale=initial_scale,
        converged=converged,
        iteration_count=iteration_count,
        function_evaluation_count=function_evaluation_count,
        simplex_shrink_count=simplex_shrink_count,
    )


def optimize_single_log_parameter(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
) -> _DiscreteOptimizationRun:
    from bijux_phylogenetics.comparative.search import (
        run_bounded_golden_section_maximization,
    )

    search_result = run_bounded_golden_section_maximization(
        lower_bound=DISCRETE_LOG_PARAMETER_LOWER_BOUND,
        upper_bound=DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        evaluate=lambda parameter: evaluate_single_log_parameter_candidate(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            parameter=parameter,
        ),
        optimizer_name="golden-section-search",
    )
    best_parameter = search_result.parameter_value
    return _DiscreteOptimizationRun(
        log_parameters=numpy.array([best_parameter], dtype=float),
        log_likelihood=search_result.objective_value,
        optimizer_name=search_result.diagnostics.optimizer_name,
        initial_candidate_count=1,
        initial_scale=float(math.exp(best_parameter)),
        converged=search_result.diagnostics.converged,
        iteration_count=max(search_result.diagnostics.function_evaluation_count - 2, 0),
        function_evaluation_count=search_result.diagnostics.function_evaluation_count,
        simplex_shrink_count=0,
    )


def evaluate_single_log_parameter_candidate(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
    parameter: float,
) -> tuple[float, float]:
    objective_value = evaluate_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode=root_prior_mode,
        log_parameters=numpy.array([parameter], dtype=float),
    )
    return parameter, objective_value


def evaluate_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
    log_parameters: numpy.ndarray,
) -> float:
    rate_matrix = rate_matrix_from_log_parameters(
        log_parameters,
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    return tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=None
        if root_prior_mode == "observed"
        else uniform_root_prior(len(state_order)),
        root_prior_mode=root_prior_mode,
    )


def regularize_plateau_log_parameters(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
    log_parameters: numpy.ndarray,
    reference_log_likelihood: float,
) -> numpy.ndarray:
    regularized = log_parameters.copy()
    for parameter_index in range(regularized.size):
        current_value = float(regularized[parameter_index])
        if current_value <= (DISCRETE_LOG_PARAMETER_LOWER_BOUND + 1e-9):
            continue
        current_reference_log_likelihood = evaluate_log_likelihood(
            tree,
            states_by_taxon,
            state_order=state_order,
            model=model,
            state_ordering=state_ordering,
            allowed_transition_pairs=allowed_transition_pairs,
            root_prior_mode=root_prior_mode,
            log_parameters=regularized,
        )
        if current_reference_log_likelihood < (
            reference_log_likelihood - DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE
        ):
            current_reference_log_likelihood = reference_log_likelihood
        low = DISCRETE_LOG_PARAMETER_LOWER_BOUND
        high = current_value
        for _ in range(80):
            midpoint = (low + high) / 2.0
            candidate = regularized.copy()
            candidate[parameter_index] = midpoint
            candidate_log_likelihood = evaluate_log_likelihood(
                tree,
                states_by_taxon,
                state_order=state_order,
                model=model,
                state_ordering=state_ordering,
                allowed_transition_pairs=allowed_transition_pairs,
                root_prior_mode=root_prior_mode,
                log_parameters=candidate,
            )
            if (
                current_reference_log_likelihood - candidate_log_likelihood
                <= DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE
            ):
                high = midpoint
            else:
                low = midpoint
        regularized[parameter_index] = high
    return regularized


def tree_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
) -> float:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def visit(node) -> tuple[numpy.ndarray, float]:
        if node.is_leaf():
            likelihood = numpy.zeros(len(state_order), dtype=float)
            likelihood[state_index[states_by_taxon[node.name]]] = 1.0
            return likelihood, 0.0
        partial = numpy.ones(len(state_order), dtype=float)
        log_scale = 0.0
        for child in node.children:
            child_partial, child_scale = visit(child)
            branch_transition = transition(branch_length(child))
            partial *= branch_transition @ child_partial
            log_scale += child_scale
        scale = float(partial.sum())
        if scale <= 0.0:
            return partial, float("-inf")
        partial /= scale
        return partial, log_scale + math.log(scale)

    root_partial, subtree_log_scale = visit(tree.root)
    if root_prior_mode == "observed":
        observed_root_total = float(root_partial.sum())
        if observed_root_total <= 0.0:
            return float("-inf")
        root_scale = float((root_partial @ root_partial) / observed_root_total)
        if root_scale <= 0.0:
            return float("-inf")
        return subtree_log_scale + math.log(root_scale)
    if root_prior is None:
        raise ValueError("root_prior is required unless root_prior_mode is 'observed'")
    root_weight = root_prior * root_partial
    root_scale = float(root_weight.sum())
    if root_scale <= 0.0:
        return float("-inf")
    return subtree_log_scale + math.log(root_scale)


def estimate_marginal_state_probabilities(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    state_index = {state: index for index, state in enumerate(state_order)}
    posterior_by_node: dict[str, numpy.ndarray] = {}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(branch_length)
        if cached is None:
            cached = transition_probability_matrix(rate_matrix, branch_length)
            transition_cache[branch_length] = cached
        return cached

    def postorder(node) -> numpy.ndarray:
        signature = node_signature(node)
        if node.is_leaf():
            partial = numpy.zeros(len(state_order), dtype=float)
            partial[state_index[states_by_taxon[node.name]]] = 1.0
            posterior_by_node[signature] = partial
            return partial
        partial = numpy.ones(len(state_order), dtype=float)
        for child in node.children:
            child_partial = postorder(child)
            partial *= transition(branch_length(child)) @ child_partial
        partial = normalize_array(partial)
        posterior_by_node[signature] = partial
        return partial

    postorder(tree.root)
    root_signature = node_signature(tree.root)
    posterior_by_node[root_signature] = normalize_array(
        root_prior * posterior_by_node[root_signature]
    )

    def preorder(node) -> None:
        parent_signature = node_signature(node)
        if node.is_leaf():
            return
        parent_probabilities = posterior_by_node[parent_signature]
        for child in node.children:
            if child.is_leaf():
                continue
            child_signature = node_signature(child)
            branch_transition = transition(branch_length(child))
            child_probabilities = posterior_by_node[child_signature]
            denominator = child_probabilities @ branch_transition
            updated = (parent_probabilities / denominator) @ branch_transition
            posterior_by_node[child_signature] = normalize_array(
                updated * child_probabilities
            )
            preorder(child)

    preorder(tree.root)
    return {
        node: {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                normalize_array(probabilities),
                strict=True,
            )
        }
        for node, probabilities in posterior_by_node.items()
    }


def transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    if branch_length <= 0.0:
        return numpy.eye(rate_matrix.shape[0], dtype=float)
    eigenvalues, eigenvectors = numpy.linalg.eig(rate_matrix)
    inverse_vectors = numpy.linalg.inv(eigenvectors)
    diagonal = numpy.diag(numpy.exp(eigenvalues * branch_length))
    transition = eigenvectors @ diagonal @ inverse_vectors
    transition = numpy.real_if_close(transition, tol=1000).astype(float)
    transition[transition < 0.0] = 0.0
    row_sums = transition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return transition / row_sums


def branch_length(node) -> float:
    if node.branch_length is None:
        return 1.0
    return max(float(node.branch_length), 0.0)

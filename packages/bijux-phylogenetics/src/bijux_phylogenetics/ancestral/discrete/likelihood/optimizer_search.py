from __future__ import annotations

from dataclasses import dataclass
import math
import random

import numpy

from ..policy import uniform_root_prior
from .likelihood_math import tree_log_likelihood
from .rate_matrix import rate_matrix_from_log_parameters

DISCRETE_LOG_PARAMETER_LOWER_BOUND = -10.0
DISCRETE_LOG_PARAMETER_UPPER_BOUND = 10.0
DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE = 1e-10


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
    ascertainment_policy: str,
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
            ascertainment_policy=ascertainment_policy,
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
            ascertainment_policy=ascertainment_policy,
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
                ascertainment_policy=ascertainment_policy,
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
            ascertainment_policy=ascertainment_policy,
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
                    ascertainment_policy=ascertainment_policy,
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
    ascertainment_policy: str,
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
            ascertainment_policy=ascertainment_policy,
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
    ascertainment_policy: str,
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
        ascertainment_policy=ascertainment_policy,
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
    ascertainment_policy: str,
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
        ascertainment_policy=ascertainment_policy,
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
    ascertainment_policy: str,
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
            ascertainment_policy=ascertainment_policy,
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
                ascertainment_policy=ascertainment_policy,
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

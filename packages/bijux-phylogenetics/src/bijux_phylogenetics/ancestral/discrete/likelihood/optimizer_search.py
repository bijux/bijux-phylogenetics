from __future__ import annotations

from dataclasses import dataclass
import math

import numpy

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from ..policy import uniform_root_prior
from .likelihood_math import tree_log_likelihood
from .parameter_layout import (
    broadcast_discrete_mk_log_parameters,
    build_discrete_mk_parameter_layout,
)
from .rate_matrix import rate_matrix_from_log_parameters

DISCRETE_LOG_PARAMETER_LOWER_BOUND = -10.0
DISCRETE_LOG_PARAMETER_UPPER_BOUND = 10.0
DISCRETE_PLATEAU_REGULARIZATION_TOLERANCE = 1e-10


def build_discrete_initial_candidates(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    parameter_count: int,
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
    ascertainment_policy: str,
) -> list[tuple[float, numpy.ndarray]]:
    uniform_scales = (0.1, 1.0, 3.0)
    candidates = _deduplicated_discrete_initial_candidates(
        [
            (
                scale,
                numpy.full(parameter_count, math.log(scale), dtype=float),
            )
            for scale in uniform_scales
        ]
    )
    if parameter_count <= 1:
        return candidates
    equal_layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model="equal-rates",
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    target_layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    equal_run = optimize_single_log_parameter(
        tree,
        states_by_taxon,
        state_order=state_order,
        model="equal-rates",
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode=root_prior_mode,
        ascertainment_policy=ascertainment_policy,
    )
    equal_seed = numpy.array(
        broadcast_discrete_mk_log_parameters(
            equal_layout,
            target_layout,
            tuple(float(value) for value in equal_run.log_parameters),
        ),
        dtype=float,
    )
    candidates.append(
        (
            _scale_from_log_parameters(equal_seed),
            equal_seed,
        )
    )
    if model == "symmetric":
        candidates.extend(_perturbed_discrete_seed_candidates(equal_seed, limit=4))
        return _deduplicated_discrete_initial_candidates(candidates)
    symmetric_seed = _fit_symmetric_warm_start(
        tree,
        states_by_taxon,
        state_order=state_order,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode=root_prior_mode,
        ascertainment_policy=ascertainment_policy,
        equal_seed=equal_seed,
        target_layout=target_layout,
    )
    if symmetric_seed is not None:
        candidates.append(
            (
                _scale_from_log_parameters(symmetric_seed),
                symmetric_seed,
            )
        )
        candidates.extend(_perturbed_discrete_seed_candidates(symmetric_seed, limit=2))
    candidates.append(
        (
            _scale_from_log_parameters(equal_seed),
            _state_count_directional_seed(
                equal_seed,
                state_order=state_order,
                states_by_taxon=states_by_taxon,
                layout=target_layout,
            ),
        )
    )
    candidates.extend(_perturbed_discrete_seed_candidates(equal_seed, limit=2))
    candidates = _deduplicated_discrete_initial_candidates(candidates)
    return candidates


def _scale_from_log_parameters(log_parameters: numpy.ndarray) -> float:
    return float(math.exp(float(numpy.mean(log_parameters))))


def _deduplicated_discrete_initial_candidates(
    candidates: list[tuple[float, numpy.ndarray]],
) -> list[tuple[float, numpy.ndarray]]:
    unique_candidates: list[tuple[float, numpy.ndarray]] = []
    seen: set[tuple[float, ...]] = set()
    for scale, log_parameters in candidates:
        clipped = numpy.clip(
            log_parameters,
            DISCRETE_LOG_PARAMETER_LOWER_BOUND,
            DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        )
        key = tuple(round(float(value), 12) for value in clipped.tolist())
        if key in seen:
            continue
        seen.add(key)
        unique_candidates.append((scale, clipped))
    return unique_candidates


def _perturbed_discrete_seed_candidates(
    log_parameters: numpy.ndarray,
    *,
    limit: int | None = None,
) -> list[tuple[float, numpy.ndarray]]:
    if log_parameters.size <= 1:
        return []
    position_span = numpy.linspace(-1.0, 1.0, log_parameters.size, dtype=float)
    alternating_span = numpy.array(
        [1.0 if index % 2 == 0 else -1.0 for index in range(log_parameters.size)],
        dtype=float,
    )
    perturbations = (
        0.35 * position_span,
        -0.35 * position_span,
        0.45 * alternating_span,
        -0.45 * alternating_span,
    )
    candidates: list[tuple[float, numpy.ndarray]] = []
    for offset in perturbations:
        candidate = numpy.clip(
            log_parameters + offset,
            DISCRETE_LOG_PARAMETER_LOWER_BOUND,
            DISCRETE_LOG_PARAMETER_UPPER_BOUND,
        )
        candidates.append(
            (
                _scale_from_log_parameters(candidate),
                candidate,
            )
        )
    if limit is not None:
        return candidates[:limit]
    return candidates


def _fit_symmetric_warm_start(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
    root_prior_mode: str,
    ascertainment_policy: str,
    equal_seed: numpy.ndarray,
    target_layout,
) -> numpy.ndarray | None:
    reciprocal_pairs = {
        pair
        for pair in allowed_transition_pairs
        if (pair[1], pair[0]) in allowed_transition_pairs
    }
    if not reciprocal_pairs or not _reciprocal_transition_graph_is_connected(
        state_count=len(state_order),
        reciprocal_pairs=reciprocal_pairs,
    ):
        return None
    symmetric_layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model="symmetric",
        state_ordering=state_ordering,
        allowed_transition_pairs=reciprocal_pairs,
    )
    ard_layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model="all-rates-different",
        state_ordering=state_ordering,
        allowed_transition_pairs=reciprocal_pairs,
    )
    symmetric_seed = numpy.array(
        broadcast_discrete_mk_log_parameters(
            build_discrete_mk_parameter_layout(
                state_order=state_order,
                model="equal-rates",
                state_ordering=state_ordering,
                allowed_transition_pairs=reciprocal_pairs,
            ),
            symmetric_layout,
            (float(equal_seed[0]),),
        ),
        dtype=float,
    )
    candidate_seeds = [symmetric_seed]
    candidate_seeds.extend(
        log_parameters
        for _, log_parameters in _perturbed_discrete_seed_candidates(
            symmetric_seed,
            limit=2,
        )
    )
    candidate_runs: list[_DiscreteOptimizationRun] = []
    for seed in candidate_seeds:
        try:
            candidate_runs.append(
                optimize_log_parameters(
                    tree,
                    states_by_taxon,
                    state_order=state_order,
                    model="symmetric",
                    state_ordering=state_ordering,
                    allowed_transition_pairs=reciprocal_pairs,
                    root_prior_mode=root_prior_mode,
                    ascertainment_policy=ascertainment_policy,
                    initial_log_parameters=seed,
                    initial_scale=_scale_from_log_parameters(seed),
                )
            )
        except PhylogeneticsError:
            continue
    if not candidate_runs:
        return None
    best_run = max(candidate_runs, key=lambda run: run.log_likelihood)
    reciprocal_seed = numpy.array(
        broadcast_discrete_mk_log_parameters(
            symmetric_layout,
            ard_layout,
            tuple(float(value) for value in best_run.log_parameters),
        ),
        dtype=float,
    )
    if target_layout.allowed_transition_pairs == ard_layout.allowed_transition_pairs:
        return reciprocal_seed
    merged_seed = equal_seed.copy()
    for arc in ard_layout.arcs:
        merged_seed[
            target_layout.parameter_index_for_pair(
                arc.source_index,
                arc.target_index,
            )
        ] = reciprocal_seed[arc.parameter_index]
    return merged_seed


def _reciprocal_transition_graph_is_connected(
    *,
    state_count: int,
    reciprocal_pairs: set[tuple[int, int]],
) -> bool:
    if state_count == 0:
        return False
    adjacency = {state_index: set() for state_index in range(state_count)}
    for source_index, target_index in reciprocal_pairs:
        adjacency[source_index].add(target_index)
    if any(not neighbors for neighbors in adjacency.values()):
        return False
    visited: set[int] = set()
    frontier = [0]
    while frontier:
        source_index = frontier.pop()
        if source_index in visited:
            continue
        visited.add(source_index)
        frontier.extend(
            target_index
            for target_index in adjacency[source_index]
            if target_index not in visited
        )
    return len(visited) == state_count


def _state_count_directional_seed(
    baseline_seed: numpy.ndarray,
    *,
    state_order: list[str],
    states_by_taxon: dict[str, str],
    layout,
) -> numpy.ndarray:
    state_counts = {
        state: sum(
            1 for observed_state in states_by_taxon.values() if observed_state == state
        )
        for state in state_order
    }
    directional_seed = baseline_seed.copy()
    for arc in layout.arcs:
        source_count = state_counts.get(arc.source_state, 0)
        target_count = state_counts.get(arc.target_state, 0)
        directional_seed[arc.parameter_index] += 0.5 * math.log(
            (target_count + 1.0) / (source_count + 1.0)
        )
    return numpy.clip(
        directional_seed,
        DISCRETE_LOG_PARAMETER_LOWER_BOUND,
        DISCRETE_LOG_PARAMETER_UPPER_BOUND,
    )


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

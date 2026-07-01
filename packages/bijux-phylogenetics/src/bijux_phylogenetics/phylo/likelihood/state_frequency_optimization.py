from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from typing import TypeVar

import numpy

from bijux_phylogenetics.core.categorical_probability_vectors import (
    build_categorical_probability_vector,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FiniteStateFrequencyOptimizationReport,
    LikelihoodOptimizationBoundaryWarning,
    StateFrequencyOptimizationRow,
)
from bijux_phylogenetics.phylo.likelihood.simplex_search import (
    run_simplex_likelihood_search,
)

_PayloadT = TypeVar("_PayloadT")
_BOUNDARY_TOLERANCE = 1e-9


def build_equal_input_rate_matrix(
    *,
    state_labels: Sequence[str],
    stationary_frequencies: Mapping[str, float] | Sequence[float],
    owner_name: str,
) -> numpy.ndarray:
    """Build one normalized equal-input CTMC generator from one stationary simplex."""
    validated_state_labels = tuple(state_labels)
    if isinstance(stationary_frequencies, Mapping):
        probability_vector = build_categorical_probability_vector(
            stationary_frequencies,
            expected_states=validated_state_labels,
        )
        stationary_vector = numpy.asarray(probability_vector.probabilities, dtype=float)
    else:
        stationary_vector = numpy.asarray(stationary_frequencies, dtype=float)
        if stationary_vector.shape != (len(validated_state_labels),):
            raise ValueError(
                f"{owner_name} requires exactly {len(validated_state_labels)} stationary-frequency values"
            )
        probability_vector = build_categorical_probability_vector(
            dict(
                zip(
                    validated_state_labels,
                    stationary_vector.tolist(),
                    strict=True,
                )
            ),
            expected_states=validated_state_labels,
        )
        stationary_vector = numpy.asarray(probability_vector.probabilities, dtype=float)
    rate_matrix = numpy.tile(stationary_vector, (len(validated_state_labels), 1))
    numpy.fill_diagonal(rate_matrix, 0.0)
    numpy.fill_diagonal(rate_matrix, -rate_matrix.sum(axis=1))
    return normalize_ctmc_rate_matrix_by_expected_substitution_rate(
        rate_matrix,
        stationary_vector,
        state_labels=validated_state_labels,
    )


def optimize_finite_state_stationary_frequencies(
    *,
    model_name: str,
    state_labels: Sequence[str],
    initial_values_by_name: Mapping[str, float],
    evaluate: Callable[[dict[str, float]], tuple[_PayloadT, float]],
    minimum_state_frequency: float = 1e-6,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
    nonconvergence_warning: str | None = None,
) -> FiniteStateFrequencyOptimizationReport:
    """Optimize one finite-state stationary simplex while keeping every component valid."""
    search_result = run_simplex_likelihood_search(
        initial_values_by_name=initial_values_by_name,
        expected_component_names=tuple(state_labels),
        owner_name=f"{model_name} stationary-frequency optimization",
        evaluate=evaluate,
        minimum_component_value=minimum_state_frequency,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    initial_mapping = search_result.initial_parameterization.constrained_mapping()
    _, initial_log_likelihood = evaluate(initial_mapping)
    optimized_mapping = search_result.optimized_parameterization.constrained_mapping()
    state_rows = [
        StateFrequencyOptimizationRow(
            state_label=state_label,
            initial_value=initial_mapping[state_label],
            optimized_value=optimized_mapping[state_label],
            hit_lower_bound=(
                abs(optimized_mapping[state_label] - minimum_state_frequency)
                <= _BOUNDARY_TOLERANCE
            ),
        )
        for state_label in state_labels
    ]
    warnings = []
    if not search_result.converged and nonconvergence_warning is not None:
        warnings.append(nonconvergence_warning)
    boundary_warning_rows = _build_state_frequency_boundary_warning_rows(
        state_rows,
        minimum_state_frequency=minimum_state_frequency,
    )
    warnings.extend(
        warning.message
        for warning in boundary_warning_rows
    )
    optimized_payload = search_result.payload
    return FiniteStateFrequencyOptimizationReport(
        model_name=model_name,
        taxa=list(getattr(optimized_payload, "taxa")),
        site_count=int(getattr(optimized_payload, "site_count")),
        pattern_count=int(getattr(optimized_payload, "pattern_count")),
        tree_newick=str(getattr(optimized_payload, "tree_newick")),
        state_count=int(getattr(optimized_payload, "state_count")),
        state_frequency_source="optimized",
        initial_log_likelihood=float(initial_log_likelihood),
        optimized_log_likelihood=float(search_result.objective_value),
        function_evaluation_count=search_result.function_evaluation_count + 1,
        optimization_pass_count=search_result.optimization_pass_count,
        converged=search_result.converged,
        minimum_state_frequency=minimum_state_frequency,
        state_rows=state_rows,
        warnings=warnings,
        boundary_warning_rows=boundary_warning_rows,
        parameter_diagnostics=dict(search_result.parameter_diagnostics),
    )


def _build_state_frequency_boundary_warning_rows(
    state_rows: Sequence[StateFrequencyOptimizationRow],
    *,
    minimum_state_frequency: float,
) -> list[LikelihoodOptimizationBoundaryWarning]:
    warnings: list[LikelihoodOptimizationBoundaryWarning] = []
    for row in state_rows:
        if row.hit_lower_bound:
            warnings.append(
                LikelihoodOptimizationBoundaryWarning(
                    warning_kind="frequency-boundary",
                    affected_parameter=f"state_frequency[{row.state_label}]",
                    boundary_side="lower",
                    observed_value=row.optimized_value,
                    lower_bound=minimum_state_frequency,
                    upper_bound=1.0,
                    affected_branch_id=None,
                    affected_branch_clade_id=None,
                    message=(
                        f"state_frequency[{row.state_label}] hit lower frequency boundary"
                    ),
                )
            )
    return warnings


__all__ = [
    "build_equal_input_rate_matrix",
    "optimize_finite_state_stationary_frequencies",
]

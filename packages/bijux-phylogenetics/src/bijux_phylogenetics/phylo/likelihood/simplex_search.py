from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
import math
from typing import Generic, TypeVar

import numpy

from bijux_phylogenetics.core.bounded_scalar_search import (
    BoundedScalarSearchDiagnostics,
)
from bijux_phylogenetics.core.categorical_probability_vectors import (
    build_categorical_probability_vector,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.simplex_coordinates import (
    SimplexCoordinateParameterization,
    parameterize_named_simplex,
    resolve_named_simplex_from_unconstrained,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_PayloadT = TypeVar("_PayloadT")


@dataclass(slots=True)
class SimplexLikelihoodSearchResult(Generic[_PayloadT]):
    """Best result from one simplex-constrained likelihood search."""

    initial_parameterization: SimplexCoordinateParameterization
    optimized_parameterization: SimplexCoordinateParameterization
    payload: _PayloadT
    objective_value: float
    function_evaluation_count: int
    optimization_pass_count: int
    converged: bool
    minimum_component_value: float
    coordinate_lower_bound: float
    coordinate_upper_bound: float
    parameter_diagnostics: dict[str, BoundedScalarSearchDiagnostics]


def resolve_simplex_coordinate_search_bounds(
    component_count: int,
    *,
    minimum_component_value: float,
    owner_name: str,
) -> tuple[float, float]:
    """Return symmetric additive log-ratio bounds induced by one simplex floor."""
    if component_count < 2:
        raise PhylogeneticsError(
            f"{owner_name} requires at least two simplex components",
            code="simplex_search_component_count_invalid",
            details={"component_count": component_count},
        )
    validated_minimum_component_value = _validate_minimum_component_value(
        component_count=component_count,
        minimum_component_value=minimum_component_value,
        owner_name=owner_name,
    )
    maximum_component_value = 1.0 - (
        (component_count - 1) * validated_minimum_component_value
    )
    coordinate_bound = math.log(
        maximum_component_value / validated_minimum_component_value
    )
    return (-coordinate_bound, coordinate_bound)


def run_simplex_likelihood_search(
    *,
    initial_values_by_name: Mapping[str, float],
    expected_component_names: Sequence[str],
    owner_name: str,
    evaluate: Callable[[dict[str, float]], tuple[_PayloadT, float]],
    reference_component_name: str | None = None,
    normalization_tolerance: float = 1e-9,
    minimum_component_value: float = 1e-6,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> SimplexLikelihoodSearchResult[_PayloadT]:
    """Maximize one likelihood-like objective over one strict simplex."""
    initial_parameterization = _resolve_initial_simplex_parameterization(
        initial_values_by_name,
        expected_component_names=expected_component_names,
        owner_name=owner_name,
        reference_component_name=reference_component_name,
        normalization_tolerance=normalization_tolerance,
        minimum_component_value=minimum_component_value,
    )
    coordinate_component_names = tuple(
        component_name
        for component_name in initial_parameterization.component_names
        if component_name != initial_parameterization.reference_component_name
    )
    coordinate_lower_bound, coordinate_upper_bound = (
        resolve_simplex_coordinate_search_bounds(
            len(initial_parameterization.component_names),
            minimum_component_value=minimum_component_value,
            owner_name=owner_name,
        )
    )
    search_result = run_bounded_coordinate_likelihood_search(
        initial_values=dict(
            zip(
                coordinate_component_names,
                initial_parameterization.unconstrained_values,
                strict=True,
            )
        ),
        bounds_by_name=dict.fromkeys(
            coordinate_component_names, (coordinate_lower_bound, coordinate_upper_bound)
        ),
        evaluate=lambda coordinate_values: _evaluate_simplex_coordinates(
            coordinate_values,
            component_names=initial_parameterization.component_names,
            reference_component_name=initial_parameterization.reference_component_name,
            owner_name=owner_name,
            normalization_tolerance=normalization_tolerance,
            minimum_component_value=minimum_component_value,
            evaluate=evaluate,
        ),
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    optimized_parameterization = _lift_parameterization_to_minimum_component_value(
        resolve_named_simplex_from_unconstrained(
            [
                search_result.parameter_values[coordinate_name]
                for coordinate_name in coordinate_component_names
            ],
            component_names=initial_parameterization.component_names,
            owner_name=owner_name,
            reference_component_name=initial_parameterization.reference_component_name,
            normalization_tolerance=normalization_tolerance,
        ),
        minimum_component_value=minimum_component_value,
        owner_name=owner_name,
    )
    return SimplexLikelihoodSearchResult(
        initial_parameterization=initial_parameterization,
        optimized_parameterization=optimized_parameterization,
        payload=search_result.payload,
        objective_value=search_result.objective_value,
        function_evaluation_count=search_result.function_evaluation_count,
        optimization_pass_count=search_result.optimization_pass_count,
        converged=search_result.converged,
        minimum_component_value=minimum_component_value,
        coordinate_lower_bound=coordinate_lower_bound,
        coordinate_upper_bound=coordinate_upper_bound,
        parameter_diagnostics=dict(search_result.parameter_diagnostics),
    )


def _resolve_initial_simplex_parameterization(
    initial_values_by_name: Mapping[str, float],
    *,
    expected_component_names: Sequence[str],
    owner_name: str,
    reference_component_name: str | None,
    normalization_tolerance: float,
    minimum_component_value: float,
) -> SimplexCoordinateParameterization:
    component_names = tuple(expected_component_names)
    validated_minimum_component_value = _validate_minimum_component_value(
        component_count=len(component_names),
        minimum_component_value=minimum_component_value,
        owner_name=owner_name,
    )
    validated_vector = build_categorical_probability_vector(
        initial_values_by_name,
        expected_states=component_names,
        normalization_tolerance=normalization_tolerance,
    )
    initial_probabilities = numpy.asarray(validated_vector.probabilities, dtype=float)
    if numpy.all(initial_probabilities >= validated_minimum_component_value):
        resolved_probabilities = initial_probabilities
    else:
        residual_mass = 1.0 - (len(component_names) * validated_minimum_component_value)
        resolved_probabilities = validated_minimum_component_value + (
            residual_mass * initial_probabilities
        )
    return parameterize_named_simplex(
        dict(zip(component_names, resolved_probabilities.tolist(), strict=True)),
        expected_component_names=component_names,
        owner_name=owner_name,
        reference_component_name=reference_component_name,
        normalization_tolerance=normalization_tolerance,
    )


def _evaluate_simplex_coordinates(
    coordinate_values: Mapping[str, float],
    *,
    component_names: tuple[str, ...],
    reference_component_name: str,
    owner_name: str,
    normalization_tolerance: float,
    minimum_component_value: float,
    evaluate: Callable[[dict[str, float]], tuple[_PayloadT, float]],
) -> tuple[_PayloadT, float]:
    coordinate_component_names = tuple(
        component_name
        for component_name in component_names
        if component_name != reference_component_name
    )
    parameterization = _lift_parameterization_to_minimum_component_value(
        resolve_named_simplex_from_unconstrained(
            [
                float(coordinate_values[coordinate_name])
                for coordinate_name in coordinate_component_names
            ],
            component_names=component_names,
            owner_name=owner_name,
            reference_component_name=reference_component_name,
            normalization_tolerance=normalization_tolerance,
        ),
        minimum_component_value=minimum_component_value,
        owner_name=owner_name,
    )
    return evaluate(parameterization.constrained_mapping())


def _lift_parameterization_to_minimum_component_value(
    parameterization: SimplexCoordinateParameterization,
    *,
    minimum_component_value: float,
    owner_name: str,
) -> SimplexCoordinateParameterization:
    constrained_values = numpy.asarray(parameterization.constrained_values, dtype=float)
    if numpy.all(constrained_values >= minimum_component_value):
        return parameterization
    residual_mass = 1.0 - (
        len(parameterization.component_names) * minimum_component_value
    )
    adjusted_values = minimum_component_value + (residual_mass * constrained_values)
    adjusted_parameterization = parameterize_named_simplex(
        dict(
            zip(
                parameterization.component_names,
                adjusted_values.tolist(),
                strict=True,
            )
        ),
        expected_component_names=parameterization.component_names,
        owner_name=owner_name,
        reference_component_name=parameterization.reference_component_name,
        normalization_tolerance=parameterization.normalization_tolerance,
    )
    for component_name, component_value in zip(
        adjusted_parameterization.component_names,
        adjusted_parameterization.constrained_values,
        strict=True,
    ):
        if component_value < minimum_component_value:
            raise PhylogeneticsError(
                f"{owner_name} simplex search could not enforce the declared minimum component value",
                code="simplex_search_component_floor_violated",
                details={
                    "component_name": component_name,
                    "component_value": component_value,
                    "minimum_component_value": minimum_component_value,
                },
            )
    return adjusted_parameterization


def _validate_minimum_component_value(
    *,
    component_count: int,
    minimum_component_value: float,
    owner_name: str,
) -> float:
    validated_minimum_component_value = float(minimum_component_value)
    if not math.isfinite(validated_minimum_component_value):
        raise PhylogeneticsError(
            f"{owner_name} minimum component value must be finite",
            code="simplex_search_component_floor_not_finite",
        )
    if validated_minimum_component_value <= 0.0:
        raise PhylogeneticsError(
            f"{owner_name} minimum component value must be positive",
            code="simplex_search_component_floor_not_positive",
            details={
                "minimum_component_value": validated_minimum_component_value,
            },
        )
    maximum_minimum_component_value = 1.0 / component_count
    if validated_minimum_component_value >= maximum_minimum_component_value:
        raise PhylogeneticsError(
            f"{owner_name} minimum component value must stay below the uniform simplex mass",
            code="simplex_search_component_floor_too_large",
            details={
                "component_count": component_count,
                "minimum_component_value": validated_minimum_component_value,
                "maximum_minimum_component_value": maximum_minimum_component_value,
            },
        )
    return validated_minimum_component_value


__all__ = [
    "SimplexLikelihoodSearchResult",
    "resolve_simplex_coordinate_search_bounds",
    "run_simplex_likelihood_search",
]

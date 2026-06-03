from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.core.categorical_probability_vectors import (
    build_categorical_probability_vector,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


@dataclass(frozen=True, slots=True)
class SimplexCoordinateParameterization:
    """One named simplex vector together with one additive log-ratio coordinate view."""

    component_names: tuple[str, ...]
    constrained_values: tuple[float, ...]
    unconstrained_values: tuple[float, ...]
    reference_component_name: str
    normalization_tolerance: float

    def constrained_mapping(self) -> dict[str, float]:
        return dict(zip(self.component_names, self.constrained_values, strict=True))


def parameterize_named_simplex(
    values_by_name: Mapping[str, float],
    *,
    expected_component_names: Sequence[str],
    owner_name: str,
    reference_component_name: str | None = None,
    normalization_tolerance: float = 1e-9,
) -> SimplexCoordinateParameterization:
    component_names = tuple(expected_component_names)
    validated_reference_component_name = _validate_reference_component_name(
        component_names=component_names,
        reference_component_name=reference_component_name,
        owner_name=owner_name,
    )
    validated_vector = build_categorical_probability_vector(
        values_by_name,
        expected_states=component_names,
        normalization_tolerance=normalization_tolerance,
    )
    constrained_values = validated_vector.probabilities
    _require_strictly_positive_components(
        component_names=component_names,
        constrained_values=constrained_values,
        owner_name=owner_name,
    )
    reference_value = constrained_values[
        component_names.index(validated_reference_component_name)
    ]
    unconstrained_values = tuple(
        math.log(constrained_values[index] / reference_value)
        for index, component_name in enumerate(component_names)
        if component_name != validated_reference_component_name
    )
    return SimplexCoordinateParameterization(
        component_names=component_names,
        constrained_values=constrained_values,
        unconstrained_values=unconstrained_values,
        reference_component_name=validated_reference_component_name,
        normalization_tolerance=normalization_tolerance,
    )


def resolve_named_simplex_from_unconstrained(
    unconstrained_values: Sequence[float],
    *,
    component_names: Sequence[str],
    owner_name: str,
    reference_component_name: str | None = None,
    normalization_tolerance: float = 1e-9,
) -> SimplexCoordinateParameterization:
    resolved_component_names = tuple(component_names)
    validated_reference_component_name = _validate_reference_component_name(
        component_names=resolved_component_names,
        reference_component_name=reference_component_name,
        owner_name=owner_name,
    )
    expected_dimension = len(resolved_component_names) - 1
    if len(unconstrained_values) != expected_dimension:
        raise PhylogeneticsError(
            f"{owner_name} requires exactly {expected_dimension} unconstrained simplex coordinates",
            code="simplex_coordinate_dimension_mismatch",
            details={
                "expected_dimension": expected_dimension,
                "observed_dimension": len(unconstrained_values),
            },
        )
    log_weights: list[float] = []
    unconstrained_index = 0
    for component_name in resolved_component_names:
        if component_name == validated_reference_component_name:
            log_weights.append(0.0)
            continue
        coordinate = float(unconstrained_values[unconstrained_index])
        if not math.isfinite(coordinate):
            raise PhylogeneticsError(
                f"{owner_name} unconstrained simplex coordinates must all be finite",
                code="simplex_coordinate_not_finite",
                details={
                    "component_name": component_name,
                    "coordinate_index": unconstrained_index,
                    "coordinate_value": coordinate,
                },
            )
        log_weights.append(coordinate)
        unconstrained_index += 1
    maximum_log_weight = max(log_weights)
    scaled_weights = [math.exp(weight - maximum_log_weight) for weight in log_weights]
    normalization = math.fsum(scaled_weights)
    constrained_values = tuple(weight / normalization for weight in scaled_weights)
    return parameterize_named_simplex(
        dict(zip(resolved_component_names, constrained_values, strict=True)),
        expected_component_names=resolved_component_names,
        owner_name=owner_name,
        reference_component_name=validated_reference_component_name,
        normalization_tolerance=normalization_tolerance,
    )


def _validate_reference_component_name(
    *,
    component_names: tuple[str, ...],
    reference_component_name: str | None,
    owner_name: str,
) -> str:
    if not component_names:
        raise PhylogeneticsError(
            f"{owner_name} requires at least two simplex components",
            code="simplex_coordinate_component_list_empty",
        )
    validated_reference_component_name = (
        component_names[-1]
        if reference_component_name is None
        else reference_component_name
    )
    if validated_reference_component_name not in component_names:
        raise PhylogeneticsError(
            f"{owner_name} reference component must belong to the simplex component set",
            code="simplex_coordinate_reference_invalid",
            details={
                "reference_component_name": validated_reference_component_name,
                "component_names": list(component_names),
            },
        )
    if len(component_names) < 2:
        raise PhylogeneticsError(
            f"{owner_name} requires at least two simplex components",
            code="simplex_coordinate_component_count_invalid",
            details={"component_count": len(component_names)},
        )
    return validated_reference_component_name


def _require_strictly_positive_components(
    *,
    component_names: tuple[str, ...],
    constrained_values: tuple[float, ...],
    owner_name: str,
) -> None:
    for component_name, constrained_value in zip(
        component_names,
        constrained_values,
        strict=True,
    ):
        if constrained_value <= 0.0:
            raise PhylogeneticsError(
                f"{owner_name} requires strictly positive simplex components",
                code="simplex_coordinate_component_not_positive",
                details={
                    "component_name": component_name,
                    "component_value": constrained_value,
                },
            )


__all__ = [
    "SimplexCoordinateParameterization",
    "parameterize_named_simplex",
    "resolve_named_simplex_from_unconstrained",
]

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math

from bijux_phylogenetics.runtime.errors import PhylogeneticsError

CATEGORICAL_MISSING_STATE_POLICIES = (
    "fill-zero",
    "reject",
)


@dataclass(frozen=True, slots=True)
class CategoricalProbabilityVector:
    """One validated categorical probability vector over named states."""

    states: tuple[str, ...]
    probabilities: tuple[float, ...]
    missing_state_policy: str
    normalization_tolerance: float

    def probability_for(self, state: str) -> float:
        for state_label, probability in zip(self.states, self.probabilities, strict=True):
            if state_label == state:
                return probability
        raise KeyError(state)

    @property
    def total_probability(self) -> float:
        return math.fsum(self.probabilities)

    def as_mapping(self) -> dict[str, float]:
        return dict(zip(self.states, self.probabilities, strict=True))


def build_categorical_probability_vector(
    probabilities_by_state: Mapping[str, float],
    *,
    expected_states: Sequence[str] | None = None,
    missing_state_policy: str = "reject",
    normalization_tolerance: float = 1e-9,
) -> CategoricalProbabilityVector:
    validated_policy = _validate_missing_state_policy(missing_state_policy)
    validated_tolerance = _validate_normalization_tolerance(normalization_tolerance)
    if not probabilities_by_state:
        raise PhylogeneticsError(
            "categorical probability vector requires at least one state",
            code="categorical_probability_vector_empty",
        )

    state_order = _resolve_state_order(
        observed_states=tuple(probabilities_by_state),
        expected_states=expected_states,
        missing_state_policy=validated_policy,
    )
    probabilities = tuple(
        _validate_state_probability(
            state=state,
            probability=probabilities_by_state.get(state, 0.0),
        )
        for state in state_order
    )
    total_probability = math.fsum(probabilities)
    if not math.isclose(
        total_probability,
        1.0,
        rel_tol=0.0,
        abs_tol=validated_tolerance,
    ):
        raise PhylogeneticsError(
            "categorical probability vector must sum to one within tolerance",
            code="categorical_probability_vector_not_normalized",
            details={
                "total_probability": total_probability,
                "expected_total": 1.0,
                "absolute_error": abs(total_probability - 1.0),
                "normalization_tolerance": validated_tolerance,
            },
        )
    return CategoricalProbabilityVector(
        states=state_order,
        probabilities=probabilities,
        missing_state_policy=validated_policy,
        normalization_tolerance=validated_tolerance,
    )


def _resolve_state_order(
    *,
    observed_states: tuple[str, ...],
    expected_states: Sequence[str] | None,
    missing_state_policy: str,
) -> tuple[str, ...]:
    if expected_states is None:
        return observed_states

    expected_state_order = tuple(expected_states)
    if not expected_state_order:
        raise PhylogeneticsError(
            "expected state list must not be empty",
            code="categorical_probability_vector_expected_states_empty",
        )
    if len(set(expected_state_order)) != len(expected_state_order):
        raise PhylogeneticsError(
            "expected state list must not contain duplicates",
            code="categorical_probability_vector_expected_states_duplicate",
            details={"expected_states": list(expected_state_order)},
        )
    unexpected_states = [
        state for state in observed_states if state not in set(expected_state_order)
    ]
    if unexpected_states:
        raise PhylogeneticsError(
            "categorical probability vector contains unexpected states",
            code="categorical_probability_vector_unexpected_states",
            details={"unexpected_states": unexpected_states},
        )
    missing_states = [
        state for state in expected_state_order if state not in set(observed_states)
    ]
    if missing_states and missing_state_policy == "reject":
        raise PhylogeneticsError(
            "categorical probability vector is missing expected states",
            code="categorical_probability_vector_missing_states",
            details={
                "missing_states": missing_states,
                "missing_state_policy": missing_state_policy,
            },
        )
    return expected_state_order


def _validate_missing_state_policy(missing_state_policy: str) -> str:
    if missing_state_policy not in CATEGORICAL_MISSING_STATE_POLICIES:
        raise PhylogeneticsError(
            "categorical probability vector missing-state policy is unsupported",
            code="categorical_probability_vector_missing_state_policy_invalid",
            details={
                "missing_state_policy": missing_state_policy,
                "allowed_policies": list(CATEGORICAL_MISSING_STATE_POLICIES),
            },
        )
    return missing_state_policy


def _validate_normalization_tolerance(normalization_tolerance: float) -> float:
    if (
        math.isnan(normalization_tolerance)
        or math.isinf(normalization_tolerance)
        or normalization_tolerance < 0.0
    ):
        raise PhylogeneticsError(
            "categorical probability vector normalization tolerance must be finite and non-negative",
            code="categorical_probability_vector_tolerance_invalid",
            details={"normalization_tolerance": normalization_tolerance},
        )
    return float(normalization_tolerance)


def _validate_state_probability(*, state: str, probability: float) -> float:
    if math.isnan(probability) or math.isinf(probability):
        raise PhylogeneticsError(
            "categorical probability vector values must be finite",
            code="categorical_probability_vector_value_not_finite",
            details={"state": state, "probability": probability},
        )
    if probability < 0.0:
        raise PhylogeneticsError(
            "categorical probability vector values must be non-negative",
            code="categorical_probability_vector_value_negative",
            details={"state": state, "probability": probability},
        )
    return float(probability)


__all__ = [
    "CATEGORICAL_MISSING_STATE_POLICIES",
    "CategoricalProbabilityVector",
    "build_categorical_probability_vector",
]

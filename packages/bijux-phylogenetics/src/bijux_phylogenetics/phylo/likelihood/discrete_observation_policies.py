from __future__ import annotations

import numpy

from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

DISCRETE_AMBIGUITY_SEPARATOR = "|"
DISCRETE_MISSING_OBSERVATION_SYMBOLS = frozenset({"?", "-"})
_SUPPORTED_DISCRETE_OBSERVATION_POLICIES = frozenset(
    {"reject", "treat-as-missing", "ambiguity-vector"}
)


def validate_discrete_observation_policy(
    observation_policy: str,
    *,
    owner_name: str,
) -> str:
    normalized_policy = observation_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_DISCRETE_OBSERVATION_POLICIES:
        raise ValueError(
            f"{owner_name} observation_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_DISCRETE_OBSERVATION_POLICIES))
        )
    return normalized_policy


def normalize_discrete_observation_token(observation: str) -> str:
    return observation.strip()


def is_missing_discrete_observation_token(observation: str) -> bool:
    return (
        not observation
        or normalize_discrete_observation_token(observation)
        in DISCRETE_MISSING_OBSERVATION_SYMBOLS
    )


def parse_discrete_ambiguity_token(observation: str) -> tuple[str, ...] | None:
    normalized_observation = normalize_discrete_observation_token(observation)
    if DISCRETE_AMBIGUITY_SEPARATOR not in normalized_observation:
        return None
    parts = tuple(part.strip() for part in normalized_observation.split("|"))
    if not parts or any(not part for part in parts):
        raise InvalidAlignmentError(
            f"discrete likelihood encountered malformed ambiguity token '{normalized_observation}'"
        )
    unique_parts = tuple(dict.fromkeys(parts))
    if len(unique_parts) < 2:
        raise InvalidAlignmentError(
            f"discrete likelihood ambiguity token '{normalized_observation}' must contain at least two distinct states"
        )
    return unique_parts


def resolve_discrete_observation_leaf_vector(
    observation: str,
    *,
    model_name: str,
    node_name: str,
    state_order: list[str],
    observation_policy: str,
) -> numpy.ndarray:
    validated_policy = validate_discrete_observation_policy(
        observation_policy,
        owner_name=f"{model_name} likelihood",
    )
    normalized_observation = normalize_discrete_observation_token(observation)
    state_index = {state: index for index, state in enumerate(state_order)}
    if normalized_observation in state_index:
        vector = numpy.zeros(len(state_order), dtype=float)
        vector[state_index[normalized_observation]] = 1.0
        return vector
    if is_missing_discrete_observation_token(normalized_observation):
        if validated_policy == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood observation policy 'reject' does not allow missing discrete state '{normalized_observation}' at taxon '{node_name}'"
            )
        return numpy.ones(len(state_order), dtype=float)
    ambiguous_states = parse_discrete_ambiguity_token(normalized_observation)
    if ambiguous_states is not None:
        if validated_policy == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood observation policy 'reject' does not allow ambiguous discrete state '{normalized_observation}' at taxon '{node_name}'"
            )
        if validated_policy == "treat-as-missing":
            return numpy.ones(len(state_order), dtype=float)
        unknown_states = [
            state for state in ambiguous_states if state not in state_index
        ]
        if unknown_states:
            raise InvalidAlignmentError(
                f"{model_name} likelihood ambiguity token '{normalized_observation}' contains undeclared states {', '.join(sorted(unknown_states))}"
            )
        vector = numpy.zeros(len(state_order), dtype=float)
        for state in ambiguous_states:
            vector[state_index[state]] = 1.0
        return vector
    raise InvalidAlignmentError(
        f"{model_name} likelihood encountered unsupported discrete state '{normalized_observation}' at taxon '{node_name}'"
    )


__all__ = [
    "DISCRETE_AMBIGUITY_SEPARATOR",
    "DISCRETE_MISSING_OBSERVATION_SYMBOLS",
    "is_missing_discrete_observation_token",
    "normalize_discrete_observation_token",
    "parse_discrete_ambiguity_token",
    "resolve_discrete_observation_leaf_vector",
    "validate_discrete_observation_policy",
]

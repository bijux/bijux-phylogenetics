from __future__ import annotations

from .._statistics import _round_float


def _normalize_discrete_states(states: list[str]) -> tuple[str, ...]:
    unique_states = tuple(dict.fromkeys(state for state in states if state))
    if len(unique_states) < 2:
        raise ValueError("states must contain at least two distinct non-empty states")
    return unique_states


def _normalize_root_state_probabilities(
    *,
    states: tuple[str, ...],
    root_state: str | None,
    root_state_probabilities: dict[str, float] | None,
) -> dict[str, float]:
    if root_state is not None and root_state_probabilities is not None:
        raise ValueError(
            "root_state and root_state_probabilities cannot be supplied together"
        )
    state_set = set(states)
    if root_state is not None:
        if root_state not in state_set:
            raise ValueError(f"root_state '{root_state}' is not present in states")
        return {state: 1.0 if state == root_state else 0.0 for state in states}
    if root_state_probabilities is None:
        probability = 1.0 / len(states)
        return {state: _round_float(probability) for state in states}
    unknown_states = set(root_state_probabilities).difference(state_set)
    if unknown_states:
        unknown_state = sorted(unknown_states)[0]
        raise ValueError(
            f"root_state_probabilities contains unknown state '{unknown_state}'"
        )
    probabilities = {
        state: float(root_state_probabilities.get(state, 0.0)) for state in states
    }
    if any(value < 0.0 for value in probabilities.values()):
        raise ValueError("root_state_probabilities cannot contain negative values")
    total = sum(probabilities.values())
    if total <= 0.0:
        raise ValueError("root_state_probabilities must sum to a positive value")
    return {
        state: _round_float(value / total) for state, value in probabilities.items()
    }

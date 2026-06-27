from __future__ import annotations

from collections.abc import Iterable

import numpy


def build_state_selection_vector(
    *,
    state_count: int,
    selected_indices: Iterable[int],
) -> numpy.ndarray:
    """Materialize one explicit partial-likelihood vector over selected states."""
    if state_count <= 0:
        raise ValueError("state_count must be positive")
    vector = numpy.zeros(state_count, dtype=float)
    selected_any_state = False
    for index in selected_indices:
        if index < 0 or index >= state_count:
            raise ValueError(
                f"selected state index {index} is outside the supported state space"
            )
        vector[index] = 1.0
        selected_any_state = True
    if not selected_any_state:
        raise ValueError("selected_indices must contain at least one state index")
    return vector


def build_all_states_observation_vector(
    *,
    state_count: int,
) -> numpy.ndarray:
    """Materialize one explicit observation vector that allows every state."""
    if state_count <= 0:
        raise ValueError("state_count must be positive")
    return numpy.ones(state_count, dtype=float)


__all__ = [
    "build_all_states_observation_vector",
    "build_state_selection_vector",
]

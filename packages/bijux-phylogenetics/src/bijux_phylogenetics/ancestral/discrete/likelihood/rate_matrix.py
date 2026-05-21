from __future__ import annotations

import math

import numpy

from ..models import DiscreteTransitionRateRow
from ..policy import transition_allowed


def rate_matrix_from_log_parameters(
    log_parameters: numpy.ndarray,
    *,
    state_order: list[str],
    model: str,
    state_ordering: str,
    allowed_transition_pairs: set[tuple[int, int]],
) -> numpy.ndarray:
    state_count = len(state_order)
    rate_matrix = numpy.zeros((state_count, state_count), dtype=float)
    parameter_index = 0
    if model == "equal-rates":
        rate = math.exp(float(log_parameters[0]))
        for left_index in range(state_count):
            for right_index in range(state_count):
                if transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    rate_matrix[left_index, right_index] = rate
    elif model == "symmetric":
        for left_index in range(state_count):
            for right_index in range(left_index + 1, state_count):
                if not transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                if not transition_allowed(
                    right_index,
                    left_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                rate = math.exp(float(log_parameters[parameter_index]))
                parameter_index += 1
                rate_matrix[left_index, right_index] = rate
                rate_matrix[right_index, left_index] = rate
    else:
        for left_index in range(state_count):
            for right_index in range(state_count):
                if not transition_allowed(
                    left_index,
                    right_index,
                    state_count=state_count,
                    state_ordering=state_ordering,
                    allowed_transition_pairs=allowed_transition_pairs,
                ):
                    continue
                rate_matrix[left_index, right_index] = math.exp(
                    float(log_parameters[parameter_index])
                )
                parameter_index += 1
    for state_index in range(state_count):
        rate_matrix[state_index, state_index] = -float(
            numpy.sum(rate_matrix[state_index, :])
        )
    return rate_matrix


def build_transition_rate_rows(
    *,
    state_order: list[str],
    state_ordering: str,
    rate_matrix: numpy.ndarray,
    allowed_transition_pairs: set[tuple[int, int]],
) -> list[DiscreteTransitionRateRow]:
    rows: list[DiscreteTransitionRateRow] = []
    state_count = len(state_order)
    for left_index, source_state in enumerate(state_order):
        for right_index, target_state in enumerate(state_order):
            if left_index == right_index:
                continue
            rows.append(
                DiscreteTransitionRateRow(
                    source_state=source_state,
                    target_state=target_state,
                    transition_allowed=transition_allowed(
                        left_index,
                        right_index,
                        state_count=state_count,
                        state_ordering=state_ordering,
                        allowed_transition_pairs=allowed_transition_pairs,
                    ),
                    step_distance=abs(left_index - right_index),
                    rate=float(rate_matrix[left_index, right_index]),
                )
            )
    return rows

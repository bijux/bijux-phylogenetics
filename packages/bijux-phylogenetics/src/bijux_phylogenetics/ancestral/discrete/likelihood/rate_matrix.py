from __future__ import annotations

import math

import numpy

from ..models import DiscreteTransitionRateRow
from .parameter_layout import build_discrete_mk_parameter_layout


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
    layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model=model,
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    for arc in layout.arcs:
        rate_matrix[arc.source_index, arc.target_index] = math.exp(
            float(log_parameters[arc.parameter_index])
        )
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
    layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model="all-rates-different",
        state_ordering=state_ordering,
        allowed_transition_pairs=allowed_transition_pairs,
    )
    allowed_pairs = {(arc.source_index, arc.target_index) for arc in layout.arcs}
    rows: list[DiscreteTransitionRateRow] = []
    for left_index, source_state in enumerate(state_order):
        for right_index, target_state in enumerate(state_order):
            if left_index == right_index:
                continue
            rows.append(
                DiscreteTransitionRateRow(
                    source_state=source_state,
                    target_state=target_state,
                    transition_allowed=(left_index, right_index) in allowed_pairs,
                    step_distance=abs(left_index - right_index),
                    rate=float(rate_matrix[left_index, right_index]),
                )
            )
    return rows

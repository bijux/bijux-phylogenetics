from __future__ import annotations

from collections.abc import Mapping, Sequence
import math

import numpy


def two_tip_log_likelihood_for_labels(
    *,
    observed_pairs: Sequence[tuple[str, str]],
    state_index: Mapping[str, int],
    root_prior: numpy.ndarray,
    left_transition: numpy.ndarray,
    right_transition: numpy.ndarray,
) -> float:
    """Return one exact two-tip log-likelihood from explicit state labels."""
    probability = 1.0
    for left_state, right_state in observed_pairs:
        probability *= two_tip_site_probability(
            left_index=state_index[left_state],
            right_index=state_index[right_state],
            root_prior=root_prior,
            left_transition=left_transition,
            right_transition=right_transition,
        )
    return math.log(probability)


def two_tip_site_probability(
    *,
    left_index: int,
    right_index: int,
    root_prior: numpy.ndarray,
    left_transition: numpy.ndarray,
    right_transition: numpy.ndarray,
) -> float:
    """Return one explicit two-tip site probability under a fixed root prior."""
    return float(
        numpy.sum(
            root_prior
            * left_transition[:, left_index]
            * right_transition[:, right_index]
        )
    )

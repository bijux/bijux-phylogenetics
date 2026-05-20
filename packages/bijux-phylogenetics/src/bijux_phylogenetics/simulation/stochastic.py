from __future__ import annotations

from math import exp
import random


def _poisson_count(expected_changes: float, rng: random.Random) -> int:
    if expected_changes <= 0.0:
        return 0
    threshold = exp(-expected_changes)
    product = 1.0
    changes = 0
    while product > threshold:
        changes += 1
        product *= rng.random()
    return changes - 1

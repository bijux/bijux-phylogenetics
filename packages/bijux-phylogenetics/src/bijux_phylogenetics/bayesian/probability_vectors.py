from __future__ import annotations

from dataclasses import dataclass
import math

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


__all__ = [
    "CATEGORICAL_MISSING_STATE_POLICIES",
    "CategoricalProbabilityVector",
]

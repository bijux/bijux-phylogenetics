from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, TypeGuard


class BayesianPriorModelProtocol(Protocol):
    """Shared structural interface for one Bayesian prior model."""

    family: str

    def parameter_values(self) -> dict[str, float]:
        """Return the explicit numeric parameters that define the prior."""


@dataclass(frozen=True, slots=True)
class BayesianPriorModelSummary:
    """One shared summary of one Bayesian prior model."""

    family: str
    parameter_values: dict[str, float]
    parameter_count: int


def summarize_bayesian_prior_model(
    prior_model: BayesianPriorModelProtocol,
) -> BayesianPriorModelSummary:
    """Build one protocol-level summary of one prior model."""
    parameter_values = dict(prior_model.parameter_values())
    return BayesianPriorModelSummary(
        family=prior_model.family,
        parameter_values=parameter_values,
        parameter_count=len(parameter_values),
    )


def supports_bayesian_prior_model_protocol(
    value: object,
) -> TypeGuard[BayesianPriorModelProtocol]:
    """Check whether one value exposes the shared Bayesian prior interface."""
    return isinstance(getattr(value, "family", None), str) and callable(
        getattr(value, "parameter_values", None)
    )


__all__ = [
    "BayesianPriorModelProtocol",
    "BayesianPriorModelSummary",
    "summarize_bayesian_prior_model",
    "supports_bayesian_prior_model_protocol",
]

from __future__ import annotations

import math
from typing import Collection

from bijux_phylogenetics.phylo.likelihood.models import (
    NucleotideLikelihoodSearchConvergenceDecision,
)


def validate_nucleotide_likelihood_search_improvement_tolerance(
    tolerance: float,
) -> float:
    """Validate one native likelihood tree-search acceptance tolerance."""
    if not math.isfinite(tolerance) or tolerance < 0.0:
        raise ValueError(
            "native likelihood search improvement tolerance must be finite and nonnegative"
        )
    return float(tolerance)


def resolve_nucleotide_likelihood_search_convergence_decision(
    *,
    best_improving_delta: float | None,
    improvement_tolerance: float,
    budget_stopping_reason: str | None = None,
    candidate_topology_fingerprint: str | None = None,
    seen_topology_fingerprints: Collection[str] | None = None,
    failure_detected: bool = False,
) -> NucleotideLikelihoodSearchConvergenceDecision:
    """Resolve whether one native likelihood topology search must stop."""
    validated_tolerance = validate_nucleotide_likelihood_search_improvement_tolerance(
        improvement_tolerance
    )
    if failure_detected:
        return NucleotideLikelihoodSearchConvergenceDecision(
            should_stop=True,
            stopping_reason="search-failure",
        )
    if budget_stopping_reason is not None:
        return NucleotideLikelihoodSearchConvergenceDecision(
            should_stop=True,
            stopping_reason=budget_stopping_reason,
        )
    if best_improving_delta is None or best_improving_delta <= 0.0:
        return NucleotideLikelihoodSearchConvergenceDecision(
            should_stop=True,
            stopping_reason="no-improving-neighbor",
        )
    if best_improving_delta <= validated_tolerance:
        return NucleotideLikelihoodSearchConvergenceDecision(
            should_stop=True,
            stopping_reason="improvement-within-tolerance",
        )
    if (
        candidate_topology_fingerprint is not None
        and seen_topology_fingerprints is not None
        and candidate_topology_fingerprint in seen_topology_fingerprints
    ):
        return NucleotideLikelihoodSearchConvergenceDecision(
            should_stop=True,
            stopping_reason="repeated-topology",
        )
    return NucleotideLikelihoodSearchConvergenceDecision(
        should_stop=False,
        stopping_reason=None,
    )

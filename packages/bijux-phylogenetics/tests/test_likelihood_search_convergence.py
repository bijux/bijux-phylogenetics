from __future__ import annotations

import pytest

from bijux_phylogenetics.phylo.likelihood.search_convergence import (
    resolve_nucleotide_likelihood_search_convergence_decision,
    validate_nucleotide_likelihood_search_improvement_tolerance,
)


def test_validate_nucleotide_likelihood_search_improvement_tolerance_rejects_negative_values() -> None:
    with pytest.raises(ValueError, match="must be finite and nonnegative"):
        validate_nucleotide_likelihood_search_improvement_tolerance(-1e-9)


def test_likelihood_search_convergence_reports_no_improving_neighbor() -> None:
    decision = resolve_nucleotide_likelihood_search_convergence_decision(
        best_improving_delta=None,
        improvement_tolerance=1e-9,
    )

    assert decision.should_stop is True
    assert decision.stopping_reason == "no-improving-neighbor"


def test_likelihood_search_convergence_reports_budget_stop_reason() -> None:
    decision = resolve_nucleotide_likelihood_search_convergence_decision(
        best_improving_delta=1.0,
        improvement_tolerance=1e-9,
        budget_stopping_reason="candidate-budget-exhausted",
    )

    assert decision.should_stop is True
    assert decision.stopping_reason == "candidate-budget-exhausted"


def test_likelihood_search_convergence_reports_improvement_within_tolerance() -> None:
    decision = resolve_nucleotide_likelihood_search_convergence_decision(
        best_improving_delta=0.25,
        improvement_tolerance=0.5,
    )

    assert decision.should_stop is True
    assert decision.stopping_reason == "improvement-within-tolerance"


def test_likelihood_search_convergence_reports_repeated_topology() -> None:
    decision = resolve_nucleotide_likelihood_search_convergence_decision(
        best_improving_delta=1.0,
        improvement_tolerance=1e-9,
        candidate_topology_fingerprint="((A,B),(C,D))",
        seen_topology_fingerprints={"((A,B),(C,D))"},
    )

    assert decision.should_stop is True
    assert decision.stopping_reason == "repeated-topology"


def test_likelihood_search_convergence_reports_failure() -> None:
    decision = resolve_nucleotide_likelihood_search_convergence_decision(
        best_improving_delta=1.0,
        improvement_tolerance=1e-9,
        failure_detected=True,
    )

    assert decision.should_stop is True
    assert decision.stopping_reason == "search-failure"


def test_likelihood_search_convergence_allows_continuation_for_real_improvement() -> None:
    decision = resolve_nucleotide_likelihood_search_convergence_decision(
        best_improving_delta=1.0,
        improvement_tolerance=1e-9,
        candidate_topology_fingerprint="((A,B),(C,D))",
        seen_topology_fingerprints={"((A,C),(B,D))"},
    )

    assert decision.should_stop is False
    assert decision.stopping_reason is None

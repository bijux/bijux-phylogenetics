from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodSearchConvergenceDecision,
    resolve_nucleotide_likelihood_search_convergence_decision,
    validate_nucleotide_likelihood_search_improvement_tolerance,
)


def test_public_likelihood_exports_search_convergence_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodSearchConvergenceDecision
        is NucleotideLikelihoodSearchConvergenceDecision
    )
    assert (
        likelihood_api.resolve_nucleotide_likelihood_search_convergence_decision
        is resolve_nucleotide_likelihood_search_convergence_decision
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_search_improvement_tolerance
        is validate_nucleotide_likelihood_search_improvement_tolerance
    )

from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodOptimizationEquivalenceReport,
    check_nucleotide_likelihood_optimization_equivalence,
    check_nucleotide_likelihood_optimization_equivalence_from_alignment,
    validate_nucleotide_likelihood_optimization_equivalence_tolerances,
)


def test_public_runtime_exports_optimization_equivalence_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodOptimizationEquivalenceReport
        is NucleotideLikelihoodOptimizationEquivalenceReport
    )
    assert (
        likelihood_api.check_nucleotide_likelihood_optimization_equivalence
        is check_nucleotide_likelihood_optimization_equivalence
    )
    assert (
        likelihood_api.check_nucleotide_likelihood_optimization_equivalence_from_alignment
        is check_nucleotide_likelihood_optimization_equivalence_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_optimization_equivalence_tolerances
        is validate_nucleotide_likelihood_optimization_equivalence_tolerances
    )

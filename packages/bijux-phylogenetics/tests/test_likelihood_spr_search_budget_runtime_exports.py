from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodSprSearchBudget,
    search_nucleotide_likelihood_spr,
    search_nucleotide_likelihood_spr_from_alignment,
    validate_nucleotide_likelihood_spr_search_budget,
    write_nucleotide_likelihood_spr_artifacts,
)


def test_public_runtime_exports_spr_search_budget_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodSprSearchBudget
        is NucleotideLikelihoodSprSearchBudget
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_spr
        is search_nucleotide_likelihood_spr
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_spr_from_alignment
        is search_nucleotide_likelihood_spr_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_spr_search_budget
        is validate_nucleotide_likelihood_spr_search_budget
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_spr_artifacts
        is write_nucleotide_likelihood_spr_artifacts
    )

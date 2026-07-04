from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodNniCandidateRow,
    search_nucleotide_likelihood_nni,
    search_nucleotide_likelihood_nni_from_alignment,
    write_nucleotide_likelihood_nni_artifacts,
    write_nucleotide_likelihood_nni_candidate_table,
)


def test_public_runtime_exports_nni_best_improvement_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodNniCandidateRow
        is NucleotideLikelihoodNniCandidateRow
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_nni
        is search_nucleotide_likelihood_nni
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_nni_from_alignment
        is search_nucleotide_likelihood_nni_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_nni_artifacts
        is write_nucleotide_likelihood_nni_artifacts
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_nni_candidate_table
        is write_nucleotide_likelihood_nni_candidate_table
    )

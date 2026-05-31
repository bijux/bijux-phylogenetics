from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    NucleotideLikelihoodTbrSearchReport,
    NucleotideLikelihoodTbrTraceRow,
    search_nucleotide_likelihood_tbr,
    search_nucleotide_likelihood_tbr_from_alignment,
    write_nucleotide_likelihood_tbr_artifacts,
)


def test_public_likelihood_exports_tbr_search_surface() -> None:
    assert (
        likelihood_api.NucleotideLikelihoodTbrSearchReport
        is NucleotideLikelihoodTbrSearchReport
    )
    assert (
        likelihood_api.NucleotideLikelihoodTbrTraceRow
        is NucleotideLikelihoodTbrTraceRow
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_tbr
        is search_nucleotide_likelihood_tbr
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_tbr_from_alignment
        is search_nucleotide_likelihood_tbr_from_alignment
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_tbr_artifacts
        is write_nucleotide_likelihood_tbr_artifacts
    )

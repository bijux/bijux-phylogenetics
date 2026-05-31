from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    search_nucleotide_likelihood_nni,
    search_nucleotide_likelihood_nni_from_alignment,
    validate_nucleotide_likelihood_nni_improvement_policy,
    write_nucleotide_likelihood_nni_artifacts,
)


def test_public_likelihood_exports_nni_first_improvement_surface() -> None:
    assert (
        likelihood_api.search_nucleotide_likelihood_nni
        is search_nucleotide_likelihood_nni
    )
    assert (
        likelihood_api.search_nucleotide_likelihood_nni_from_alignment
        is search_nucleotide_likelihood_nni_from_alignment
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_nni_improvement_policy
        is validate_nucleotide_likelihood_nni_improvement_policy
    )
    assert (
        likelihood_api.write_nucleotide_likelihood_nni_artifacts
        is write_nucleotide_likelihood_nni_artifacts
    )

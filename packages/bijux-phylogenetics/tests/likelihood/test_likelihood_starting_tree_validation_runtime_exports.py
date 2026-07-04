from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    validate_nucleotide_likelihood_starting_tree,
    validate_nucleotide_likelihood_starting_tree_from_alignment,
)


def test_public_runtime_exports_likelihood_starting_tree_validation_surface() -> None:
    assert (
        likelihood_api.validate_nucleotide_likelihood_starting_tree
        is validate_nucleotide_likelihood_starting_tree
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_starting_tree_from_alignment
        is validate_nucleotide_likelihood_starting_tree_from_alignment
    )

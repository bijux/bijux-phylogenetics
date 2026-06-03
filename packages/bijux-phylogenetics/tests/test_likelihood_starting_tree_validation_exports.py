from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    validate_nucleotide_likelihood_starting_tree,
    validate_nucleotide_likelihood_starting_tree_from_alignment,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_validation import (
    validate_nucleotide_likelihood_starting_tree as validate_nucleotide_likelihood_starting_tree_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_validation import (
    validate_nucleotide_likelihood_starting_tree_from_alignment as validate_nucleotide_likelihood_starting_tree_from_alignment_impl,
)


def test_likelihood_exports_starting_tree_validation_surface() -> None:
    assert (
        validate_nucleotide_likelihood_starting_tree
        is validate_nucleotide_likelihood_starting_tree_impl
    )
    assert (
        validate_nucleotide_likelihood_starting_tree_from_alignment
        is validate_nucleotide_likelihood_starting_tree_from_alignment_impl
    )

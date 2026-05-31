from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    select_nucleotide_likelihood_starting_tree_pool,
    validate_nucleotide_likelihood_starting_tree_selection_count,
    validate_nucleotide_likelihood_starting_tree_selection_policy,
    validate_nucleotide_likelihood_starting_tree_strategy_priority,
)


def test_public_runtime_exports_likelihood_starting_tree_selection_surface() -> None:
    assert (
        likelihood_api.select_nucleotide_likelihood_starting_tree_pool
        is select_nucleotide_likelihood_starting_tree_pool
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_starting_tree_selection_policy
        is validate_nucleotide_likelihood_starting_tree_selection_policy
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_starting_tree_selection_count
        is validate_nucleotide_likelihood_starting_tree_selection_count
    )
    assert (
        likelihood_api.validate_nucleotide_likelihood_starting_tree_strategy_priority
        is validate_nucleotide_likelihood_starting_tree_strategy_priority
    )

from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    select_nucleotide_likelihood_starting_tree_pool,
    validate_nucleotide_likelihood_starting_tree_selection_count,
    validate_nucleotide_likelihood_starting_tree_selection_policy,
    validate_nucleotide_likelihood_starting_tree_strategy_priority,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_selection import (
    select_nucleotide_likelihood_starting_tree_pool as select_nucleotide_likelihood_starting_tree_pool_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_selection import (
    validate_nucleotide_likelihood_starting_tree_selection_count as validate_nucleotide_likelihood_starting_tree_selection_count_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_selection import (
    validate_nucleotide_likelihood_starting_tree_selection_policy as validate_nucleotide_likelihood_starting_tree_selection_policy_impl,
)
from bijux_phylogenetics.phylo.likelihood.starting_tree_selection import (
    validate_nucleotide_likelihood_starting_tree_strategy_priority as validate_nucleotide_likelihood_starting_tree_strategy_priority_impl,
)


def test_likelihood_exports_starting_tree_selection_surface() -> None:
    assert (
        select_nucleotide_likelihood_starting_tree_pool
        is select_nucleotide_likelihood_starting_tree_pool_impl
    )
    assert (
        validate_nucleotide_likelihood_starting_tree_selection_policy
        is validate_nucleotide_likelihood_starting_tree_selection_policy_impl
    )
    assert (
        validate_nucleotide_likelihood_starting_tree_selection_count
        is validate_nucleotide_likelihood_starting_tree_selection_count_impl
    )
    assert (
        validate_nucleotide_likelihood_starting_tree_strategy_priority
        is validate_nucleotide_likelihood_starting_tree_strategy_priority_impl
    )

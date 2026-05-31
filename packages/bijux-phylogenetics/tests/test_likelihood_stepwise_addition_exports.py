from __future__ import annotations

from bijux_phylogenetics.phylo.likelihood import (
    build_likelihood_stepwise_addition_tree,
    build_likelihood_stepwise_addition_tree_from_alignment,
    validate_likelihood_stepwise_addition_model,
)
from bijux_phylogenetics.phylo.likelihood.stepwise_addition import (
    build_likelihood_stepwise_addition_tree as build_likelihood_stepwise_addition_tree_impl,
)
from bijux_phylogenetics.phylo.likelihood.stepwise_addition import (
    build_likelihood_stepwise_addition_tree_from_alignment as build_likelihood_stepwise_addition_tree_from_alignment_impl,
)
from bijux_phylogenetics.phylo.likelihood.stepwise_addition import (
    validate_likelihood_stepwise_addition_model as validate_likelihood_stepwise_addition_model_impl,
)


def test_likelihood_exports_stepwise_addition_surface() -> None:
    assert (
        build_likelihood_stepwise_addition_tree
        is build_likelihood_stepwise_addition_tree_impl
    )
    assert (
        build_likelihood_stepwise_addition_tree_from_alignment
        is build_likelihood_stepwise_addition_tree_from_alignment_impl
    )
    assert (
        validate_likelihood_stepwise_addition_model
        is validate_likelihood_stepwise_addition_model_impl
    )

from __future__ import annotations

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    build_likelihood_stepwise_addition_tree,
    build_likelihood_stepwise_addition_tree_from_alignment,
    validate_likelihood_stepwise_addition_model,
)


def test_public_runtime_exports_likelihood_stepwise_addition_surface() -> None:
    assert (
        likelihood_api.build_likelihood_stepwise_addition_tree
        is build_likelihood_stepwise_addition_tree
    )
    assert (
        likelihood_api.build_likelihood_stepwise_addition_tree_from_alignment
        is build_likelihood_stepwise_addition_tree_from_alignment
    )
    assert (
        likelihood_api.validate_likelihood_stepwise_addition_model
        is validate_likelihood_stepwise_addition_model
    )

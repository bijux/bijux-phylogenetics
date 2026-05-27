from __future__ import annotations

import numpy

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.phylo.likelihood.pruning import (
    postorder_conditional_likelihoods,
    transition_probability_matrix,
)

from ..policy import normalize_array
from .likelihood_math import branch_length


def estimate_marginal_state_probabilities(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(cached_branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(cached_branch_length)
        if cached is None:
            cached = transition_probability_matrix(rate_matrix, cached_branch_length)
            transition_cache[cached_branch_length] = cached
        return cached

    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(state_order),
        leaf_likelihood=lambda node: _leaf_likelihood_vector(
            states_by_taxon,
            state_index=state_index,
            state_count=len(state_order),
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: transition(branch_length(child)),
    )
    posterior_by_node = {
        node_signature(node): pruning_pass.conditional_for_node(node)
        for node in tree.iter_nodes(order="preorder")
    }
    root_signature = node_signature(tree.root)
    posterior_by_node[root_signature] = normalize_array(
        root_prior * posterior_by_node[root_signature]
    )

    def preorder(node) -> None:
        parent_signature = node_signature(node)
        if node.is_leaf():
            return
        parent_probabilities = posterior_by_node[parent_signature]
        for child in node.children:
            if child.is_leaf():
                continue
            child_signature = node_signature(child)
            branch_transition = transition(branch_length(child))
            child_probabilities = posterior_by_node[child_signature]
            denominator = child_probabilities @ branch_transition
            updated = (parent_probabilities / denominator) @ branch_transition
            posterior_by_node[child_signature] = normalize_array(
                updated * child_probabilities
            )
            preorder(child)

    preorder(tree.root)
    return {
        node: {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                normalize_array(probabilities),
                strict=True,
            )
        }
        for node, probabilities in posterior_by_node.items()
    }


def _leaf_likelihood_vector(
    states_by_taxon: dict[str, str],
    *,
    state_index: dict[str, int],
    state_count: int,
    node_name: str | None,
) -> numpy.ndarray:
    if node_name is None:
        raise ValueError("leaf nodes in a discrete likelihood tree must be named")
    likelihood = numpy.zeros(state_count, dtype=float)
    likelihood[state_index[states_by_taxon[node_name]]] = 1.0
    return likelihood

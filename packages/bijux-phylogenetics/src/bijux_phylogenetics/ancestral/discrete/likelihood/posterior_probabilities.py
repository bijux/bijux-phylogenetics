from __future__ import annotations

import numpy

from bijux_phylogenetics.ancestral.common import node_signature

from ..policy import normalize_array
from .likelihood_math import branch_length
from .likelihood_math import transition_probability_matrix


def estimate_marginal_state_probabilities(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    state_index = {state: index for index, state in enumerate(state_order)}
    posterior_by_node: dict[str, numpy.ndarray] = {}
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(cached_branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(cached_branch_length)
        if cached is None:
            cached = transition_probability_matrix(rate_matrix, cached_branch_length)
            transition_cache[cached_branch_length] = cached
        return cached

    def postorder(node) -> numpy.ndarray:
        signature = node_signature(node)
        if node.is_leaf():
            partial = numpy.zeros(len(state_order), dtype=float)
            partial[state_index[states_by_taxon[node.name]]] = 1.0
            posterior_by_node[signature] = partial
            return partial
        partial = numpy.ones(len(state_order), dtype=float)
        for child in node.children:
            child_partial = postorder(child)
            partial *= transition(branch_length(child)) @ child_partial
        partial = normalize_array(partial)
        posterior_by_node[signature] = partial
        return partial

    postorder(tree.root)
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

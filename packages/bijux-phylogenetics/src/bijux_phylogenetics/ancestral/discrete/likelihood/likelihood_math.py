from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
    transition_probability_matrix,
)


def tree_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
) -> float:
    transition_cache: dict[float, numpy.ndarray] = {}

    def transition(cached_branch_length: float) -> numpy.ndarray:
        cached = transition_cache.get(cached_branch_length)
        if cached is None:
            cached = transition_probability_matrix(rate_matrix, cached_branch_length)
            transition_cache[cached_branch_length] = cached
        return cached

    state_index = {state: index for index, state in enumerate(state_order)}
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
    root_partial = pruning_pass.conditional_for_node(tree.root)
    if root_prior_mode == "observed":
        observed_root_total = float(root_partial.sum())
        if observed_root_total <= 0.0:
            return float("-inf")
        root_scale = float((root_partial @ root_partial) / observed_root_total)
        if root_scale <= 0.0:
            return float("-inf")
        subtree_log_scale = pruning_pass.subtree_log_scaling_for_node(tree.root)
        if math.isinf(subtree_log_scale) and subtree_log_scale < 0.0:
            return float("-inf")
        return subtree_log_scale + math.log(root_scale)
    if root_prior is None:
        raise ValueError("root_prior is required unless root_prior_mode is 'observed'")
    return log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=root_prior,
    )


def branch_length(node) -> float:
    if node.branch_length is None:
        return 1.0
    return max(float(node.branch_length), 0.0)


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

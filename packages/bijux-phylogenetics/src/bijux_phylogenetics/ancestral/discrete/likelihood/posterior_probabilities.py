from __future__ import annotations

import numpy

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    postorder_conditional_likelihoods,
)

from .likelihood_math import branch_length


def estimate_marginal_state_probabilities(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
) -> dict[str, dict[str, float]]:
    """Estimate governed internal-node state probabilities for one discrete Mk fit.

    This follows the `ape::ace(..., type='discrete')` `lik.anc` contract that the
    repository governs against for discrete ancestral reconstruction parity.
    """
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(state_order),
        leaf_likelihood=lambda node: _leaf_likelihood_vector(
            states_by_taxon,
            state_index=state_index,
            state_count=len(state_order),
            node_name=node.name,
        ),
        transition_matrix_for_child=lambda child: (
            transition_evaluator.transition_probability_matrix(branch_length(child))
        ),
        normalize_internal_conditionals=True,
    )
    root_node = tree.root
    if root_node.node_id is None:
        raise ValueError("tree root is missing a stable node_id")
    internal_nodes = [
        current_node
        for current_node in tree.iter_nodes(order="postorder")
        if not current_node.is_leaf()
    ]
    probability_by_node_id = {
        current_node.node_id: pruning_pass.conditional_for_node(current_node).copy()
        for current_node in internal_nodes
        if current_node.node_id is not None
    }
    probability_by_node_id[root_node.node_id] = _normalize_probability_row(
        root_prior * probability_by_node_id[root_node.node_id]
    )

    # Match ape::ace discrete `lik.anc` by applying the backward adjustment in
    # reverse postorder over internal branches after the normalized downpass.
    for current_node in reversed(internal_nodes):
        current_node_id = current_node.node_id
        if current_node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        parent_row = probability_by_node_id[current_node_id]
        for child in current_node.children:
            if child.is_leaf():
                continue
            child_id = child.node_id
            if child_id is None:
                raise ValueError("tree child node is missing a stable node_id")
            child_row = probability_by_node_id[child_id]
            branch_transition = transition_evaluator.transition_probability_matrix(
                branch_length(child)
            )
            denominator = child_row @ branch_transition
            if numpy.any(denominator <= 0.0) or not numpy.all(
                numpy.isfinite(denominator)
            ):
                raise ValueError(
                    "discrete ancestral backward adjustment encountered a non-positive branch message"
                )
            adjustment = parent_row / denominator
            probability_by_node_id[child_id] = (adjustment @ branch_transition) * (
                child_row
            )
        for node_id, probability_row in probability_by_node_id.items():
            probability_by_node_id[node_id] = _normalize_probability_row(
                probability_row
            )

    return {
        node_signature(current_node): {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                (
                    probability_by_node_id[current_node.node_id]
                    if not current_node.is_leaf()
                    else _leaf_likelihood_vector(
                        states_by_taxon,
                        state_index=state_index,
                        state_count=len(state_order),
                        node_name=current_node.name,
                    )
                ),
                strict=True,
            )
        }
        for current_node in tree.iter_nodes(order="preorder")
        if current_node.node_id is not None
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


def _normalize_probability_row(probability_row: numpy.ndarray) -> numpy.ndarray:
    candidate = numpy.asarray(probability_row, dtype=float)
    if not numpy.all(numpy.isfinite(candidate)):
        raise ValueError(
            "discrete ancestral backward adjustment produced non-finite probabilities"
        )
    if numpy.any(candidate < 0.0):
        raise ValueError(
            "discrete ancestral backward adjustment produced negative probabilities"
        )
    total = float(candidate.sum())
    if total <= 0.0:
        raise ValueError(
            "discrete ancestral backward adjustment produced a zero-sum probability row"
        )
    return candidate / total

from __future__ import annotations

import numpy

from bijux_phylogenetics.ancestral.common import node_signature
from bijux_phylogenetics.phylo.likelihood.posteriors import (
    compute_marginal_state_posteriors,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
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
    state_index = {state: index for index, state in enumerate(state_order)}
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)

    posterior_pass = compute_marginal_state_posteriors(
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
        root_prior=root_prior,
    )
    return {
        node_signature(current_node): {
            state: float(format(probability, ".15g"))
            for state, probability in zip(
                state_order,
                posterior_pass.posterior_for_node(current_node),
                strict=True,
            )
        }
        for current_node in tree.iter_nodes(order="preorder")
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

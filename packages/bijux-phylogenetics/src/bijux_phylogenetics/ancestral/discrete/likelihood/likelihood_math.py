from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood.discrete_observation_policies import (
    resolve_discrete_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.logspace import logsumexp
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    transition_probability_matrix as finite_state_transition_probability_matrix,
)


def transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    return finite_state_transition_probability_matrix(rate_matrix, branch_length)


def tree_log_likelihood(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
    ascertainment_policy: str = "none",
    observation_policy: str = "reject",
) -> float:
    base_log_likelihood = _tree_log_likelihood_without_ascertainment(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode=root_prior_mode,
        observation_policy=observation_policy,
    )
    if ascertainment_policy == "none":
        return base_log_likelihood
    if ascertainment_policy != "lewis-variable-only":
        raise ValueError(
            "unsupported discrete likelihood ascertainment policy: "
            f"{ascertainment_policy}"
        )
    conditioning_log_probability = variable_pattern_log_probability(
        tree,
        taxa=sorted(states_by_taxon),
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode=root_prior_mode,
    )
    if conditioning_log_probability == float("-inf"):
        return float("-inf")
    return base_log_likelihood - conditioning_log_probability


def variable_pattern_log_probability(
    tree,
    *,
    taxa: list[str],
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
) -> float:
    invariant_log_probability = invariant_pattern_log_probability(
        tree,
        taxa=taxa,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        root_prior_mode=root_prior_mode,
    )
    if invariant_log_probability == float("-inf"):
        return 0.0
    invariant_probability = math.exp(invariant_log_probability)
    if invariant_probability >= 1.0:
        return float("-inf")
    return math.log1p(-invariant_probability)


def invariant_pattern_log_probability(
    tree,
    *,
    taxa: list[str],
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
) -> float:
    return logsumexp(
        [
            _tree_log_likelihood_without_ascertainment(
                tree,
                dict.fromkeys(taxa, state),
                state_order=state_order,
                rate_matrix=rate_matrix,
                root_prior=root_prior,
                root_prior_mode=root_prior_mode,
            )
            for state in state_order
        ]
    )


def _tree_log_likelihood_without_ascertainment(
    tree,
    states_by_taxon: dict[str, str],
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    root_prior_mode: str = "given",
    observation_policy: str = "reject",
) -> float:
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(state_order),
        leaf_likelihood=lambda node: _leaf_likelihood_vector(
            states_by_taxon,
            state_order=state_order,
            node_name=node.name,
            observation_policy=observation_policy,
        ),
        transition_matrix_for_child=lambda child: (
            transition_evaluator.transition_probability_matrix(branch_length(child))
        ),
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
    state_order: list[str],
    node_name: str | None,
    observation_policy: str,
) -> numpy.ndarray:
    if node_name is None:
        raise ValueError("leaf nodes in a discrete likelihood tree must be named")
    return resolve_discrete_observation_leaf_vector(
        states_by_taxon[node_name],
        model_name="discrete Mk",
        node_name=node_name,
        state_order=state_order,
        observation_policy=observation_policy,
    )

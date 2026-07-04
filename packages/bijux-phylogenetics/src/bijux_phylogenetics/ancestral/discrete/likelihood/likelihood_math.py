from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CachedFiniteStateSiteEvaluation,
    CompiledFiniteStateLikelihoodContext,
    FiniteStateSiteEvaluationKey,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_site_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.logspace import logsumexp
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
    context = _compile_discrete_mk_context(
        tree,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        observation_policy="reject",
    )
    return logsumexp(
        [
            _site_log_likelihood_from_evaluation(
                tree,
                _evaluate_discrete_site_pattern(
                    context,
                    states_by_taxon=dict.fromkeys(taxa, state),
                ),
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
    context = _compile_discrete_mk_context(
        tree,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        observation_policy=observation_policy,
    )
    evaluation = _evaluate_discrete_site_pattern(
        context,
        states_by_taxon=states_by_taxon,
    )
    return _site_log_likelihood_from_evaluation(
        tree,
        evaluation,
        root_prior_mode=root_prior_mode,
    )


def branch_length(node) -> float:
    if node.branch_length is None:
        return 1.0
    return max(float(node.branch_length), 0.0)


def _compile_discrete_mk_context(
    tree,
    *,
    state_order: list[str],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | None,
    observation_policy: str,
) -> CompiledFiniteStateLikelihoodContext:
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior_policy=("provided" if root_prior is not None else "equal"),
        root_prior=root_prior,
        observation_policy=observation_policy,
        model_name="discrete Mk",
    )
    return compile_finite_state_likelihood_context(tree, surface)


def _evaluate_discrete_site_pattern(
    context: CompiledFiniteStateLikelihoodContext,
    *,
    states_by_taxon: dict[str, str],
) -> CachedFiniteStateSiteEvaluation:
    taxon_order = tuple(sorted(states_by_taxon))
    states = tuple(states_by_taxon[taxon] for taxon in taxon_order)
    evaluate_finite_state_site_log_likelihood(
        context,
        states,
        taxon_order=list(taxon_order),
    )
    evaluation = context.site_evaluation_by_key.get(
        FiniteStateSiteEvaluationKey(
            taxon_order=taxon_order,
            states=states,
        )
    )
    if evaluation is None:
        raise RuntimeError(
            "finite-state site evaluation cache did not retain the requested discrete Mk pattern"
        )
    return evaluation


def _site_log_likelihood_from_evaluation(
    tree,
    evaluation: CachedFiniteStateSiteEvaluation,
    *,
    root_prior_mode: str,
) -> float:
    if root_prior_mode == "observed":
        root_partial = evaluation.pruning_pass.conditional_for_node(tree.root)
        observed_root_total = float(root_partial.sum())
        if observed_root_total <= 0.0:
            return float("-inf")
        root_scale = float((root_partial @ root_partial) / observed_root_total)
        if root_scale <= 0.0:
            return float("-inf")
        subtree_log_scale = evaluation.pruning_pass.subtree_log_scaling_for_node(
            tree.root
        )
        if math.isinf(subtree_log_scale) and subtree_log_scale < 0.0:
            return float("-inf")
        return subtree_log_scale + math.log(root_scale)
    return evaluation.log_likelihood

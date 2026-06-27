from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.likelihood.discrete_observation_policies import (
    resolve_discrete_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def build_discrete_mk_finite_state_surface(
    *,
    state_order: list[str] | tuple[str, ...],
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    observation_policy: str = "reject",
    model_name: str = "discrete Mk",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state surface for one discrete Mk rate matrix."""
    validated_state_order = tuple(state_order)
    validated_root_prior = _validate_root_prior(
        root_prior,
        state_count=len(validated_state_order),
        model_name=model_name,
    )
    state_index = {
        state: index for index, state in enumerate(validated_state_order)
    }
    return FiniteStateLikelihoodSurface(
        model_name=model_name,
        state_labels=validated_state_order,
        root_prior=validated_root_prior,
        root_prior_source="provided",
        observation_policy=observation_policy,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_discrete_observation_leaf_vector(
                observed_state,
                model_name=model_name,
                node_name=taxon_name,
                state_order=list(validated_state_order),
                observation_policy=observation_policy,
            )
        ),
        parameter_values={"state_count": float(len(state_index))},
    )


def _validate_root_prior(
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    *,
    state_count: int,
    model_name: str,
) -> numpy.ndarray:
    vector = numpy.asarray(root_prior, dtype=float)
    if vector.shape != (state_count,):
        raise ValueError(
            f"{model_name} finite-state surface requires {state_count} root-prior entries"
        )
    if not numpy.all(numpy.isfinite(vector)):
        raise ValueError(f"{model_name} finite-state surface requires finite root prior")
    if numpy.any(vector < 0.0):
        raise ValueError(
            f"{model_name} finite-state surface requires nonnegative root prior"
        )
    total = float(vector.sum())
    if total <= 0.0:
        raise ValueError(
            f"{model_name} finite-state surface requires positive root-prior mass"
        )
    return vector / total

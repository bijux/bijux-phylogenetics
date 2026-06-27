from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    dna_observation_state_order,
    resolve_dna_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
)


def build_selected_nucleotide_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    owner_name: str,
    observation_policy: str = "reject",
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state DNA likelihood surface."""
    specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name=model_name,
        owner_name=owner_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return FiniteStateLikelihoodSurface(
        model_name=specification.model_name,
        state_labels=tuple(
            dna_observation_state_order(
                observation_policy=specification.observation_policy
            )
        ),
        root_prior=specification.root_prior,
        root_prior_source=specification.root_prior_source,
        observation_policy=specification.observation_policy,
        transition_matrix_for_branch_length=(
            specification.transition_matrix_for_branch_length
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_dna_observation_leaf_vector(
                {taxon_name: observed_state},
                model_name=specification.model_name,
                node_name=taxon_name,
                observation_policy=specification.observation_policy,
            )
        ),
        parameter_values=specification.parameter_values,
    )

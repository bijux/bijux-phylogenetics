from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.codon_observation_policies import (
    resolve_codon_observation_leaf_vector,
    validate_codon_observation_policy,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix

from .codon import normalize_codon_likelihood_records


def build_codon_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state codon CTMC likelihood surface."""
    normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    state_space, rate_matrix, frequencies, frequency_source = (
        build_equal_rate_codon_ctmc_rate_matrix(
            codon_frequencies,
            genetic_code=genetic_code,
        )
    )
    validated_observation_policy = validate_codon_observation_policy(
        observation_policy,
        owner_name="codon CTMC likelihood",
    )
    return FiniteStateLikelihoodSurface(
        model_name="codon CTMC",
        state_labels=tuple(state_space.state_order),
        root_prior=frequencies,
        root_prior_source=frequency_source,
        observation_policy=validated_observation_policy,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_codon_observation_leaf_vector(
                observed_state,
                model_name="codon CTMC",
                node_name=taxon_name,
                observation_policy=validated_observation_policy,
                state_order=state_space.state_order,
                state_index=state_space.state_index,
            )
        ),
    )

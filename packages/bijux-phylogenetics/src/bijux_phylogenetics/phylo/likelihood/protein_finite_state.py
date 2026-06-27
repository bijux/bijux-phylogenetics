from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
)
from bijux_phylogenetics.phylo.likelihood.poisson import (
    protein_poisson_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    UNIFORM_PROTEIN_ROOT_PRIOR,
    normalize_protein_likelihood_records,
    protein_leaf_likelihood_vector,
    validate_empirical_protein_rate_matrix,
    validate_protein_root_prior,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix


def build_protein_poisson_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state protein Poisson likelihood surface."""
    normalize_protein_likelihood_records(
        records,
        model_name="protein Poisson",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    rate_matrix = protein_poisson_rate_matrix()
    return FiniteStateLikelihoodSurface(
        model_name="protein Poisson",
        state_labels=PROTEIN_STATE_ORDER,
        root_prior=UNIFORM_PROTEIN_ROOT_PRIOR,
        root_prior_source="uniform",
        observation_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            protein_leaf_likelihood_vector(
                {taxon_name: observed_state},
                model_name="protein Poisson",
                node_name=taxon_name,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
                ambiguity_policy=ambiguity_policy,
            )
        ),
    )


def build_empirical_protein_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state empirical protein likelihood surface."""
    normalize_protein_likelihood_records(
        records,
        model_name="empirical protein matrix",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix",
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix",
        )
        root_prior_source = "provided"
    return FiniteStateLikelihoodSurface(
        model_name=f"empirical protein matrix ({matrix_label})",
        state_labels=PROTEIN_STATE_ORDER,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        observation_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(validated_rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            protein_leaf_likelihood_vector(
                {taxon_name: observed_state},
                model_name="empirical protein matrix",
                node_name=taxon_name,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
                ambiguity_policy=ambiguity_policy,
            )
        ),
    )


def _protein_observation_policy_label(
    *,
    gap_policy: str,
    missing_policy: str,
    ambiguity_policy: str,
) -> str:
    return (
        f"gap={gap_policy};missing={missing_policy};ambiguity={ambiguity_policy}"
    )

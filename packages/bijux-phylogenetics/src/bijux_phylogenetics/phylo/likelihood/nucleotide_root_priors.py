from __future__ import annotations

from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    UNIFORM_DNA_ROOT_PRIOR,
    estimate_empirical_dna_base_frequencies,
    fixed_state_dna_root_prior,
    validate_dna_root_prior,
)

_SUPPORTED_NUCLEOTIDE_ROOT_PRIOR_POLICIES = frozenset(
    {"equal", "empirical", "stationary", "fixed-state", "provided", "user-supplied"}
)


@dataclass(frozen=True, slots=True)
class ResolvedNucleotideRootPrior:
    """One validated nucleotide root-prior vector plus its durable source label."""

    root_prior: numpy.ndarray
    root_prior_source: str


def validate_nucleotide_root_prior_policy(
    root_prior_policy: str,
    *,
    owner_name: str,
) -> str:
    normalized_policy = root_prior_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_NUCLEOTIDE_ROOT_PRIOR_POLICIES:
        raise ValueError(
            f"{owner_name} root_prior_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_NUCLEOTIDE_ROOT_PRIOR_POLICIES))
        )
    if normalized_policy == "user-supplied":
        return "provided"
    return normalized_policy


def resolve_nucleotide_root_prior(
    records: list[AlignmentRecord],
    *,
    owner_name: str,
    default_policy: str,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    stationary_frequencies: numpy.ndarray | None = None,
) -> ResolvedNucleotideRootPrior:
    normalized_policy = (
        validate_nucleotide_root_prior_policy(
            root_prior_policy,
            owner_name=owner_name,
        )
        if root_prior_policy is not None
        else validate_nucleotide_root_prior_policy(
            default_policy,
            owner_name=owner_name,
        )
    )
    if normalized_policy == "equal":
        _reject_extra_root_prior_inputs(
            owner_name,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideRootPrior(
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
            root_prior_source="equal",
        )
    if normalized_policy == "empirical":
        _reject_extra_root_prior_inputs(
            owner_name,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideRootPrior(
            root_prior=estimate_empirical_dna_base_frequencies(records),
            root_prior_source="empirical",
        )
    if normalized_policy == "stationary":
        _reject_extra_root_prior_inputs(
            owner_name,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        if stationary_frequencies is None:
            raise ValueError(
                f"{owner_name} does not provide stationary frequencies for root_prior_policy 'stationary'"
            )
        return ResolvedNucleotideRootPrior(
            root_prior=validate_dna_root_prior(
                stationary_frequencies,
                owner_name=owner_name,
            ),
            root_prior_source="stationary",
        )
    if normalized_policy == "fixed-state":
        if root_prior is not None:
            raise ValueError(
                f"{owner_name} does not accept 'root_prior' when root_prior_policy is 'fixed-state'"
            )
        if fixed_root_state is None:
            raise ValueError(
                f"{owner_name} requires 'fixed_root_state' when root_prior_policy is 'fixed-state'"
            )
        return ResolvedNucleotideRootPrior(
            root_prior=fixed_state_dna_root_prior(
                fixed_root_state,
                owner_name=owner_name,
            ),
            root_prior_source="fixed-state",
        )
    if fixed_root_state is not None:
        raise ValueError(
            f"{owner_name} does not accept 'fixed_root_state' when root_prior_policy is 'provided'"
        )
    if root_prior is None:
        raise ValueError(
            f"{owner_name} requires 'root_prior' when root_prior_policy is 'provided'"
        )
    return ResolvedNucleotideRootPrior(
        root_prior=validate_dna_root_prior(
            root_prior,
            owner_name=owner_name,
        ),
        root_prior_source="provided",
    )


def _reject_extra_root_prior_inputs(
    owner_name: str,
    *,
    root_prior: object,
    fixed_root_state: str | None,
) -> None:
    if root_prior is not None:
        raise ValueError(
            f"{owner_name} does not accept 'root_prior' for this root_prior_policy"
        )
    if fixed_root_state is not None:
        raise ValueError(
            f"{owner_name} does not accept 'fixed_root_state' for this root_prior_policy"
        )

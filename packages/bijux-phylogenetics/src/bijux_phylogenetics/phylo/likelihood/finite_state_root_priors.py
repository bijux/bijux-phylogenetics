from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

_SUPPORTED_FINITE_STATE_ROOT_PRIOR_POLICIES = frozenset(
    {"equal", "empirical", "stationary", "fixed-state", "provided", "user-supplied"}
)


@dataclass(frozen=True, slots=True)
class ResolvedFiniteStateRootPrior:
    """One validated finite-state root-prior vector plus its durable source label."""

    root_prior: numpy.ndarray
    root_prior_source: str


def validate_finite_state_root_prior_policy(
    root_prior_policy: str,
    *,
    owner_name: str,
) -> str:
    normalized_policy = root_prior_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_FINITE_STATE_ROOT_PRIOR_POLICIES:
        raise ValueError(
            f"{owner_name} root_prior_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_FINITE_STATE_ROOT_PRIOR_POLICIES))
        )
    if normalized_policy == "user-supplied":
        return "provided"
    return normalized_policy


def validate_finite_state_root_prior(
    root_prior: Mapping[str, float] | Sequence[float] | numpy.ndarray,
    *,
    state_labels: Sequence[str],
    owner_name: str,
) -> numpy.ndarray:
    resolved_state_labels = tuple(state_labels)
    if isinstance(root_prior, Mapping):
        unexpected_states = sorted(
            state for state in root_prior if state not in set(resolved_state_labels)
        )
        missing_states = [
            state for state in resolved_state_labels if state not in root_prior
        ]
        if unexpected_states or missing_states:
            details: list[str] = []
            if unexpected_states:
                details.append(
                    "unexpected states: " + ", ".join(unexpected_states[:8])
                )
            if missing_states:
                details.append("missing states: " + ", ".join(missing_states[:8]))
            raise InvalidAlignmentError(
                f"{owner_name} requires probabilities for exactly the declared finite-state order"
                + (f" ({'; '.join(details)})" if details else "")
            )
        vector = numpy.asarray(
            [float(root_prior[state]) for state in resolved_state_labels],
            dtype=float,
        )
    else:
        vector = numpy.asarray(root_prior, dtype=float)
    return _normalize_finite_state_root_prior_vector(
        vector,
        state_labels=resolved_state_labels,
        owner_name=owner_name,
    )


def fixed_state_finite_state_root_prior(
    fixed_root_state: str,
    *,
    state_labels: Sequence[str],
    owner_name: str,
) -> numpy.ndarray:
    resolved_state_labels = tuple(state_labels)
    state_index = {state: index for index, state in enumerate(resolved_state_labels)}
    if fixed_root_state not in state_index:
        raise InvalidAlignmentError(
            f"{owner_name} fixed_root_state must be one of {list(resolved_state_labels)}"
        )
    vector = numpy.zeros(len(resolved_state_labels), dtype=float)
    vector[state_index[fixed_root_state]] = 1.0
    return vector


def equal_finite_state_root_prior(
    state_labels: Sequence[str],
) -> numpy.ndarray:
    resolved_state_labels = tuple(state_labels)
    return numpy.full(
        len(resolved_state_labels),
        1.0 / len(resolved_state_labels),
        dtype=float,
    )


def estimate_empirical_finite_state_root_prior(
    records: list[AlignmentRecord],
    *,
    state_labels: Sequence[str],
    owner_name: str,
    sequence_tokenizer: Callable[[str], Iterable[str]],
) -> numpy.ndarray:
    resolved_state_labels = tuple(state_labels)
    state_index = {state: index for index, state in enumerate(resolved_state_labels)}
    counts = numpy.zeros(len(resolved_state_labels), dtype=float)
    for record in records:
        for token in sequence_tokenizer(record.sequence):
            index = state_index.get(token)
            if index is not None:
                counts[index] += 1.0
    total = float(counts.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            f"{owner_name} empirical root_prior_policy requires at least one observed unambiguous state"
        )
    return counts / total


def resolve_finite_state_root_prior(
    records: list[AlignmentRecord],
    *,
    state_labels: Sequence[str],
    owner_name: str,
    default_policy: str,
    root_prior_policy: str | None = None,
    root_prior: Mapping[str, float] | Sequence[float] | numpy.ndarray | None = None,
    fixed_root_state: str | None = None,
    stationary_frequencies: Mapping[str, float]
    | Sequence[float]
    | numpy.ndarray
    | None = None,
    sequence_tokenizer: Callable[[str], Iterable[str]] | None = None,
) -> ResolvedFiniteStateRootPrior:
    normalized_policy = (
        validate_finite_state_root_prior_policy(
            root_prior_policy,
            owner_name=owner_name,
        )
        if root_prior_policy is not None
        else validate_finite_state_root_prior_policy(
            default_policy,
            owner_name=owner_name,
        )
    )
    resolved_state_labels = tuple(state_labels)
    if normalized_policy == "equal":
        _reject_extra_root_prior_inputs(
            owner_name,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedFiniteStateRootPrior(
            root_prior=equal_finite_state_root_prior(resolved_state_labels),
            root_prior_source="equal",
        )
    if normalized_policy == "empirical":
        _reject_extra_root_prior_inputs(
            owner_name,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        if sequence_tokenizer is None:
            raise ValueError(
                f"{owner_name} does not provide empirical state tokenization for root_prior_policy 'empirical'"
            )
        return ResolvedFiniteStateRootPrior(
            root_prior=estimate_empirical_finite_state_root_prior(
                records,
                state_labels=resolved_state_labels,
                owner_name=owner_name,
                sequence_tokenizer=sequence_tokenizer,
            ),
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
        return ResolvedFiniteStateRootPrior(
            root_prior=validate_finite_state_root_prior(
                stationary_frequencies,
                state_labels=resolved_state_labels,
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
        return ResolvedFiniteStateRootPrior(
            root_prior=fixed_state_finite_state_root_prior(
                fixed_root_state,
                state_labels=resolved_state_labels,
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
    return ResolvedFiniteStateRootPrior(
        root_prior=validate_finite_state_root_prior(
            root_prior,
            state_labels=resolved_state_labels,
            owner_name=owner_name,
        ),
        root_prior_source="provided",
    )


def _normalize_finite_state_root_prior_vector(
    vector: numpy.ndarray,
    *,
    state_labels: Sequence[str],
    owner_name: str,
) -> numpy.ndarray:
    resolved_state_labels = tuple(state_labels)
    if vector.shape != (len(resolved_state_labels),):
        raise InvalidAlignmentError(
            f"{owner_name} requires exactly {len(resolved_state_labels)} root-prior entries in the declared state order"
        )
    if not numpy.all(numpy.isfinite(vector)):
        raise InvalidAlignmentError(
            f"{owner_name} root prior must contain only finite values"
        )
    if numpy.any(vector < 0.0):
        raise InvalidAlignmentError(f"{owner_name} root prior must be nonnegative")
    total = float(vector.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            f"{owner_name} root prior must sum to a positive value"
        )
    return vector / total


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

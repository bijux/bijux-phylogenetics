from __future__ import annotations

from collections.abc import Mapping
import math

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_INDEX,
    DNA_STATE_ORDER,
    validate_dna_base_frequencies,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_root_priors import (
    ResolvedNucleotideRootPrior,
    resolve_nucleotide_root_prior,
)
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
    PhylogeneticsError,
)

_FIFTH_STATE_GAP_SYMBOL = "-"
_EXPLICIT_MISSING_SYMBOLS = frozenset({"?"})
_SUPPORTED_DNA_OBSERVATION_POLICIES = frozenset(
    {"reject", "treat-as-missing", "ambiguity-vector", "fifth-state"}
)
_DNA_AMBIGUITY_STATES_BY_SYMBOL: Mapping[str, tuple[str, ...]] = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("C", "G"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": DNA_STATE_ORDER,
}


def validate_dna_observation_policy(
    observation_policy: str,
    *,
    owner_name: str,
) -> str:
    normalized_policy = observation_policy.strip().lower()
    if normalized_policy not in _SUPPORTED_DNA_OBSERVATION_POLICIES:
        raise ValueError(
            f"{owner_name} observation_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_DNA_OBSERVATION_POLICIES))
        )
    return normalized_policy


def normalize_dna_likelihood_records(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    observation_policy: str = "reject",
) -> list[AlignmentRecord]:
    validated_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name=f"{model_name} likelihood",
    )
    normalized_records: list[AlignmentRecord] = []
    allowed_states = (
        set(DNA_STATE_ORDER)
        | set(_DNA_AMBIGUITY_STATES_BY_SYMBOL)
        | {_FIFTH_STATE_GAP_SYMBOL}
        | _EXPLICIT_MISSING_SYMBOLS
    )
    for record in records:
        normalized_sequence = record.sequence.upper()
        invalid_states = sorted(
            {state for state in normalized_sequence if state not in allowed_states}
        )
        if invalid_states:
            joined_states = ", ".join(invalid_states)
            raise InvalidAlignmentError(
                f"{model_name} likelihood does not recognize DNA observation states {joined_states}"
            )
        if validated_policy == "reject":
            rejected_states = sorted(
                {state for state in normalized_sequence if state not in DNA_STATE_INDEX}
            )
            if rejected_states:
                joined_states = ", ".join(rejected_states)
                raise InvalidAlignmentError(
                    f"{model_name} likelihood observation policy 'reject' requires A, C, G, and T only; "
                    f"record '{record.identifier}' contains {joined_states}"
                )
        normalized_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=normalized_sequence,
            )
        )
    return normalized_records


def dna_observation_state_order(
    *,
    observation_policy: str,
) -> tuple[str, ...]:
    validated_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name="DNA likelihood observation policy",
    )
    if validated_policy == "fifth-state":
        return (*DNA_STATE_ORDER, _FIFTH_STATE_GAP_SYMBOL)
    return DNA_STATE_ORDER


def resolve_dna_observation_leaf_vector(
    states_by_taxon: dict[str, str],
    *,
    model_name: str,
    node_name: str | None,
    observation_policy: str,
) -> numpy.ndarray:
    if node_name is None:
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires named tree tips for alignment lookup"
        )
    validated_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name=f"{model_name} likelihood",
    )
    state_order = dna_observation_state_order(observation_policy=validated_policy)
    symbol = states_by_taxon[node_name].upper()
    vector = numpy.zeros(len(state_order), dtype=float)
    if symbol in DNA_STATE_INDEX:
        vector[DNA_STATE_INDEX[symbol]] = 1.0
        return vector
    if symbol in _DNA_AMBIGUITY_STATES_BY_SYMBOL:
        if validated_policy == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood observation policy 'reject' does not allow ambiguity code '{symbol}'"
            )
        if validated_policy == "treat-as-missing":
            vector[:] = 1.0
            return vector
        for state in _DNA_AMBIGUITY_STATES_BY_SYMBOL[symbol]:
            vector[DNA_STATE_INDEX[state]] = 1.0
        return vector
    if symbol == _FIFTH_STATE_GAP_SYMBOL:
        if validated_policy == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood observation policy 'reject' does not allow gap state '-'"
            )
        if validated_policy == "fifth-state":
            vector[-1] = 1.0
            return vector
        vector[:] = 1.0
        return vector
    if symbol in _EXPLICIT_MISSING_SYMBOLS:
        if validated_policy == "reject":
            raise InvalidAlignmentError(
                f"{model_name} likelihood observation policy 'reject' does not allow missing-data state '{symbol}'"
            )
        vector[:] = 1.0
        return vector
    raise InvalidAlignmentError(
        f"{model_name} likelihood does not recognize DNA observation state '{symbol}'"
    )


def estimate_empirical_dna_base_frequencies_from_records(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    observation_policy: str,
) -> numpy.ndarray:
    validated_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name=f"{model_name} likelihood",
    )
    counts = numpy.zeros(len(DNA_STATE_ORDER), dtype=float)
    for record in records:
        for symbol in record.sequence.upper():
            if symbol in DNA_STATE_INDEX:
                counts[DNA_STATE_INDEX[symbol]] += 1.0
                continue
            if symbol in _DNA_AMBIGUITY_STATES_BY_SYMBOL:
                if validated_policy == "treat-as-missing":
                    continue
                ambiguity_states = _DNA_AMBIGUITY_STATES_BY_SYMBOL[symbol]
                fractional_weight = 1.0 / len(ambiguity_states)
                for state in ambiguity_states:
                    counts[DNA_STATE_INDEX[state]] += fractional_weight
            elif (
                symbol in _EXPLICIT_MISSING_SYMBOLS or symbol == _FIFTH_STATE_GAP_SYMBOL
            ):
                continue
            else:
                raise InvalidAlignmentError(
                    f"{model_name} likelihood does not recognize DNA observation state '{symbol}'"
                )
    total = float(counts.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} likelihood base-frequency estimation requires at least one resolved or partially resolved nucleotide"
        )
    return counts / total


def estimate_empirical_gap_state_frequency(
    records: list[AlignmentRecord],
    *,
    model_name: str,
) -> float:
    informative_symbol_count = 0.0
    gap_symbol_count = 0.0
    for record in records:
        for symbol in record.sequence.upper():
            if symbol in _EXPLICIT_MISSING_SYMBOLS:
                continue
            if symbol == _FIFTH_STATE_GAP_SYMBOL:
                informative_symbol_count += 1.0
                gap_symbol_count += 1.0
                continue
            if symbol in DNA_STATE_INDEX or symbol in _DNA_AMBIGUITY_STATES_BY_SYMBOL:
                informative_symbol_count += 1.0
                continue
            raise InvalidAlignmentError(
                f"{model_name} likelihood does not recognize DNA observation state '{symbol}'"
            )
    if informative_symbol_count <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} fifth-state observation policy requires at least one non-missing observation"
        )
    observed_gap_frequency = gap_symbol_count / informative_symbol_count
    return min(max(observed_gap_frequency, 1e-6), 1.0 - 1e-6)


def extend_dna_stationary_frequencies_with_gap_state(
    nucleotide_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...],
    *,
    gap_state_frequency: float,
    model_name: str,
) -> numpy.ndarray:
    validated_nucleotide_frequencies = validate_dna_base_frequencies(
        nucleotide_frequencies,
        model_name=model_name,
    )
    if not math.isfinite(gap_state_frequency) or not (0.0 < gap_state_frequency < 1.0):
        raise InvalidAlignmentError(
            f"{model_name} gap-state frequency must lie strictly between zero and one"
        )
    nucleotide_scale = 1.0 - gap_state_frequency
    return numpy.concatenate(
        (
            validated_nucleotide_frequencies * nucleotide_scale,
            numpy.array([gap_state_frequency], dtype=float),
        )
    )


def augment_dna_rate_matrix_with_gap_state(
    nucleotide_rate_matrix: numpy.ndarray,
    *,
    nucleotide_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...],
    gap_state_frequency: float,
    model_name: str,
    gap_exchangeability: float = 1.0,
) -> numpy.ndarray:
    validated_nucleotide_frequencies = validate_dna_base_frequencies(
        nucleotide_frequencies,
        model_name=model_name,
    )
    extended_stationary_frequencies = extend_dna_stationary_frequencies_with_gap_state(
        validated_nucleotide_frequencies,
        gap_state_frequency=gap_state_frequency,
        model_name=model_name,
    )
    candidate_rate_matrix = numpy.asarray(nucleotide_rate_matrix, dtype=float)
    if candidate_rate_matrix.shape != (len(DNA_STATE_ORDER), len(DNA_STATE_ORDER)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood rate matrix must be 4x4 before gap-state augmentation"
        )
    if not numpy.all(numpy.isfinite(candidate_rate_matrix)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood rate matrix must contain only finite values before gap-state augmentation"
        )
    if not math.isfinite(gap_exchangeability) or gap_exchangeability <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} likelihood gap-state exchangeability must be a finite positive value"
        )
    augmented_rate_matrix = numpy.zeros((5, 5), dtype=float)
    for left_index in range(len(DNA_STATE_ORDER)):
        for right_index in range(left_index + 1, len(DNA_STATE_ORDER)):
            exchangeability = 0.5 * (
                candidate_rate_matrix[left_index, right_index]
                / validated_nucleotide_frequencies[right_index]
                + candidate_rate_matrix[right_index, left_index]
                / validated_nucleotide_frequencies[left_index]
            )
            augmented_rate_matrix[left_index, right_index] = (
                exchangeability * extended_stationary_frequencies[right_index]
            )
            augmented_rate_matrix[right_index, left_index] = (
                exchangeability * extended_stationary_frequencies[left_index]
            )
    gap_index = len(DNA_STATE_ORDER)
    for nucleotide_index in range(len(DNA_STATE_ORDER)):
        augmented_rate_matrix[nucleotide_index, gap_index] = (
            gap_exchangeability * extended_stationary_frequencies[gap_index]
        )
        augmented_rate_matrix[gap_index, nucleotide_index] = (
            gap_exchangeability * extended_stationary_frequencies[nucleotide_index]
        )
    for row_index in range(augmented_rate_matrix.shape[0]):
        augmented_rate_matrix[row_index, row_index] = -float(
            numpy.sum(augmented_rate_matrix[row_index, :])
        )
    try:
        return normalize_ctmc_rate_matrix_by_expected_substitution_rate(
            augmented_rate_matrix,
            extended_stationary_frequencies,
            state_labels=dna_observation_state_order(observation_policy="fifth-state"),
        )
    except PhylogeneticsError as error:
        raise InvalidAlignmentError(
            f"{model_name} fifth-state rate augmentation requires a positive finite expected substitution rate"
        ) from error


def resolve_default_dna_root_prior_for_observation_policy(
    records: list[AlignmentRecord],
    *,
    owner_name: str,
    default_policy: str,
    root_prior_policy: str | None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None,
    fixed_root_state: str | None,
    stationary_frequencies: numpy.ndarray,
    observation_policy: str,
) -> ResolvedNucleotideRootPrior:
    validated_policy = validate_dna_observation_policy(
        observation_policy,
        owner_name=owner_name,
    )
    if validated_policy != "fifth-state":
        return resolve_nucleotide_root_prior(
            records,
            owner_name=owner_name,
            default_policy=default_policy,
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            stationary_frequencies=stationary_frequencies,
        )
    if (
        root_prior_policy is not None
        or root_prior is not None
        or fixed_root_state is not None
    ):
        raise ValueError(
            f"{owner_name} does not yet support explicit root-prior overrides with observation_policy 'fifth-state'"
        )
    gap_state_frequency = estimate_empirical_gap_state_frequency(
        records,
        model_name=owner_name,
    )
    if default_policy == "equal":
        return ResolvedNucleotideRootPrior(
            root_prior=numpy.full(5, 0.2, dtype=float),
            root_prior_source="equal",
        )
    if default_policy == "stationary":
        return ResolvedNucleotideRootPrior(
            root_prior=extend_dna_stationary_frequencies_with_gap_state(
                stationary_frequencies,
                gap_state_frequency=gap_state_frequency,
                model_name=owner_name,
            ),
            root_prior_source="stationary",
        )
    raise ValueError(
        f"{owner_name} does not support default root-prior policy '{default_policy}' with observation_policy 'fifth-state'"
    )


__all__ = [
    "augment_dna_rate_matrix_with_gap_state",
    "dna_observation_state_order",
    "estimate_empirical_dna_base_frequencies_from_records",
    "estimate_empirical_gap_state_frequency",
    "extend_dna_stationary_frequencies_with_gap_state",
    "normalize_dna_likelihood_records",
    "resolve_default_dna_root_prior_for_observation_policy",
    "resolve_dna_observation_leaf_vector",
    "validate_dna_observation_policy",
]

from __future__ import annotations

import math
from typing import cast

import numpy

from bijux_phylogenetics.core.categorical_probability_vectors import (
    build_categorical_probability_vector,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    normalize_ctmc_rate_matrix_by_expected_substitution_rate,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
    PhylogeneticsError,
)

DNA_STATE_ORDER = ("A", "C", "G", "T")
DNA_STATE_INDEX = {state: index for index, state in enumerate(DNA_STATE_ORDER)}
DNA_EXCHANGEABILITY_ORDER = (
    ("A", "C"),
    ("A", "G"),
    ("A", "T"),
    ("C", "G"),
    ("C", "T"),
    ("G", "T"),
)
DNA_TRANSITION_PAIRS = frozenset(
    {
        ("A", "G"),
        ("G", "A"),
        ("C", "T"),
        ("T", "C"),
    }
)
UNIFORM_DNA_ROOT_PRIOR = numpy.full(4, 0.25, dtype=float)


def normalize_unambiguous_dna_records(
    records: list[AlignmentRecord],
    *,
    model_name: str,
) -> list[AlignmentRecord]:
    normalized_records: list[AlignmentRecord] = []
    for record in records:
        normalized_sequence = record.sequence.upper()
        invalid_states = sorted(
            {
                state
                for state in normalized_sequence
                if state not in DNA_STATE_INDEX
            }
        )
        if invalid_states:
            joined_states = ", ".join(invalid_states)
            raise InvalidAlignmentError(
                f"{model_name} likelihood currently requires unambiguous DNA states A, C, G, and T only; "
                f"record '{record.identifier}' contains {joined_states}"
            )
        normalized_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=normalized_sequence,
            )
        )
    return normalized_records


def estimate_empirical_dna_base_frequencies(
    records: list[AlignmentRecord],
) -> numpy.ndarray:
    counts = numpy.zeros(len(DNA_STATE_ORDER), dtype=float)
    for record in records:
        for symbol in record.sequence.upper():
            counts[DNA_STATE_INDEX[symbol]] += 1.0
    total = float(counts.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            "DNA likelihood base-frequency estimation requires at least one resolved nucleotide"
        )
    return counts / total


def validate_dna_base_frequencies(
    base_frequencies: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...],
    *,
    model_name: str,
) -> numpy.ndarray:
    if isinstance(base_frequencies, dict):
        if set(base_frequencies) != set(DNA_STATE_ORDER):
            raise InvalidAlignmentError(
                f"{model_name} likelihood requires base frequencies for exactly A, C, G, and T"
            )
        vector = numpy.array(
            [float(base_frequencies[state]) for state in DNA_STATE_ORDER],
            dtype=float,
        )
    else:
        vector = numpy.asarray(base_frequencies, dtype=float)
    if vector.shape != (len(DNA_STATE_ORDER),):
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires exactly four base frequencies in A/C/G/T order"
        )
    try:
        validated_vector = build_categorical_probability_vector(
            dict(zip(DNA_STATE_ORDER, vector.tolist(), strict=True)),
            expected_states=DNA_STATE_ORDER,
        )
    except PhylogeneticsError as error:
        raise _dna_probability_error(
            owner_name=f"{model_name} likelihood",
            parameter_name="base frequencies",
            error=error,
        ) from error
    return numpy.array(validated_vector.probabilities, dtype=float)


def validate_dna_root_prior(
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...],
    *,
    owner_name: str,
) -> numpy.ndarray:
    if isinstance(root_prior, dict):
        if set(root_prior) != set(DNA_STATE_ORDER):
            raise InvalidAlignmentError(
                f"{owner_name} requires root-prior probabilities for exactly A, C, G, and T"
            )
        vector = numpy.array(
            [float(root_prior[state]) for state in DNA_STATE_ORDER],
            dtype=float,
        )
    else:
        vector = numpy.asarray(root_prior, dtype=float)
    if vector.shape != (len(DNA_STATE_ORDER),):
        raise InvalidAlignmentError(
            f"{owner_name} requires exactly four root-prior probabilities in A/C/G/T order"
        )
    try:
        validated_vector = build_categorical_probability_vector(
            dict(zip(DNA_STATE_ORDER, vector.tolist(), strict=True)),
            expected_states=DNA_STATE_ORDER,
        )
    except PhylogeneticsError as error:
        raise _dna_probability_error(
            owner_name=owner_name,
            parameter_name="root-prior probabilities",
            error=error,
        ) from error
    return numpy.array(validated_vector.probabilities, dtype=float)


def fixed_state_dna_root_prior(
    state: str,
    *,
    owner_name: str,
) -> numpy.ndarray:
    normalized_state = state.strip().upper()
    if normalized_state not in DNA_STATE_INDEX:
        raise InvalidAlignmentError(
            f"{owner_name} fixed root state must be one of A, C, G, and T"
        )
    root_prior = numpy.zeros(len(DNA_STATE_ORDER), dtype=float)
    root_prior[DNA_STATE_INDEX[normalized_state]] = 1.0
    return root_prior


def validate_dna_exchangeabilities(
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    *,
    model_name: str,
) -> numpy.ndarray:
    if isinstance(exchangeabilities, dict):
        if set(exchangeabilities) == set(DNA_EXCHANGEABILITY_ORDER):
            paired_exchangeabilities = cast(dict[tuple[str, str], float], exchangeabilities)
            ordered_exchangeabilities = [
                float(paired_exchangeabilities[pair])
                for pair in DNA_EXCHANGEABILITY_ORDER
            ]
            vector = numpy.asarray(ordered_exchangeabilities, dtype=float)
        elif set(exchangeabilities) == {
            "".join(pair) for pair in DNA_EXCHANGEABILITY_ORDER
        }:
            named_exchangeabilities = cast(dict[str, float], exchangeabilities)
            ordered_exchangeabilities = [
                float(named_exchangeabilities["".join(pair)])
                for pair in DNA_EXCHANGEABILITY_ORDER
            ]
            vector = numpy.asarray(ordered_exchangeabilities, dtype=float)
        else:
            raise InvalidAlignmentError(
                f"{model_name} likelihood requires exchangeabilities for exactly AC, AG, AT, CG, CT, and GT"
            )
    else:
        vector = numpy.asarray(exchangeabilities, dtype=float)
    if vector.shape != (len(DNA_EXCHANGEABILITY_ORDER),):
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires exactly six exchangeabilities in AC/AG/AT/CG/CT/GT order"
        )
    if not numpy.all(numpy.isfinite(vector)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood exchangeabilities must all be finite"
        )
    if numpy.any(vector < 0.0):
        raise InvalidAlignmentError(
            f"{model_name} likelihood exchangeabilities must be nonnegative"
        )
    if float(vector.sum()) <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} likelihood exchangeabilities must contain at least one positive value"
        )
    return vector.astype(float, copy=True)


def normalize_dna_exchangeabilities_by_anchor(
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    *,
    model_name: str,
    anchor_pair: tuple[str, str] = ("A", "C"),
) -> numpy.ndarray:
    vector = validate_dna_exchangeabilities(
        exchangeabilities,
        model_name=model_name,
    )
    try:
        anchor_index = DNA_EXCHANGEABILITY_ORDER.index(anchor_pair)
    except ValueError as error:
        raise ValueError(f"unsupported DNA exchangeability anchor {anchor_pair}") from error
    anchor_value = float(vector[anchor_index])
    if anchor_value <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires a positive {''.join(anchor_pair)} exchangeability anchor"
        )
    return vector / anchor_value


def normalize_dna_rate_matrix(
    off_diagonal_rates: numpy.ndarray,
    *,
    stationary_frequencies: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...],
    model_name: str,
) -> numpy.ndarray:
    candidate = numpy.asarray(off_diagonal_rates, dtype=float)
    if candidate.shape != (len(DNA_STATE_ORDER), len(DNA_STATE_ORDER)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood rate matrix must be 4x4 in A/C/G/T order"
        )
    if not numpy.all(numpy.isfinite(candidate)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood rate matrix must contain only finite values"
        )
    if numpy.any(candidate < 0.0):
        raise InvalidAlignmentError(
            f"{model_name} likelihood rate matrix must not contain negative off-diagonal rates"
        )
    rate_matrix = candidate.copy()
    numpy.fill_diagonal(rate_matrix, 0.0)
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    stationary = validate_dna_base_frequencies(
        stationary_frequencies,
        model_name=model_name,
    )
    try:
        return normalize_ctmc_rate_matrix_by_expected_substitution_rate(
            rate_matrix,
            stationary,
            state_labels=DNA_STATE_ORDER,
        )
    except PhylogeneticsError as error:
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires a positive finite expected substitution rate"
        ) from error


def is_dna_transition(left_state: str, right_state: str) -> bool:
    return (left_state, right_state) in DNA_TRANSITION_PAIRS


def validate_positive_kappa(
    kappa: float,
    *,
    model_name: str,
) -> float:
    if not math.isfinite(kappa) or kappa <= 0.0:
        raise ValueError(f"{model_name} kappa must be a finite positive value")
    return float(kappa)


def one_hot_dna_leaf_vector(
    states_by_taxon: dict[str, str],
    *,
    model_name: str,
    node_name: str | None,
) -> numpy.ndarray:
    if node_name is None:
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires named tree tips for alignment lookup"
        )
    vector = numpy.zeros(4, dtype=float)
    vector[DNA_STATE_INDEX[states_by_taxon[node_name]]] = 1.0
    return vector


def evaluate_fixed_topology_dna_site_log_likelihood(
    tree: PhyloTree,
    states: tuple[str, ...],
    *,
    taxon_order: list[str],
    model_name: str,
    observation_policy: str,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    transition_matrix_for_child,
) -> float:
    """Evaluate one DNA site log likelihood on one fixed topology."""
    from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
        dna_observation_state_order,
        resolve_dna_observation_leaf_vector,
    )

    validated_root_prior = numpy.asarray(root_prior, dtype=float)
    states_by_taxon = dict(zip(taxon_order, states, strict=True))
    state_count = len(dna_observation_state_order(observation_policy=observation_policy))
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=state_count,
        leaf_likelihood=lambda node: resolve_dna_observation_leaf_vector(
            states_by_taxon,
            model_name=model_name,
            node_name=node.name,
            observation_policy=observation_policy,
        ),
        transition_matrix_for_child=transition_matrix_for_child,
    )
    return log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=validated_root_prior,
    )


def _dna_probability_error(
    *,
    owner_name: str,
    parameter_name: str,
    error: PhylogeneticsError,
) -> InvalidAlignmentError:
    if error.code == "categorical_probability_vector_value_not_finite":
        return InvalidAlignmentError(
            f"{owner_name} {parameter_name} must all be finite"
        )
    if error.code == "categorical_probability_vector_value_negative":
        return InvalidAlignmentError(
            f"{owner_name} {parameter_name} must be non-negative"
        )
    if error.code == "categorical_probability_vector_not_normalized":
        return InvalidAlignmentError(
            f"{owner_name} {parameter_name} must sum to one within tolerance"
        )
    if error.code in {
        "categorical_probability_vector_missing_states",
        "categorical_probability_vector_unexpected_states",
    }:
        return InvalidAlignmentError(
            f"{owner_name} requires {parameter_name} for exactly A, C, G, and T"
        )
    return InvalidAlignmentError(str(error))

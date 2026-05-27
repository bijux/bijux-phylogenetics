from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import AlignmentTaxonMismatchError
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

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
    if not numpy.all(numpy.isfinite(vector)):
        raise InvalidAlignmentError(
            f"{model_name} likelihood base frequencies must all be finite"
        )
    if numpy.any(vector < 0.0):
        raise InvalidAlignmentError(
            f"{model_name} likelihood base frequencies must be nonnegative"
        )
    total = float(vector.sum())
    if total <= 0.0:
        raise InvalidAlignmentError(
            f"{model_name} likelihood base frequencies must sum to a positive value"
        )
    return vector / total


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
            vector = numpy.array(
                [
                    float(exchangeabilities[pair])
                    for pair in DNA_EXCHANGEABILITY_ORDER
                ],
                dtype=float,
            )
        elif set(exchangeabilities) == {
            "".join(pair) for pair in DNA_EXCHANGEABILITY_ORDER
        }:
            vector = numpy.array(
                [
                    float(exchangeabilities["".join(pair)])
                    for pair in DNA_EXCHANGEABILITY_ORDER
                ],
                dtype=float,
            )
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
    expected_rate = -float(numpy.sum(stationary * numpy.diag(rate_matrix)))
    if expected_rate <= 0.0 or not math.isfinite(expected_rate):
        raise InvalidAlignmentError(
            f"{model_name} likelihood requires a positive finite expected substitution rate"
        )
    return rate_matrix / expected_rate


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
    root_prior: numpy.ndarray | list[float] | tuple[float, ...],
    transition_matrix_for_child,
) -> float:
    """Evaluate one DNA site log likelihood on one fixed topology."""
    states_by_taxon = dict(zip(taxon_order, states, strict=True))
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(DNA_STATE_ORDER),
        leaf_likelihood=lambda node: one_hot_dna_leaf_vector(
            states_by_taxon,
            model_name=model_name,
            node_name=node.name,
        ),
        transition_matrix_for_child=transition_matrix_for_child,
    )
    return log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=root_prior,
    )

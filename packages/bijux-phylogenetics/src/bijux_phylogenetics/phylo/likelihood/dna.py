from __future__ import annotations

import math

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.patterns import CompressedAlignmentSitePatterns
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import AlignmentTaxonMismatchError
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

DNA_STATE_ORDER = ("A", "C", "G", "T")
DNA_STATE_INDEX = {state: index for index, state in enumerate(DNA_STATE_ORDER)}
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


def validate_explicit_branch_lengths(
    tree: PhyloTree,
    *,
    model_name: str,
) -> None:
    for _parent, child in tree.iter_edges():
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                f"{model_name} fixed-topology likelihood requires explicit branch lengths on every edge"
            )
        if child.branch_length < 0.0:
            raise InvalidBranchLengthError(
                f"{model_name} likelihood does not accept negative branch lengths"
            )


def validate_tree_taxa_against_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    model_name: str,
) -> None:
    tree_taxa = [leaf.name for leaf in tree.iter_leaves()]
    if any(name is None for name in tree_taxa):
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires every tree tip to have a matching alignment identifier"
        )
    observed_tree_taxa = [name for name in tree_taxa if name is not None]
    if len(set(observed_tree_taxa)) != len(observed_tree_taxa):
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires uniquely named tree tips"
        )
    expected_taxa = compressed_patterns.taxon_order
    if set(observed_tree_taxa) != set(expected_taxa):
        missing_from_alignment = sorted(set(observed_tree_taxa) - set(expected_taxa))
        missing_from_tree = sorted(set(expected_taxa) - set(observed_tree_taxa))
        details: list[str] = []
        if missing_from_alignment:
            details.append(f"tree-only taxa: {', '.join(missing_from_alignment)}")
        if missing_from_tree:
            details.append(f"alignment-only taxa: {', '.join(missing_from_tree)}")
        detail_suffix = f" ({'; '.join(details)})" if details else ""
        raise AlignmentTaxonMismatchError(
            f"{model_name} likelihood requires identical tree and alignment taxon sets"
            f"{detail_suffix}"
        )


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

from __future__ import annotations

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.patterns import CompressedAlignmentSitePatterns
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import AlignmentTaxonMismatchError
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

DNA_STATE_ORDER = ("A", "C", "G", "T")
DNA_STATE_INDEX = {state: index for index, state in enumerate(DNA_STATE_ORDER)}
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

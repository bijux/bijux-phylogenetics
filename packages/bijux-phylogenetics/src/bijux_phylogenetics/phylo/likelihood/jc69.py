from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import Jc69TreeLikelihoodReport
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    sum_compressed_site_pattern_log_likelihoods,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import AlignmentTaxonMismatchError
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

_JC69_STATE_ORDER = ("A", "C", "G", "T")
_JC69_STATE_INDEX = {state: index for index, state in enumerate(_JC69_STATE_ORDER)}
_JC69_ROOT_PRIOR = numpy.full(4, 0.25, dtype=float)


def jc69_rate_matrix() -> numpy.ndarray:
    """Return the normalized JC69 rate matrix with expected rate one."""
    rate_matrix = numpy.full((4, 4), 1.0 / 3.0, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def jc69_transition_probability_matrix(branch_length: float) -> numpy.ndarray:
    """Return the native closed-form JC69 transition matrix for one branch."""
    if branch_length <= 0.0:
        return numpy.eye(4, dtype=float)
    decay = math.exp((-4.0 * branch_length) / 3.0)
    same_probability = 0.25 + (0.75 * decay)
    different_probability = 0.25 - (0.25 * decay)
    transition = numpy.full((4, 4), different_probability, dtype=float)
    numpy.fill_diagonal(transition, same_probability)
    return transition


def evaluate_jc69_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
) -> Jc69TreeLikelihoodReport:
    """Evaluate one fixed-topology JC69 likelihood from aligned DNA records."""
    normalized_records = _normalized_jc69_records(records)
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    return _evaluate_jc69_tree_likelihood_from_patterns(tree, compressed_patterns)


def evaluate_jc69_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
) -> Jc69TreeLikelihoodReport:
    """Evaluate one fixed-topology JC69 likelihood from one tree path and alignment."""
    return evaluate_jc69_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
    )


def _evaluate_jc69_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
) -> Jc69TreeLikelihoodReport:
    _validate_explicit_branch_lengths(tree)
    _validate_tree_taxa_against_patterns(tree, compressed_patterns)

    transition_by_node_id = {
        child.node_id: jc69_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0)
        )
        for _parent, child in tree.iter_edges()
    }

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        states_by_taxon = dict(zip(compressed_patterns.taxon_order, states, strict=True))
        pruning_pass = postorder_conditional_likelihoods(
            tree,
            state_count=4,
            leaf_likelihood=lambda node: _one_hot_leaf_vector(
                states_by_taxon,
                node_name=node.name,
            ),
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )
        return log_likelihood_from_root_prior(
            tree,
            pruning_pass,
            root_prior=_JC69_ROOT_PRIOR,
        )

    log_likelihood = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )
    return Jc69TreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        log_likelihood=log_likelihood,
    )


def _normalized_jc69_records(records: list[AlignmentRecord]) -> list[AlignmentRecord]:
    normalized_records: list[AlignmentRecord] = []
    for record in records:
        normalized_sequence = record.sequence.upper()
        invalid_states = sorted(
            {
                state
                for state in normalized_sequence
                if state not in _JC69_STATE_INDEX
            }
        )
        if invalid_states:
            joined_states = ", ".join(invalid_states)
            raise InvalidAlignmentError(
                "JC69 likelihood currently requires unambiguous DNA states A, C, G, and T only; "
                f"record '{record.identifier}' contains {joined_states}"
            )
        normalized_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=normalized_sequence,
            )
        )
    return normalized_records


def _validate_explicit_branch_lengths(tree: PhyloTree) -> None:
    for _parent, child in tree.iter_edges():
        if child.branch_length is None:
            raise InvalidBranchLengthError(
                "JC69 fixed-topology likelihood requires explicit branch lengths on every edge"
            )
        if child.branch_length < 0.0:
            raise InvalidBranchLengthError(
                "JC69 likelihood does not accept negative branch lengths"
            )


def _validate_tree_taxa_against_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
) -> None:
    tree_taxa = [leaf.name for leaf in tree.iter_leaves()]
    if any(name is None for name in tree_taxa):
        raise AlignmentTaxonMismatchError(
            "JC69 likelihood requires every tree tip to have a matching alignment identifier"
        )
    observed_tree_taxa = [name for name in tree_taxa if name is not None]
    if len(set(observed_tree_taxa)) != len(observed_tree_taxa):
        raise AlignmentTaxonMismatchError(
            "JC69 likelihood requires uniquely named tree tips"
        )
    expected_taxa = compressed_patterns.taxon_order
    if set(observed_tree_taxa) != set(expected_taxa):
        missing_from_alignment = sorted(set(observed_tree_taxa) - set(expected_taxa))
        missing_from_tree = sorted(set(expected_taxa) - set(observed_tree_taxa))
        details: list[str] = []
        if missing_from_alignment:
            details.append(
                f"tree-only taxa: {', '.join(missing_from_alignment)}"
            )
        if missing_from_tree:
            details.append(f"alignment-only taxa: {', '.join(missing_from_tree)}")
        raise AlignmentTaxonMismatchError(
            "JC69 likelihood requires identical tree and alignment taxon sets"
            + (f" ({'; '.join(details)})" if details else "")
        )


def _one_hot_leaf_vector(
    states_by_taxon: dict[str, str],
    *,
    node_name: str | None,
) -> numpy.ndarray:
    if node_name is None:
        raise AlignmentTaxonMismatchError(
            "JC69 likelihood requires named tree tips for alignment lookup"
        )
    vector = numpy.zeros(4, dtype=float)
    vector[_JC69_STATE_INDEX[states_by_taxon[node_name]]] = 1.0
    return vector

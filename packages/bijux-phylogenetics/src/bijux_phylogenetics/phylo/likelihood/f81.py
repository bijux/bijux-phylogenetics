from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    estimate_empirical_dna_base_frequencies,
    normalize_unambiguous_dna_records,
    one_hot_dna_leaf_vector,
    validate_dna_base_frequencies,
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.likelihood.models import F81TreeLikelihoodReport
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
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


def f81_rate_matrix(base_frequencies: dict[str, float] | numpy.ndarray) -> numpy.ndarray:
    """Return the normalized F81 rate matrix with expected rate one."""
    stationary = validate_dna_base_frequencies(base_frequencies, model_name="F81")
    variability = _f81_variability(stationary)
    if variability <= 0.0:
        raise InvalidAlignmentError(
            "F81 likelihood requires at least two nucleotides with nonzero stationary frequency"
        )
    rate_matrix = numpy.zeros((4, 4), dtype=float)
    for row_index in range(len(DNA_STATE_ORDER)):
        row_total = 0.0
        for column_index in range(len(DNA_STATE_ORDER)):
            if row_index == column_index:
                continue
            rate = stationary[column_index] / variability
            rate_matrix[row_index, column_index] = rate
            row_total += rate
        rate_matrix[row_index, row_index] = -row_total
    return rate_matrix


def f81_transition_probability_matrix(
    branch_length: float,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray,
) -> numpy.ndarray:
    """Return the native closed-form F81 transition matrix for one branch."""
    stationary = validate_dna_base_frequencies(base_frequencies, model_name="F81")
    if branch_length <= 0.0:
        return numpy.eye(4, dtype=float)
    variability = _f81_variability(stationary)
    if variability <= 0.0:
        raise InvalidAlignmentError(
            "F81 likelihood requires at least two nucleotides with nonzero stationary frequency"
        )
    decay = math.exp(-branch_length / variability)
    transition_matrix = numpy.zeros((4, 4), dtype=float)
    for row_index in range(len(DNA_STATE_ORDER)):
        for column_index in range(len(DNA_STATE_ORDER)):
            if row_index == column_index:
                transition_matrix[row_index, column_index] = (
                    stationary[column_index]
                    + ((1.0 - stationary[column_index]) * decay)
                )
            else:
                transition_matrix[row_index, column_index] = (
                    stationary[column_index] * (1.0 - decay)
                )
    return transition_matrix


def evaluate_f81_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
) -> F81TreeLikelihoodReport:
    """Evaluate one fixed-topology F81 likelihood from aligned DNA records."""
    normalized_records = normalize_unambiguous_dna_records(records, model_name="F81")
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    if base_frequencies is None:
        stationary = estimate_empirical_dna_base_frequencies(normalized_records)
        source = "estimated"
    else:
        stationary = validate_dna_base_frequencies(base_frequencies, model_name="F81")
        source = "provided"
    return _evaluate_f81_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        stationary_frequencies=stationary,
        base_frequency_source=source,
    )


def evaluate_f81_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
) -> F81TreeLikelihoodReport:
    """Evaluate one fixed-topology F81 likelihood from one tree path and alignment."""
    return evaluate_f81_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        base_frequencies=base_frequencies,
    )


def _evaluate_f81_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    stationary_frequencies: numpy.ndarray,
    base_frequency_source: str,
) -> F81TreeLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name="F81")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="F81",
    )
    validated_frequencies = validate_dna_base_frequencies(
        stationary_frequencies,
        model_name="F81",
    )
    transition_by_node_id = {
        child.node_id: f81_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0),
            base_frequencies=validated_frequencies,
        )
        for _parent, child in tree.iter_edges()
    }

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        states_by_taxon = dict(zip(compressed_patterns.taxon_order, states, strict=True))
        pruning_pass = postorder_conditional_likelihoods(
            tree,
            state_count=4,
            leaf_likelihood=lambda node: one_hot_dna_leaf_vector(
                states_by_taxon,
                model_name="F81",
                node_name=node.name,
            ),
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )
        return log_likelihood_from_root_prior(
            tree,
            pruning_pass,
            root_prior=validated_frequencies,
        )

    log_likelihood = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )
    return F81TreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        base_frequency_source=base_frequency_source,
        base_frequency_a=float(validated_frequencies[0]),
        base_frequency_c=float(validated_frequencies[1]),
        base_frequency_g=float(validated_frequencies[2]),
        base_frequency_t=float(validated_frequencies[3]),
        parameter_count=3,
        log_likelihood=log_likelihood,
        aic=(-2.0 * log_likelihood) + (2.0 * 3.0),
    )


def _f81_variability(stationary_frequencies: numpy.ndarray) -> float:
    return 1.0 - float(numpy.sum(stationary_frequencies * stationary_frequencies))

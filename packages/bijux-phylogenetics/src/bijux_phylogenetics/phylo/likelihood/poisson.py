from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    ProteinPoissonTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    UNIFORM_PROTEIN_ROOT_PRIOR,
    evaluate_fixed_topology_protein_site_log_likelihood,
    normalize_protein_likelihood_records,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
    validate_site_log_likelihood_reconstruction,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def protein_poisson_rate_matrix() -> numpy.ndarray:
    """Return the normalized 20-state equal-input protein Poisson rate matrix."""
    state_count = len(PROTEIN_STATE_ORDER)
    off_diagonal_rate = 1.0 / (state_count - 1)
    rate_matrix = numpy.full((state_count, state_count), off_diagonal_rate, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for row_index in range(state_count):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def protein_poisson_transition_probability_matrix(
    branch_length: float,
) -> numpy.ndarray:
    """Return the native closed-form 20-state protein Poisson transition matrix."""
    state_count = len(PROTEIN_STATE_ORDER)
    if branch_length <= 0.0:
        return numpy.eye(state_count, dtype=float)
    decay = math.exp((-state_count * branch_length) / (state_count - 1))
    same_probability = (1.0 / state_count) + (
        ((state_count - 1.0) / state_count) * decay
    )
    different_probability = (1.0 / state_count) - ((1.0 / state_count) * decay)
    transition = numpy.full(
        (state_count, state_count), different_probability, dtype=float
    )
    numpy.fill_diagonal(transition, same_probability)
    return transition


def evaluate_protein_poisson_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> ProteinPoissonTreeLikelihoodReport:
    """Evaluate one fixed-topology protein Poisson likelihood from aligned amino acids."""
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="protein Poisson",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    return _evaluate_protein_poisson_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def evaluate_protein_poisson_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> ProteinPoissonTreeLikelihoodReport:
    """Evaluate one fixed-topology protein Poisson likelihood from paths."""
    return evaluate_protein_poisson_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def _evaluate_protein_poisson_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    gap_policy: str,
    missing_policy: str,
    ambiguity_policy: str,
) -> ProteinPoissonTreeLikelihoodReport:
    transition_by_node_id = {
        child.node_id: protein_poisson_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0)
        )
        for _parent, child in tree.iter_edges()
    }
    site_log_likelihoods, log_likelihood = (
        expanded_site_log_likelihood_rows_from_patterns(
            compressed_patterns,
            site_log_likelihood=lambda states: (
                evaluate_fixed_topology_protein_site_log_likelihood(
                    tree,
                    states,
                    taxon_order=compressed_patterns.taxon_order,
                    model_name="protein Poisson",
                    root_prior=UNIFORM_PROTEIN_ROOT_PRIOR,
                    transition_matrix_for_child=lambda child: transition_by_node_id[
                        child.node_id or ""
                    ],
                    gap_policy=gap_policy,
                    missing_policy=missing_policy,
                    ambiguity_policy=ambiguity_policy,
                )
            ),
        )
    )
    validate_site_log_likelihood_reconstruction(
        site_log_likelihoods,
        expected_total_log_likelihood=log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name="protein Poisson likelihood",
    )
    return ProteinPoissonTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        state_count=len(PROTEIN_STATE_ORDER),
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
        log_likelihood=log_likelihood,
        site_log_likelihoods=site_log_likelihoods,
    )

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
    UNIFORM_DNA_ROOT_PRIOR,
    is_dna_transition,
    normalize_unambiguous_dna_records,
    one_hot_dna_leaf_vector,
    validate_positive_kappa,
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    K80KappaOptimizationReport,
    K80TreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_likelihood_search,
)
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

def k80_rate_matrix(kappa: float) -> numpy.ndarray:
    """Return the normalized K80 rate matrix with expected rate one."""
    validate_positive_kappa(kappa, model_name="K80")
    transversion_rate = 1.0 / (kappa + 2.0)
    transition_rate = kappa * transversion_rate
    rate_matrix = numpy.zeros((4, 4), dtype=float)
    for left_index, left_state in enumerate(DNA_STATE_ORDER):
        row_total = 0.0
        for right_index, right_state in enumerate(DNA_STATE_ORDER):
            if left_index == right_index:
                continue
            rate = (
                transition_rate
                if is_dna_transition(left_state, right_state)
                else transversion_rate
            )
            rate_matrix[left_index, right_index] = rate
            row_total += rate
        rate_matrix[left_index, left_index] = -row_total
    return rate_matrix


def k80_transition_probability_matrix(
    branch_length: float,
    *,
    kappa: float,
) -> numpy.ndarray:
    """Return the native closed-form K80 transition matrix for one branch."""
    validate_positive_kappa(kappa, model_name="K80")
    if branch_length <= 0.0:
        return numpy.eye(4, dtype=float)
    transversion_rate = 1.0 / (kappa + 2.0)
    transition_rate = kappa * transversion_rate
    transversion_decay = math.exp(-4.0 * transversion_rate * branch_length)
    transition_decay = math.exp(
        -2.0 * (transition_rate + transversion_rate) * branch_length
    )
    same_probability = (
        0.25
        + (0.25 * transversion_decay)
        + (0.5 * transition_decay)
    )
    transition_probability = (
        0.25
        + (0.25 * transversion_decay)
        - (0.5 * transition_decay)
    )
    transversion_probability = 0.25 - (0.25 * transversion_decay)
    transition_matrix = numpy.zeros((4, 4), dtype=float)
    for left_index, left_state in enumerate(DNA_STATE_ORDER):
        for right_index, right_state in enumerate(DNA_STATE_ORDER):
            if left_index == right_index:
                transition_matrix[left_index, right_index] = same_probability
                continue
            transition_matrix[left_index, right_index] = (
                transition_probability
                if is_dna_transition(left_state, right_state)
                else transversion_probability
            )
    return transition_matrix


def evaluate_k80_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    kappa: float,
) -> K80TreeLikelihoodReport:
    """Evaluate one fixed-topology K80 likelihood from aligned DNA records."""
    normalized_records = normalize_unambiguous_dna_records(records, model_name="K80")
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    return _evaluate_k80_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        kappa=kappa,
    )


def evaluate_k80_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    kappa: float,
) -> K80TreeLikelihoodReport:
    """Evaluate one fixed-topology K80 likelihood from one tree path and alignment."""
    return evaluate_k80_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        kappa=kappa,
    )


def optimize_k80_kappa(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    initial_kappa: float = 1.0,
    lower_kappa_bound: float = 0.05,
    upper_kappa_bound: float = 20.0,
) -> K80KappaOptimizationReport:
    """Optimize one fixed-topology K80 kappa on fixed branch lengths."""
    validate_positive_kappa(initial_kappa, model_name="K80")
    validate_positive_kappa(lower_kappa_bound, model_name="K80")
    validate_positive_kappa(upper_kappa_bound, model_name="K80")
    if upper_kappa_bound <= lower_kappa_bound:
        raise ValueError("K80 kappa bounds must be strictly increasing")

    normalized_records = normalize_unambiguous_dna_records(records, model_name="K80")
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    working_tree = tree.copy()
    validate_explicit_branch_lengths(working_tree, model_name="K80")
    initial_report = _evaluate_k80_tree_likelihood_from_patterns(
        working_tree,
        compressed_patterns,
        kappa=initial_kappa,
    )

    def evaluate_candidate_kappa(
        candidate_kappa: float,
    ) -> tuple[K80TreeLikelihoodReport, float]:
        report = _evaluate_k80_tree_likelihood_from_patterns(
            working_tree,
            compressed_patterns,
            kappa=candidate_kappa,
        )
        return report, report.log_likelihood

    search_result = run_bounded_likelihood_search(
        lower_bound=lower_kappa_bound,
        upper_bound=upper_kappa_bound,
        evaluate=evaluate_candidate_kappa,
    )
    optimized_report = search_result.payload
    return K80KappaOptimizationReport(
        taxa=optimized_report.taxa,
        site_count=optimized_report.site_count,
        pattern_count=optimized_report.pattern_count,
        tree_newick=optimized_report.tree_newick,
        initial_kappa=initial_kappa,
        optimized_kappa=search_result.parameter_value,
        initial_log_likelihood=initial_report.log_likelihood,
        optimized_log_likelihood=optimized_report.log_likelihood,
        function_evaluation_count=search_result.function_evaluation_count + 1,
        converged=search_result.converged,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
    )


def optimize_k80_kappa_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    initial_kappa: float = 1.0,
    lower_kappa_bound: float = 0.05,
    upper_kappa_bound: float = 20.0,
) -> K80KappaOptimizationReport:
    """Optimize K80 kappa from one tree path and one alignment path."""
    return optimize_k80_kappa(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        initial_kappa=initial_kappa,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
    )


def _evaluate_k80_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    kappa: float,
) -> K80TreeLikelihoodReport:
    validate_positive_kappa(kappa, model_name="K80")
    validate_explicit_branch_lengths(tree, model_name="K80")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="K80",
    )

    transition_by_node_id = {
        child.node_id: k80_transition_probability_matrix(
            max(float(child.branch_length or 0.0), 0.0),
            kappa=kappa,
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
                model_name="K80",
                node_name=node.name,
            ),
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )
        return log_likelihood_from_root_prior(
            tree,
            pruning_pass,
            root_prior=UNIFORM_DNA_ROOT_PRIOR,
        )

    log_likelihood = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )
    return K80TreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        kappa=kappa,
        log_likelihood=log_likelihood,
    )

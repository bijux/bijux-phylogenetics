from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.models import (
    Jc69BranchLengthOptimizationReport,
    Jc69BranchLengthOptimizationStep,
)
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


def optimize_jc69_branch_lengths(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    lower_branch_length_bound: float = 1e-6,
    upper_branch_length_bound: float = 5.0,
    max_coordinate_passes: int = 12,
    improvement_tolerance: float = 1e-9,
) -> Jc69BranchLengthOptimizationReport:
    """Optimize one fixed topology under native JC69 branch likelihood."""
    if lower_branch_length_bound <= 0.0:
        raise InvalidBranchLengthError("JC69 branch-length lower bound must be positive")
    if upper_branch_length_bound <= lower_branch_length_bound:
        raise InvalidBranchLengthError(
            "JC69 branch-length bounds must be strictly increasing"
        )
    if max_coordinate_passes < 1:
        raise ValueError("max_coordinate_passes must be at least one")

    normalized_records = _normalized_jc69_records(records)
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    working_tree = tree.copy()
    _validate_explicit_branch_lengths(working_tree)
    edge_nodes = [child for _parent, child in working_tree.iter_edges()]
    initial_tree_newick = dumps_newick(tree)
    initial_report = _evaluate_jc69_tree_likelihood_from_patterns(
        working_tree,
        compressed_patterns,
    )
    current_report = initial_report
    function_evaluation_count = 1
    optimization_pass_count = 0
    converged = False
    steps: list[Jc69BranchLengthOptimizationStep] = []

    for optimization_pass in range(1, max_coordinate_passes + 1):
        optimization_pass_count = optimization_pass
        improved = False
        for node in edge_nodes:
            if node.node_id is None:
                raise ValueError("tree node is missing a stable node_id")
            starting_branch_length = float(node.branch_length or 0.0)
            starting_log_likelihood = current_report.log_likelihood

            def evaluate_candidate(
                branch_length: float,
            ) -> tuple[Jc69TreeLikelihoodReport, float]:
                node.branch_length = branch_length
                report = _evaluate_jc69_tree_likelihood_from_patterns(
                    working_tree,
                    compressed_patterns,
                )
                return report, report.log_likelihood

            (
                optimized_branch_length,
                optimized_report,
                optimized_log_likelihood,
                search_function_evaluation_count,
            ) = _run_bounded_jc69_branch_search(
                lower_bound=lower_branch_length_bound,
                upper_bound=upper_branch_length_bound,
                evaluate=evaluate_candidate,
            )
            function_evaluation_count += search_function_evaluation_count
            accepted = optimized_log_likelihood > (
                current_report.log_likelihood + improvement_tolerance
            )
            if accepted:
                node.branch_length = optimized_branch_length
                current_report = optimized_report
                improved = True
            else:
                node.branch_length = starting_branch_length
            steps.append(
                Jc69BranchLengthOptimizationStep(
                    optimization_pass=optimization_pass,
                    branch_id=node.node_id,
                    child_name=node.name,
                    descendant_taxa=node.descendant_taxa,
                    starting_branch_length=starting_branch_length,
                    optimized_branch_length=optimized_branch_length,
                    starting_log_likelihood=starting_log_likelihood,
                    optimized_log_likelihood=optimized_log_likelihood,
                    accepted=accepted,
                )
            )
        if not improved:
            converged = True
            break

    return Jc69BranchLengthOptimizationReport(
        taxa=current_report.taxa,
        site_count=current_report.site_count,
        pattern_count=current_report.pattern_count,
        branch_count=len(edge_nodes),
        initial_tree_newick=initial_tree_newick,
        optimized_tree_newick=dumps_newick(working_tree),
        initial_log_likelihood=initial_report.log_likelihood,
        optimized_log_likelihood=current_report.log_likelihood,
        optimization_pass_count=optimization_pass_count,
        function_evaluation_count=function_evaluation_count,
        converged=converged,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        steps=steps,
    )


def optimize_jc69_branch_lengths_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    lower_branch_length_bound: float = 1e-6,
    upper_branch_length_bound: float = 5.0,
    max_coordinate_passes: int = 12,
    improvement_tolerance: float = 1e-9,
) -> Jc69BranchLengthOptimizationReport:
    """Optimize one fixed topology from one tree path and one alignment path."""
    return optimize_jc69_branch_lengths(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
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


def _run_bounded_jc69_branch_search(
    *,
    lower_bound: float,
    upper_bound: float,
    evaluate,
    tolerance: float = 1e-9,
    max_iterations: int = 400,
) -> tuple[float, Jc69TreeLikelihoodReport, float, int]:
    if upper_bound <= lower_bound:
        raise InvalidBranchLengthError("JC69 branch-length bounds must be strictly increasing")
    if tolerance <= 0.0:
        raise ValueError("JC69 branch-length search tolerance must be positive")
    if max_iterations < 1:
        raise ValueError("JC69 branch-length search iterations must be positive")

    phi = (math.sqrt(5.0) - 1.0) / 2.0
    cache: dict[float, tuple[Jc69TreeLikelihoodReport, float]] = {}

    def evaluate_cached(branch_length: float) -> tuple[Jc69TreeLikelihoodReport, float]:
        cache_key = round(branch_length, 12)
        cached = cache.get(cache_key)
        if cached is not None:
            return cached
        candidate = evaluate(branch_length)
        cache[cache_key] = candidate
        return candidate

    left = upper_bound - (phi * (upper_bound - lower_bound))
    right = lower_bound + (phi * (upper_bound - lower_bound))
    left_report, left_objective = evaluate_cached(left)
    right_report, right_objective = evaluate_cached(right)

    iteration = 0
    while (upper_bound - lower_bound) > tolerance and iteration < max_iterations:
        if left_objective < right_objective:
            lower_bound = left
            left = right
            left_report = right_report
            left_objective = right_objective
            right = lower_bound + (phi * (upper_bound - lower_bound))
            right_report, right_objective = evaluate_cached(right)
        else:
            upper_bound = right
            right = left
            right_report = left_report
            right_objective = left_objective
            left = upper_bound - (phi * (upper_bound - lower_bound))
            left_report, left_objective = evaluate_cached(left)
        iteration += 1

    midpoint = (lower_bound + upper_bound) / 2.0
    midpoint_report, midpoint_objective = evaluate_cached(midpoint)
    ranked = sorted(
        [
            (left, left_report, left_objective),
            (right, right_report, right_objective),
            (midpoint, midpoint_report, midpoint_objective),
        ],
        key=lambda item: (-item[2], item[0]),
    )
    best_branch_length, best_report, best_objective = ranked[0]
    return best_branch_length, best_report, best_objective, len(cache)

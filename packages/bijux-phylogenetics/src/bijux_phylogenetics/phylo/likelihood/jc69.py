from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    UNIFORM_DNA_ROOT_PRIOR,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    augment_dna_rate_matrix_with_gap_state,
    dna_observation_state_order,
    estimate_empirical_gap_state_frequency,
    normalize_dna_likelihood_records,
    resolve_default_dna_root_prior_for_observation_policy,
    resolve_dna_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    Jc69BranchLengthOptimizationReport,
    Jc69BranchLengthOptimizationStep,
    Jc69TreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.parameter_bounds import (
    validate_increasing_parameter_bounds,
    validate_parameter_within_bounds,
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
    transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    sum_compressed_site_pattern_log_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError


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
    *,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> Jc69TreeLikelihoodReport:
    """Evaluate one fixed-topology JC69 likelihood from aligned DNA records."""
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="JC69",
        observation_policy=observation_policy,
    )
    gap_state_frequency = (
        estimate_empirical_gap_state_frequency(
            normalized_records,
            model_name="JC69",
        )
        if observation_policy == "fifth-state"
        else None
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    resolved_root_prior = resolve_default_dna_root_prior_for_observation_policy(
        normalized_records,
        owner_name="JC69 likelihood",
        default_policy="equal",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=UNIFORM_DNA_ROOT_PRIOR,
        observation_policy=observation_policy,
    )
    return _evaluate_jc69_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        root_prior=resolved_root_prior.root_prior,
        observation_policy=observation_policy,
        gap_state_frequency=gap_state_frequency,
    )


def evaluate_jc69_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> Jc69TreeLikelihoodReport:
    """Evaluate one fixed-topology JC69 likelihood from one tree path and alignment."""
    return evaluate_jc69_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        observation_policy=observation_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
        raise InvalidBranchLengthError(
            "JC69 branch-length lower bound must be positive"
        )
    try:
        validated_lower_bound, validated_upper_bound = (
            validate_increasing_parameter_bounds(
                parameter_name="branch length",
                lower_bound=lower_branch_length_bound,
                upper_bound=upper_branch_length_bound,
                owner_name="JC69 branch-length optimization",
            )
        )
    except ValueError as error:
        raise InvalidBranchLengthError(str(error)) from error
    if max_coordinate_passes < 1:
        raise ValueError("max_coordinate_passes must be at least one")

    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="JC69",
        observation_policy="reject",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    working_tree = tree.copy()
    validate_explicit_branch_lengths(working_tree, model_name="JC69")
    for _parent, child in working_tree.iter_edges():
        if child.branch_length is None:
            continue
        try:
            validate_parameter_within_bounds(
                parameter_name="branch length",
                value=float(child.branch_length),
                lower_bound=validated_lower_bound,
                upper_bound=validated_upper_bound,
                owner_name="JC69 branch-length optimization",
            )
        except ValueError as error:
            raise InvalidBranchLengthError(str(error)) from error
    edge_nodes = [child for _parent, child in working_tree.iter_edges()]
    initial_tree_newick = dumps_newick(tree)
    initial_report = _evaluate_jc69_tree_likelihood_from_patterns(
        working_tree,
        compressed_patterns,
        root_prior=UNIFORM_DNA_ROOT_PRIOR,
        observation_policy="reject",
        gap_state_frequency=None,
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
            current_node = node

            def evaluate_candidate(
                branch_length: float,
                current_node=current_node,
            ) -> tuple[Jc69TreeLikelihoodReport, float]:
                current_node.branch_length = branch_length
                report = _evaluate_jc69_tree_likelihood_from_patterns(
                    working_tree,
                    compressed_patterns,
                    root_prior=UNIFORM_DNA_ROOT_PRIOR,
                    observation_policy="reject",
                    gap_state_frequency=None,
                )
                return report, report.log_likelihood

            search_result = run_bounded_likelihood_search(
                lower_bound=validated_lower_bound,
                upper_bound=validated_upper_bound,
                evaluate=evaluate_candidate,
            )
            function_evaluation_count += search_result.function_evaluation_count
            optimized_branch_length = float(search_result.parameter_value)
            optimized_log_likelihood = search_result.objective_value
            accepted = optimized_log_likelihood > (
                current_report.log_likelihood + improvement_tolerance
            )
            if accepted:
                node.branch_length = optimized_branch_length
                current_report = search_result.payload
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
        lower_branch_length_bound=validated_lower_bound,
        upper_branch_length_bound=validated_upper_bound,
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
    *,
    root_prior: numpy.ndarray,
    observation_policy: str,
    gap_state_frequency: float | None,
) -> Jc69TreeLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name="JC69")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="JC69",
    )
    state_order = dna_observation_state_order(observation_policy=observation_policy)
    if observation_policy == "fifth-state":
        if gap_state_frequency is None:
            raise ValueError(
                "JC69 fifth-state observation policy requires an explicit gap_state_frequency"
            )
        augmented_rate_matrix = augment_dna_rate_matrix_with_gap_state(
            jc69_rate_matrix(),
            nucleotide_frequencies=UNIFORM_DNA_ROOT_PRIOR,
            gap_state_frequency=gap_state_frequency,
            model_name="JC69",
        )
        transition_by_node_id = {
            child.node_id: transition_probability_matrix(
                augmented_rate_matrix,
                max(float(child.branch_length or 0.0), 0.0),
            )
            for _parent, child in tree.iter_edges()
        }
    else:
        transition_by_node_id = {
            child.node_id: jc69_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0)
            )
            for _parent, child in tree.iter_edges()
        }

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        states_by_taxon = dict(
            zip(compressed_patterns.taxon_order, states, strict=True)
        )
        pruning_pass = postorder_conditional_likelihoods(
            tree,
            state_count=len(state_order),
            leaf_likelihood=lambda node: resolve_dna_observation_leaf_vector(
                states_by_taxon,
                model_name="JC69",
                node_name=node.name,
                observation_policy=observation_policy,
            ),
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )
        return log_likelihood_from_root_prior(
            tree,
            pruning_pass,
            root_prior=root_prior,
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
        state_count=len(state_order),
        observation_policy=observation_policy,
        log_likelihood=log_likelihood,
    )

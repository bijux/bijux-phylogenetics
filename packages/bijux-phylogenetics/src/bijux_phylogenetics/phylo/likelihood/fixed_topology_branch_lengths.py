from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    evaluate_fixed_topology_dna_site_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    normalize_dna_likelihood_records,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    BranchLengthOptimizationRow,
    FixedTopologyNucleotideBranchLengthOptimizationReport,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    SelectedNucleotideLikelihoodSpecification,
    resolve_selected_nucleotide_likelihood_specification,
    validate_selected_nucleotide_likelihood_model,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError


@dataclass(frozen=True, slots=True)
class BranchReoptimizationResult:
    """Optimized branch-length result for one fixed topology under one likelihood surface."""

    optimized_tree: PhyloTree
    log_likelihood: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool


def validate_fixed_topology_nucleotide_branch_length_model(model_name: str) -> str:
    """Validate one selected nucleotide model for fixed-topology branch optimization."""
    return validate_selected_nucleotide_likelihood_model(model_name)


def optimize_fixed_topology_nucleotide_branch_lengths(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    observation_policy: str = "reject",
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> FixedTopologyNucleotideBranchLengthOptimizationReport:
    """Optimize all branch lengths on one fixed nucleotide topology."""
    normalized_model_name = validate_fixed_topology_nucleotide_branch_length_model(
        model_name
    )
    if lower_branch_length_bound < 0.0:
        raise InvalidBranchLengthError(
            "fixed-topology nucleotide branch-length lower bound must be nonnegative"
        )
    if upper_branch_length_bound <= lower_branch_length_bound:
        raise InvalidBranchLengthError(
            "fixed-topology nucleotide branch-length bounds must be strictly increasing"
        )
    if max_coordinate_passes < 1:
        raise ValueError("max_coordinate_passes must be at least one")
    normalized_observation_policy = observation_policy.strip().lower()
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name=f"{normalized_model_name.upper()} fixed-topology branch optimization",
        observation_policy=normalized_observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} fixed-topology branch optimization",
        observation_policy=normalized_observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    starting_tree = tree.copy().refresh()
    for _parent, child in starting_tree.iter_edges():
        branch_length = float(child.branch_length or 0.0)
        if not (
            lower_branch_length_bound <= branch_length <= upper_branch_length_bound
        ):
            raise InvalidBranchLengthError(
                "fixed-topology nucleotide branch optimization requires every starting branch length to lie within the declared bounds"
            )
    initial_tree_newick = dumps_newick(starting_tree)
    initial_log_likelihood = evaluate_selected_nucleotide_log_likelihood_from_patterns(
        starting_tree,
        compressed_patterns,
        specification=specification,
    )
    initial_branch_lengths = _collect_branch_lengths_by_id(starting_tree)
    search_result = optimize_selected_nucleotide_branch_lengths(
        starting_tree,
        compressed_patterns,
        specification=specification,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    optimized_branch_lengths = _collect_branch_lengths_by_id(search_result.optimized_tree)
    return FixedTopologyNucleotideBranchLengthOptimizationReport(
        model_name=specification.model_name,
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        branch_count=len(initial_branch_lengths),
        initial_tree_newick=initial_tree_newick,
        optimized_tree_newick=dumps_newick(search_result.optimized_tree),
        state_count=specification.state_count,
        observation_policy=specification.observation_policy,
        root_prior_source=specification.root_prior_source,
        parameter_count=len(specification.parameter_values),
        fixed_parameter_values=dict(specification.parameter_values),
        initial_log_likelihood=initial_log_likelihood,
        optimized_log_likelihood=search_result.log_likelihood,
        optimization_pass_count=search_result.optimization_pass_count,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        branches=[
            BranchLengthOptimizationRow(
                branch_id=node.node_id or "",
                child_name=node.name,
                descendant_taxa=node.descendant_taxa,
                initial_branch_length=initial_branch_lengths[node.node_id or ""],
                optimized_branch_length=optimized_branch_lengths[node.node_id or ""],
            )
            for _parent, node in starting_tree.iter_edges()
        ],
    )


def optimize_fixed_topology_nucleotide_branch_lengths_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    observation_policy: str = "reject",
    kappa: float | None = None,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float] | numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> FixedTopologyNucleotideBranchLengthOptimizationReport:
    """Optimize all branch lengths from one tree path and one alignment path."""
    return optimize_fixed_topology_nucleotide_branch_lengths(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def optimize_selected_nucleotide_branch_lengths(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    specification: SelectedNucleotideLikelihoodSpecification,
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> BranchReoptimizationResult:
    """Reoptimize one fixed-topology nucleotide tree under one resolved likelihood surface."""
    if lower_branch_length_bound < 0.0:
        raise InvalidBranchLengthError(
            "nucleotide likelihood branch-length lower bound must be nonnegative"
        )
    if upper_branch_length_bound <= lower_branch_length_bound:
        raise InvalidBranchLengthError(
            "nucleotide likelihood branch-length bounds must be strictly increasing"
        )
    if max_coordinate_passes < 1:
        raise ValueError("max_coordinate_passes must be at least one")

    working_tree = tree.copy().refresh()
    validate_explicit_branch_lengths(working_tree, model_name=specification.model_name)
    validate_tree_taxa_against_patterns(
        working_tree,
        compressed_patterns,
        model_name=specification.model_name,
    )
    edge_nodes = [child for _parent, child in working_tree.iter_edges()]
    initial_values: dict[str, float] = {}
    bounds_by_name: dict[str, tuple[float, float]] = {}
    for node in edge_nodes:
        if node.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        branch_length = float(node.branch_length or 0.0)
        if not (
            lower_branch_length_bound <= branch_length <= upper_branch_length_bound
        ):
            raise InvalidBranchLengthError(
                "every starting branch length must lie within the declared optimization bounds"
            )
        initial_values[node.node_id] = branch_length
        bounds_by_name[node.node_id] = (
            lower_branch_length_bound,
            upper_branch_length_bound,
        )

    def evaluate_candidate(
        branch_lengths_by_id: dict[str, float],
    ) -> tuple[float, float]:
        assign_branch_lengths(working_tree, branch_lengths_by_id)
        log_likelihood = evaluate_selected_nucleotide_log_likelihood_from_patterns(
            working_tree,
            compressed_patterns,
            specification=specification,
        )
        return log_likelihood, log_likelihood

    search_result = run_bounded_coordinate_likelihood_search(
        initial_values=initial_values,
        bounds_by_name=bounds_by_name,
        evaluate=evaluate_candidate,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    assign_branch_lengths(working_tree, search_result.parameter_values)
    return BranchReoptimizationResult(
        optimized_tree=working_tree.refresh(),
        log_likelihood=float(search_result.objective_value),
        optimization_pass_count=search_result.optimization_pass_count,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
    )


def assign_branch_lengths(tree: PhyloTree, branch_lengths_by_id: dict[str, float]) -> None:
    """Assign one branch-length vector to one refreshed tree by node identity."""
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        child.branch_length = branch_lengths_by_id[child.node_id]


def evaluate_selected_nucleotide_log_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    specification: SelectedNucleotideLikelihoodSpecification,
) -> float:
    """Evaluate one resolved selected nucleotide likelihood on compressed site patterns."""
    validate_explicit_branch_lengths(tree, model_name=specification.model_name)
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name=specification.model_name,
    )
    total_log_likelihood = 0.0
    for pattern in compressed_patterns.patterns:
        total_log_likelihood += pattern.weight * evaluate_fixed_topology_dna_site_log_likelihood(
            tree,
            pattern.states,
            taxon_order=compressed_patterns.taxon_order,
            model_name=specification.model_name,
            observation_policy=specification.observation_policy,
            root_prior=specification.root_prior,
            transition_matrix_for_child=lambda child: (
                specification.transition_matrix_for_branch_length(
                    max(float(child.branch_length or 0.0), 0.0)
                )
            ),
        )
    return total_log_likelihood


def _collect_branch_lengths_by_id(tree: PhyloTree) -> dict[str, float]:
    branch_lengths: dict[str, float] = {}
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        branch_lengths[child.node_id] = float(child.branch_length or 0.0)
    return branch_lengths

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    evaluate_fixed_topology_dna_site_log_likelihood,
    normalize_unambiguous_dna_records,
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
from bijux_phylogenetics.phylo.likelihood.substitution_parameters import (
    optimize_nucleotide_substitution_parameters,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES = frozenset({"coordinate-branch-lengths"})


@dataclass(frozen=True, slots=True)
class ResolvedNucleotideTopologySearchSurface:
    """Resolved nucleotide likelihood surface used during topology search."""

    model_name: str
    substitution_parameter_policy: str
    substitution_parameter_values: dict[str, float]
    substitution_parameter_warnings: list[str]
    specification: SelectedNucleotideLikelihoodSpecification


@dataclass(frozen=True, slots=True)
class BranchReoptimizationResult:
    """Optimized branch-length result for one fixed topology under one likelihood surface."""

    optimized_tree: PhyloTree
    log_likelihood: float
    optimization_pass_count: int
    function_evaluation_count: int
    converged: bool


def resolve_nucleotide_topology_search_tree(
    tree: PhyloTree | Path,
) -> tuple[PhyloTree, Path | None]:
    """Resolve one topology-search tree input and preserve its path when present."""
    resolved_tree_path = tree if isinstance(tree, Path) else None
    resolved_tree = load_tree(tree) if isinstance(tree, Path) else tree.copy()
    return resolved_tree.refresh(), resolved_tree_path


def resolve_nucleotide_topology_search_records(
    records: list[AlignmentRecord] | Path,
) -> tuple[list[AlignmentRecord], Path | None]:
    """Resolve one topology-search alignment input and preserve its path when present."""
    resolved_alignment_path = records if isinstance(records, Path) else None
    resolved_records = (
        load_fasta_alignment(records) if isinstance(records, Path) else list(records)
    )
    return resolved_records, resolved_alignment_path


def validate_nucleotide_topology_search_tree(
    tree: PhyloTree,
    *,
    workflow_name: str,
) -> None:
    """Require one structurally valid strictly binary rooted tree."""
    errors = tree.validation_errors()
    if errors:
        raise ValueError(
            f"{workflow_name} requires a structurally valid tree: {'; '.join(errors)}"
        )
    if len(tree.root.children) != 2:
        raise ValueError(f"{workflow_name} requires a rooted binary tree")
    invalid_internal_nodes = [
        node.node_id
        for node in tree.iter_internal_nodes(order="preorder")
        if len(node.children) != 2
    ]
    if invalid_internal_nodes:
        raise ValueError(
            f"{workflow_name} requires a strictly binary tree: {invalid_internal_nodes}"
        )


def validate_branch_reoptimization_policy(policy: str) -> str:
    normalized_policy = policy.strip().lower()
    if normalized_policy not in _SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES:
        raise ValueError(
            "branch_reoptimization_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES))
        )
    return normalized_policy


def resolve_nucleotide_topology_search_surface(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
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
) -> ResolvedNucleotideTopologySearchSurface:
    """Resolve one selected nucleotide likelihood surface for topology search."""
    normalized_model_name = validate_selected_nucleotide_likelihood_model(model_name)
    normalized_records = normalize_unambiguous_dna_records(
        records,
        model_name=f"{normalized_model_name.upper()} likelihood NNI search",
    )

    if normalized_model_name == "jc69":
        specification = resolve_selected_nucleotide_likelihood_specification(
            normalized_records,
            model_name=normalized_model_name,
            owner_name="JC69 likelihood NNI search",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideTopologySearchSurface(
            model_name=specification.model_name,
            substitution_parameter_policy="fixed-from-model",
            substitution_parameter_values=specification.parameter_values,
            substitution_parameter_warnings=[],
            specification=specification,
        )

    if normalized_model_name == "f81" and base_frequencies is None:
        specification = resolve_selected_nucleotide_likelihood_specification(
            normalized_records,
            model_name=normalized_model_name,
            owner_name="F81 likelihood NNI search",
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideTopologySearchSurface(
            model_name=specification.model_name,
            substitution_parameter_policy="estimated-from-alignment",
            substitution_parameter_values=specification.parameter_values,
            substitution_parameter_warnings=[],
            specification=specification,
        )

    parameters_are_fully_provided = (
        (normalized_model_name == "k80" and kappa is not None)
        or (
            normalized_model_name == "f81"
            and base_frequencies is not None
        )
        or (
            normalized_model_name == "hky85"
            and kappa is not None
        )
        or (
            normalized_model_name == "gtr"
            and exchangeabilities is not None
        )
    )
    if parameters_are_fully_provided:
        specification = resolve_selected_nucleotide_likelihood_specification(
            normalized_records,
            model_name=normalized_model_name,
            owner_name=f"{normalized_model_name.upper()} likelihood NNI search",
            kappa=kappa,
            base_frequencies=base_frequencies,
            exchangeabilities=exchangeabilities,
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
        )
        return ResolvedNucleotideTopologySearchSurface(
            model_name=specification.model_name,
            substitution_parameter_policy="provided",
            substitution_parameter_values=specification.parameter_values,
            substitution_parameter_warnings=[],
            specification=specification,
        )

    optimization_report = optimize_nucleotide_substitution_parameters(
        tree,
        normalized_records,
        model_name=normalized_model_name,
        base_frequencies=base_frequencies,
        initial_kappa=1.0 if normalized_model_name in {"k80", "hky85"} else None,
        initial_exchangeabilities=exchangeabilities,
    )
    specification = resolve_selected_nucleotide_likelihood_specification(
        normalized_records,
        model_name=normalized_model_name,
        owner_name=f"{normalized_model_name.upper()} likelihood NNI search",
        kappa=_resolved_kappa(optimization_report),
        base_frequencies=_resolved_base_frequencies(optimization_report),
        exchangeabilities=_resolved_exchangeabilities(optimization_report),
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return ResolvedNucleotideTopologySearchSurface(
        model_name=specification.model_name,
        substitution_parameter_policy="optimized-on-start-topology",
        substitution_parameter_values=specification.parameter_values,
        substitution_parameter_warnings=list(optimization_report.warnings),
        specification=specification,
    )


def normalize_nucleotide_topology_search_records(
    records: list[AlignmentRecord],
    *,
    owner_name: str,
) -> tuple[list[AlignmentRecord], CompressedAlignmentSitePatterns]:
    """Normalize one nucleotide alignment and precompute compressed site patterns."""
    normalized_records = normalize_unambiguous_dna_records(records, model_name=owner_name)
    return normalized_records, compress_alignment_site_patterns_from_records(
        normalized_records
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


def reoptimize_nucleotide_topology_tree(
    tree: PhyloTree,
    *,
    compressed_patterns: CompressedAlignmentSitePatterns,
    resolved_surface: ResolvedNucleotideTopologySearchSurface,
    branch_reoptimization_policy: str,
    lower_branch_length_bound: float,
    upper_branch_length_bound: float,
    improvement_tolerance: float,
    max_coordinate_passes: int,
) -> BranchReoptimizationResult:
    """Reoptimize one topology-search tree under one resolved likelihood surface."""
    validated_branch_reoptimization_policy = validate_branch_reoptimization_policy(
        branch_reoptimization_policy
    )
    if validated_branch_reoptimization_policy not in _SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES:
        raise ValueError(
            "branch_reoptimization_policy must be one of "
            + ", ".join(sorted(_SUPPORTED_BRANCH_REOPTIMIZATION_POLICIES))
        )
    return optimize_selected_nucleotide_branch_lengths(
        tree,
        compressed_patterns,
        specification=resolved_surface.specification,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
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


def prefer_higher_likelihood(
    left_score: float,
    left_newick: str,
    right_score: float,
    right_newick: str,
) -> bool:
    """Prefer higher likelihoods, then deterministic Newick order, across search moves."""
    if left_score > right_score and not math.isclose(left_score, right_score):
        return True
    if right_score > left_score and not math.isclose(left_score, right_score):
        return False
    return left_newick < right_newick


def _resolved_base_frequencies(
    optimization_report,
) -> dict[str, float] | numpy.ndarray | None:
    if optimization_report.base_frequency_a is None:
        return None
    return {
        "A": optimization_report.base_frequency_a,
        "C": optimization_report.base_frequency_c or 0.0,
        "G": optimization_report.base_frequency_g or 0.0,
        "T": optimization_report.base_frequency_t or 0.0,
    }


def _resolved_kappa(optimization_report) -> float | None:
    row_by_name = {
        row.parameter_name: row.optimized_value
        for row in optimization_report.parameter_rows
    }
    return row_by_name.get("kappa")


def _resolved_exchangeabilities(
    optimization_report,
) -> dict[str, float] | None:
    if optimization_report.model_name != "GTR":
        return None
    row_by_name = {
        row.parameter_name: row.optimized_value
        for row in optimization_report.parameter_rows
    }
    return {
        "AC": optimization_report.fixed_parameter_values.get("AC", 1.0),
        "AG": row_by_name["AG"],
        "AT": row_by_name["AT"],
        "CG": row_by_name["CG"],
        "CT": row_by_name["CT"],
        "GT": row_by_name["GT"],
    }

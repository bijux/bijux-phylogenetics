from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.gamma import (
    build_discrete_gamma_rate_categories,
)
from bijux_phylogenetics.phylo.likelihood.invariant import (
    invariant_component_site_likelihood,
    invariant_proportion_boundary_warnings,
    validate_invariant_proportion,
    validate_invariant_proportion_bounds,
)
from bijux_phylogenetics.phylo.likelihood.logspace import (
    log_weighted_sum_exp,
    logsumexp,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    BranchLengthOptimizationRow,
    DiscreteGammaInvariantMixtureSiteLikelihood,
    DiscreteGammaRateCategory,
    DiscreteGammaSiteLikelihood,
    InvariantMixtureSiteLikelihood,
    ProteinEmpiricalBranchLengthOptimizationReport,
    ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport,
    ProteinEmpiricalDiscreteGammaTreeLikelihoodReport,
    ProteinEmpiricalInvariantMixtureTreeLikelihoodReport,
    ProteinEmpiricalMatrixTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
    run_bounded_likelihood_search,
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
    normalize_unambiguous_protein_records,
    validate_empirical_protein_rate_matrix,
    validate_protein_root_prior,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
    validate_site_log_likelihood_reconstruction,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

_EMPIRICAL_BRANCH_OPTIMIZATION_MODELS = frozenset(
    {
        "fixed-rate",
        "discrete-gamma",
        "invariant",
        "discrete-gamma-invariant",
    }
)


def validate_empirical_branch_optimization_model(likelihood_model: str) -> str:
    if likelihood_model not in _EMPIRICAL_BRANCH_OPTIMIZATION_MODELS:
        raise ValueError(
            "empirical protein branch optimization likelihood_model must be one of "
            f"{sorted(_EMPIRICAL_BRANCH_OPTIMIZATION_MODELS)}"
        )
    return likelihood_model


def evaluate_empirical_protein_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    """Evaluate one fixed-topology protein likelihood from one empirical rate matrix."""
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="empirical protein matrix",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix",
        )
        root_prior_source = "provided"
    return _evaluate_empirical_protein_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        rate_matrix=rate_matrix,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def evaluate_empirical_protein_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    """Evaluate one fixed-topology protein likelihood from paths and one empirical matrix."""
    return evaluate_empirical_protein_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def optimize_empirical_protein_branch_lengths(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    likelihood_model: str,
    alpha: float | None = None,
    invariant_proportion: float | None = None,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> ProteinEmpiricalBranchLengthOptimizationReport:
    """Optimize branch lengths on one fixed topology under one selected empirical protein model."""
    validated_likelihood_model = validate_empirical_branch_optimization_model(
        likelihood_model
    )
    if lower_branch_length_bound < 0.0:
        raise InvalidBranchLengthError(
            "empirical protein branch-length lower bound must be nonnegative"
        )
    if upper_branch_length_bound <= lower_branch_length_bound:
        raise InvalidBranchLengthError(
            "empirical protein branch-length bounds must be strictly increasing"
        )
    if max_coordinate_passes < 1:
        raise ValueError("max_coordinate_passes must be at least one")

    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein branch optimization",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein branch optimization",
        )
        root_prior_source = "provided"
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein branch optimization",
    )
    validated_alpha, validated_invariant_proportion = (
        _validate_empirical_branch_optimization_parameters(
            likelihood_model=validated_likelihood_model,
            alpha=alpha,
            invariant_proportion=invariant_proportion,
            category_count=category_count,
        )
    )
    categories = None
    if validated_alpha is not None:
        categories = build_discrete_gamma_rate_categories(
            alpha=validated_alpha,
            category_count=category_count,
        )

    working_tree = tree.copy()
    validate_explicit_branch_lengths(
        working_tree,
        model_name="empirical protein branch optimization",
    )
    edge_nodes = [child for _parent, child in working_tree.iter_edges()]
    initial_tree_newick = dumps_newick(tree)
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
                "empirical protein branch optimization requires every starting branch length to lie within the declared bounds"
            )
        initial_values[node.node_id] = branch_length
        bounds_by_name[node.node_id] = (
            lower_branch_length_bound,
            upper_branch_length_bound,
        )

    def evaluate_candidate(
        branch_lengths_by_id: dict[str, float],
    ) -> tuple[float, float]:
        _assign_branch_lengths(working_tree, branch_lengths_by_id)
        log_likelihood = _evaluate_empirical_protein_branch_optimization_objective(
            working_tree,
            compressed_patterns,
            validated_rate_matrix=validated_rate_matrix,
            likelihood_model=validated_likelihood_model,
            alpha=validated_alpha,
            invariant_proportion=validated_invariant_proportion,
            categories=categories,
            root_prior=validated_root_prior,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        )
        return log_likelihood, log_likelihood

    initial_log_likelihood, _ = evaluate_candidate(dict(initial_values))
    search_result = run_bounded_coordinate_likelihood_search(
        initial_values=initial_values,
        bounds_by_name=bounds_by_name,
        evaluate=evaluate_candidate,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    optimized_branch_lengths = dict(search_result.parameter_values)
    _assign_branch_lengths(working_tree, optimized_branch_lengths)
    optimized_log_likelihood = float(search_result.objective_value)
    branches = [
        BranchLengthOptimizationRow(
            branch_id=node.node_id or "",
            child_name=node.name,
            descendant_taxa=node.descendant_taxa,
            initial_branch_length=initial_values[node.node_id or ""],
            optimized_branch_length=optimized_branch_lengths[node.node_id or ""],
        )
        for node in edge_nodes
    ]
    return ProteinEmpiricalBranchLengthOptimizationReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        branch_count=len(edge_nodes),
        initial_tree_newick=initial_tree_newick,
        optimized_tree_newick=dumps_newick(working_tree),
        state_count=len(PROTEIN_STATE_ORDER),
        matrix_label=matrix_label,
        root_prior_source=root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        likelihood_model=validated_likelihood_model,
        alpha=validated_alpha,
        invariant_proportion=validated_invariant_proportion,
        initial_log_likelihood=float(initial_log_likelihood),
        optimized_log_likelihood=optimized_log_likelihood,
        optimization_pass_count=search_result.optimization_pass_count,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        branches=branches,
    )


def optimize_empirical_protein_branch_lengths_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    likelihood_model: str,
    alpha: float | None = None,
    invariant_proportion: float | None = None,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    lower_branch_length_bound: float = 0.0,
    upper_branch_length_bound: float = 5.0,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> ProteinEmpiricalBranchLengthOptimizationReport:
    """Optimize one fixed topology from one tree path and one alignment path."""
    return optimize_empirical_protein_branch_lengths(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        likelihood_model=likelihood_model,
        alpha=alpha,
        invariant_proportion=invariant_proportion,
        category_count=category_count,
        root_prior=root_prior,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        lower_branch_length_bound=lower_branch_length_bound,
        upper_branch_length_bound=upper_branch_length_bound,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def evaluate_empirical_protein_tree_likelihood_with_discrete_gamma(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalDiscreteGammaTreeLikelihoodReport:
    """Evaluate one empirical protein likelihood with discrete-gamma categories."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein matrix +G",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix +G",
        )
        root_prior_source = "provided"
    return (
        _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns(
            tree,
            compressed_patterns,
            rate_matrix=rate_matrix,
            alpha=alpha,
            category_count=category_count,
            root_prior=validated_root_prior,
            root_prior_source=root_prior_source,
            matrix_label=matrix_label,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        )
    )


def evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalDiscreteGammaTreeLikelihoodReport:
    """Evaluate one empirical protein +G likelihood from paths."""
    return evaluate_empirical_protein_tree_likelihood_with_discrete_gamma(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        alpha=alpha,
        category_count=category_count,
        root_prior=root_prior,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    initial_invariant_proportion: float = 0.1,
    lower_invariant_proportion_bound: float = 0.0,
    upper_invariant_proportion_bound: float = 0.95,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    """Fit one invariant-site mixture proportion while gamma categories remain active."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein matrix +G+I",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix +G+I",
        )
        root_prior_source = "provided"
    (
        validated_initial_invariant_proportion,
        validated_lower_invariant_proportion_bound,
        validated_upper_invariant_proportion_bound,
    ) = validate_invariant_proportion_bounds(
        initial_invariant_proportion=initial_invariant_proportion,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        model_name="empirical protein matrix +G+I",
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +G+I",
    )
    categories = build_discrete_gamma_rate_categories(
        alpha=alpha,
        category_count=category_count,
    )

    def evaluate_candidate(
        invariant_proportion: float,
    ) -> tuple[ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport, float]:
        report = _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
            tree,
            compressed_patterns,
            validated_rate_matrix=validated_rate_matrix,
            alpha=alpha,
            categories=categories,
            root_prior=validated_root_prior,
            root_prior_source=root_prior_source,
            matrix_label=matrix_label,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            invariant_proportion=invariant_proportion,
            initial_invariant_proportion=validated_initial_invariant_proportion,
            lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
            function_evaluation_count=0,
            converged=False,
            emit_boundary_warnings=False,
        )
        return report, report.log_likelihood

    initial_report, initial_log_likelihood = evaluate_candidate(
        validated_initial_invariant_proportion
    )
    search_result = run_bounded_likelihood_search(
        lower_bound=validated_lower_invariant_proportion_bound,
        upper_bound=validated_upper_invariant_proportion_bound,
        evaluate=evaluate_candidate,
    )
    optimized_report = _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
        tree,
        compressed_patterns,
        validated_rate_matrix=validated_rate_matrix,
        alpha=alpha,
        categories=categories,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        invariant_proportion=search_result.parameter_value,
        initial_invariant_proportion=validated_initial_invariant_proportion,
        lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        emit_boundary_warnings=True,
    )
    return ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport(
        taxa=optimized_report.taxa,
        site_count=optimized_report.site_count,
        pattern_count=optimized_report.pattern_count,
        compression_used=optimized_report.compression_used,
        tree_newick=optimized_report.tree_newick,
        state_count=optimized_report.state_count,
        matrix_label=optimized_report.matrix_label,
        root_prior_source=optimized_report.root_prior_source,
        gap_policy=optimized_report.gap_policy,
        missing_policy=optimized_report.missing_policy,
        alpha=optimized_report.alpha,
        category_count=optimized_report.category_count,
        category_rates=optimized_report.category_rates,
        initial_invariant_proportion=validated_initial_invariant_proportion,
        invariant_proportion=optimized_report.invariant_proportion,
        initial_log_likelihood=initial_log_likelihood,
        log_likelihood=optimized_report.log_likelihood,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
        hit_lower_invariant_proportion_boundary=(
            optimized_report.hit_lower_invariant_proportion_boundary
        ),
        hit_upper_invariant_proportion_boundary=(
            optimized_report.hit_upper_invariant_proportion_boundary
        ),
        boundary_warnings=list(optimized_report.boundary_warnings),
        site_likelihoods=optimized_report.site_likelihoods,
    )


def optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    initial_invariant_proportion: float = 0.1,
    lower_invariant_proportion_bound: float = 0.0,
    upper_invariant_proportion_bound: float = 0.95,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    """Fit one empirical protein +G+I likelihood from one tree path and alignment path."""
    return optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        alpha=alpha,
        category_count=category_count,
        root_prior=root_prior,
        initial_invariant_proportion=initial_invariant_proportion,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one empirical protein likelihood with active gamma and invariant mixture."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein matrix +G+I",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix +G+I",
        )
        root_prior_source = "provided"
    validated_invariant_proportion = validate_invariant_proportion(
        invariant_proportion,
        model_name="empirical protein matrix +G+I",
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +G+I",
    )
    categories = build_discrete_gamma_rate_categories(
        alpha=alpha,
        category_count=category_count,
    )
    return _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
        tree,
        compressed_patterns,
        validated_rate_matrix=validated_rate_matrix,
        alpha=alpha,
        categories=categories,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        invariant_proportion=validated_invariant_proportion,
        initial_invariant_proportion=validated_invariant_proportion,
        lower_invariant_proportion_bound=validated_invariant_proportion,
        upper_invariant_proportion_bound=validated_invariant_proportion,
        function_evaluation_count=1,
        converged=True,
        emit_boundary_warnings=False,
    )


def evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one empirical protein +G+I likelihood from paths."""
    return evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        alpha=alpha,
        invariant_proportion=invariant_proportion,
        category_count=category_count,
        root_prior=root_prior,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def optimize_empirical_protein_tree_likelihood_with_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    initial_invariant_proportion: float = 0.1,
    lower_invariant_proportion_bound: float = 0.0,
    upper_invariant_proportion_bound: float = 0.95,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    """Fit one invariant-site mixture proportion on one empirical protein likelihood."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein matrix +I",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix +I",
        )
        root_prior_source = "provided"
    (
        validated_initial_invariant_proportion,
        validated_lower_invariant_proportion_bound,
        validated_upper_invariant_proportion_bound,
    ) = validate_invariant_proportion_bounds(
        initial_invariant_proportion=initial_invariant_proportion,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        model_name="empirical protein matrix +I",
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +I",
    )
    transition_by_node_id = {
        child.node_id: transition_probability_matrix(
            validated_rate_matrix,
            max(float(child.branch_length or 0.0), 0.0),
        )
        for _parent, child in tree.iter_edges()
    }

    def evaluate_candidate(
        invariant_proportion: float,
    ) -> tuple[ProteinEmpiricalInvariantMixtureTreeLikelihoodReport, float]:
        report = _evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_patterns(
            tree,
            compressed_patterns,
            root_prior=validated_root_prior,
            root_prior_source=root_prior_source,
            matrix_label=matrix_label,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            invariant_proportion=invariant_proportion,
            initial_invariant_proportion=validated_initial_invariant_proportion,
            lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
            transition_by_node_id=transition_by_node_id,
            function_evaluation_count=0,
            converged=False,
        )
        return report, report.log_likelihood

    initial_report, initial_log_likelihood = evaluate_candidate(
        validated_initial_invariant_proportion
    )
    search_result = run_bounded_likelihood_search(
        lower_bound=validated_lower_invariant_proportion_bound,
        upper_bound=validated_upper_invariant_proportion_bound,
        evaluate=evaluate_candidate,
    )
    optimized_report = search_result.payload
    return ProteinEmpiricalInvariantMixtureTreeLikelihoodReport(
        taxa=optimized_report.taxa,
        site_count=optimized_report.site_count,
        pattern_count=optimized_report.pattern_count,
        compression_used=optimized_report.compression_used,
        tree_newick=optimized_report.tree_newick,
        state_count=optimized_report.state_count,
        matrix_label=optimized_report.matrix_label,
        root_prior_source=optimized_report.root_prior_source,
        gap_policy=optimized_report.gap_policy,
        missing_policy=optimized_report.missing_policy,
        initial_invariant_proportion=validated_initial_invariant_proportion,
        invariant_proportion=optimized_report.invariant_proportion,
        initial_log_likelihood=initial_log_likelihood,
        log_likelihood=optimized_report.log_likelihood,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
        site_likelihoods=optimized_report.site_likelihoods,
    )


def optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    initial_invariant_proportion: float = 0.1,
    lower_invariant_proportion_bound: float = 0.0,
    upper_invariant_proportion_bound: float = 0.95,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    """Fit one invariant-site mixture from one tree path and one alignment path."""
    return optimize_empirical_protein_tree_likelihood_with_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        initial_invariant_proportion=initial_invariant_proportion,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def evaluate_empirical_protein_tree_likelihood_with_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    invariant_proportion: float,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    """Evaluate one empirical protein likelihood under one fixed invariant-site mixture."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein matrix +I",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix +I",
        )
        root_prior_source = "provided"
    validated_invariant_proportion = validate_invariant_proportion(
        invariant_proportion,
        model_name="empirical protein matrix +I",
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +I",
    )
    transition_by_node_id = {
        child.node_id: transition_probability_matrix(
            validated_rate_matrix,
            max(float(child.branch_length or 0.0), 0.0),
        )
        for _parent, child in tree.iter_edges()
    }
    return _evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_patterns(
        tree,
        compressed_patterns,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        invariant_proportion=validated_invariant_proportion,
        initial_invariant_proportion=validated_invariant_proportion,
        lower_invariant_proportion_bound=validated_invariant_proportion,
        upper_invariant_proportion_bound=validated_invariant_proportion,
        transition_by_node_id=transition_by_node_id,
        function_evaluation_count=1,
        converged=True,
    )


def evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    invariant_proportion: float,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    """Evaluate one empirical protein +I likelihood from paths."""
    return evaluate_empirical_protein_tree_likelihood_with_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        rate_matrix=rate_matrix,
        invariant_proportion=invariant_proportion,
        root_prior=root_prior,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def _evaluate_empirical_protein_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray,
    root_prior_source: str,
    matrix_label: str,
    gap_policy: str,
    missing_policy: str,
    ambiguity_policy: str,
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix",
    )
    transition_by_node_id = {
        child.node_id: transition_probability_matrix(
            validated_rate_matrix,
            max(float(child.branch_length or 0.0), 0.0),
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
                    model_name="empirical protein matrix",
                    root_prior=root_prior,
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
        owner_name="empirical protein matrix likelihood",
    )
    return ProteinEmpiricalMatrixTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        state_count=len(PROTEIN_STATE_ORDER),
        matrix_label=matrix_label,
        root_prior_source=root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
        log_likelihood=log_likelihood,
        site_log_likelihoods=site_log_likelihoods,
    )


def _assign_branch_lengths(
    tree: PhyloTree,
    branch_lengths_by_id: dict[str, float],
) -> None:
    for _parent, child in tree.iter_edges():
        if child.node_id is None:
            raise ValueError("tree node is missing a stable node_id")
        child.branch_length = branch_lengths_by_id[child.node_id]


def _validate_empirical_branch_optimization_parameters(
    *,
    likelihood_model: str,
    alpha: float | None,
    invariant_proportion: float | None,
    category_count: int,
) -> tuple[float | None, float | None]:
    validated_alpha: float | None = None
    validated_invariant_proportion: float | None = None
    if likelihood_model in {"discrete-gamma", "discrete-gamma-invariant"}:
        if alpha is None:
            raise ValueError(
                f"likelihood_model '{likelihood_model}' requires one alpha value"
            )
        build_discrete_gamma_rate_categories(
            alpha=alpha,
            category_count=category_count,
        )
        validated_alpha = float(alpha)
    elif alpha is not None:
        raise ValueError(f"likelihood_model '{likelihood_model}' does not accept alpha")
    if likelihood_model in {"invariant", "discrete-gamma-invariant"}:
        if invariant_proportion is None:
            raise ValueError(
                f"likelihood_model '{likelihood_model}' requires one invariant_proportion value"
            )
        validated_invariant_proportion = validate_invariant_proportion(
            invariant_proportion,
            model_name="empirical protein branch optimization",
        )
    elif invariant_proportion is not None:
        raise ValueError(
            f"likelihood_model '{likelihood_model}' does not accept invariant_proportion"
        )
    return validated_alpha, validated_invariant_proportion


def _evaluate_empirical_protein_branch_optimization_objective(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    validated_rate_matrix: numpy.ndarray,
    likelihood_model: str,
    alpha: float | None,
    invariant_proportion: float | None,
    categories: list[DiscreteGammaRateCategory] | None,
    root_prior: numpy.ndarray,
    gap_policy: str,
    missing_policy: str,
) -> float:
    if likelihood_model == "fixed-rate":
        return _evaluate_empirical_protein_tree_likelihood_from_patterns(
            tree,
            compressed_patterns,
            rate_matrix=validated_rate_matrix,
            root_prior=root_prior,
            root_prior_source="provided",
            matrix_label="empirical",
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy="reject",
        ).log_likelihood
    if likelihood_model == "discrete-gamma":
        if alpha is None:
            raise ValueError("discrete-gamma objective requires alpha")
        return _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns(
            tree,
            compressed_patterns,
            rate_matrix=validated_rate_matrix,
            alpha=alpha,
            category_count=len(categories or []),
            root_prior=root_prior,
            root_prior_source="provided",
            matrix_label="empirical",
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        ).log_likelihood
    if likelihood_model == "invariant":
        if invariant_proportion is None:
            raise ValueError("invariant objective requires invariant_proportion")
        transition_by_node_id = _empirical_transition_by_node_id(
            tree,
            validated_rate_matrix,
        )
        return _evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_patterns(
            tree,
            compressed_patterns,
            root_prior=root_prior,
            root_prior_source="provided",
            matrix_label="empirical",
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            invariant_proportion=invariant_proportion,
            initial_invariant_proportion=invariant_proportion,
            lower_invariant_proportion_bound=invariant_proportion,
            upper_invariant_proportion_bound=invariant_proportion,
            transition_by_node_id=transition_by_node_id,
            function_evaluation_count=1,
            converged=True,
        ).log_likelihood
    if likelihood_model == "discrete-gamma-invariant":
        if alpha is None or invariant_proportion is None or categories is None:
            raise ValueError(
                "discrete-gamma-invariant objective requires alpha and invariant proportion"
            )
        return _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
            tree,
            compressed_patterns,
            validated_rate_matrix=validated_rate_matrix,
            alpha=alpha,
            categories=categories,
            root_prior=root_prior,
            root_prior_source="provided",
            matrix_label="empirical",
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            invariant_proportion=invariant_proportion,
            initial_invariant_proportion=invariant_proportion,
            lower_invariant_proportion_bound=invariant_proportion,
            upper_invariant_proportion_bound=invariant_proportion,
            function_evaluation_count=1,
            converged=True,
            emit_boundary_warnings=False,
        ).log_likelihood
    raise ValueError(
        f"unsupported empirical branch optimization model '{likelihood_model}'"
    )


def _empirical_transition_by_node_id(
    tree: PhyloTree,
    validated_rate_matrix: numpy.ndarray,
    *,
    rate_scale: float = 1.0,
) -> dict[str | None, numpy.ndarray]:
    return {
        child.node_id: transition_probability_matrix(
            validated_rate_matrix,
            max(float(child.branch_length or 0.0), 0.0) * rate_scale,
        )
        for _parent, child in tree.iter_edges()
    }


def _empirical_protein_site_log_likelihood(
    tree: PhyloTree,
    pattern_states: tuple[str, ...],
    *,
    taxon_order: list[str],
    model_name: str,
    root_prior: numpy.ndarray,
    transition_by_node_id: dict[str | None, numpy.ndarray],
    gap_policy: str,
    missing_policy: str,
) -> float:
    return evaluate_fixed_topology_protein_site_log_likelihood(
        tree,
        pattern_states,
        taxon_order=taxon_order,
        model_name=model_name,
        root_prior=root_prior,
        transition_matrix_for_child=lambda child: transition_by_node_id[
            child.node_id or ""
        ],
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def _empirical_protein_discrete_gamma_category_log_likelihoods(
    tree: PhyloTree,
    pattern_states: tuple[str, ...],
    *,
    taxon_order: list[str],
    validated_rate_matrix: numpy.ndarray,
    categories: list[DiscreteGammaRateCategory],
    model_name: str,
    root_prior: numpy.ndarray,
    gap_policy: str,
    missing_policy: str,
) -> list[float]:
    category_log_likelihoods: list[float] = []
    for category in categories:
        transition_by_node_id = _empirical_transition_by_node_id(
            tree,
            validated_rate_matrix,
            rate_scale=category.rate,
        )
        category_log_likelihoods.append(
            _empirical_protein_site_log_likelihood(
                tree,
                pattern_states,
                taxon_order=taxon_order,
                model_name=model_name,
                root_prior=root_prior,
                transition_by_node_id=transition_by_node_id,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
            )
        )
    return category_log_likelihoods


def _discrete_gamma_mixture_log_likelihood(
    categories: list[DiscreteGammaRateCategory],
    category_log_likelihoods: list[float],
) -> float:
    return log_weighted_sum_exp(
        category_log_likelihoods,
        weights=[
            category.weight
            for category in categories
        ],
    )


def _linear_likelihoods_from_log_values(log_values: list[float]) -> list[float]:
    return [math.exp(log_value) for log_value in log_values]


def _linear_likelihood_from_log_value(log_value: float) -> float:
    return math.exp(log_value)


def _invariant_mixture_log_likelihood(
    *,
    invariant_proportion: float,
    invariant_component_likelihood: float,
    variable_component_log_likelihood: float,
) -> float:
    if invariant_proportion <= 0.0:
        return variable_component_log_likelihood
    if invariant_proportion >= 1.0:
        if invariant_component_likelihood <= 0.0:
            return float("-inf")
        return math.log(invariant_component_likelihood)
    invariant_log_likelihood = (
        float("-inf")
        if invariant_component_likelihood <= 0.0
        else math.log(invariant_component_likelihood)
    )
    return logsumexp(
        (
            math.log(invariant_proportion) + invariant_log_likelihood,
            math.log(1.0 - invariant_proportion) + variable_component_log_likelihood,
        )
    )


def _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    category_count: int,
    root_prior: numpy.ndarray,
    root_prior_source: str,
    matrix_label: str,
    gap_policy: str,
    missing_policy: str,
) -> ProteinEmpiricalDiscreteGammaTreeLikelihoodReport:
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +G",
    )
    categories = build_discrete_gamma_rate_categories(
        alpha=alpha,
        category_count=category_count,
    )
    site_likelihoods: list[DiscreteGammaSiteLikelihood] = []
    total_log_likelihood = 0.0
    for pattern in compressed_patterns.patterns:
        category_log_likelihoods = (
            _empirical_protein_discrete_gamma_category_log_likelihoods(
                tree,
                pattern.states,
                taxon_order=compressed_patterns.taxon_order,
                validated_rate_matrix=validated_rate_matrix,
                categories=categories,
                model_name="empirical protein matrix +G",
                root_prior=root_prior,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
            )
        )
        pattern_log_likelihood = _discrete_gamma_mixture_log_likelihood(
            categories,
            category_log_likelihoods,
        )
        category_likelihoods = _linear_likelihoods_from_log_values(
            category_log_likelihoods
        )
        mixture_likelihood = _linear_likelihood_from_log_value(pattern_log_likelihood)
        for site_position in pattern.site_positions:
            site_likelihoods.append(
                DiscreteGammaSiteLikelihood(
                    pattern_id=pattern.pattern_id,
                    pattern_weight=pattern.weight,
                    site_position=site_position,
                    site_states=pattern.states,
                    category_likelihoods=list(category_likelihoods),
                    mixture_likelihood=mixture_likelihood,
                    log_likelihood=pattern_log_likelihood,
                )
            )
        total_log_likelihood += pattern.weight * pattern_log_likelihood
    validate_site_log_likelihood_reconstruction(
        site_likelihoods,
        expected_total_log_likelihood=total_log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name="empirical protein matrix +G likelihood",
    )
    return ProteinEmpiricalDiscreteGammaTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        state_count=len(PROTEIN_STATE_ORDER),
        matrix_label=matrix_label,
        root_prior_source=root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        alpha=alpha,
        category_count=category_count,
        category_rates=list(categories),
        site_likelihoods=site_likelihoods,
        log_likelihood=total_log_likelihood,
    )


def _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    validated_rate_matrix: numpy.ndarray,
    alpha: float,
    categories: list[DiscreteGammaRateCategory],
    root_prior: numpy.ndarray,
    root_prior_source: str,
    matrix_label: str,
    gap_policy: str,
    missing_policy: str,
    invariant_proportion: float,
    initial_invariant_proportion: float,
    lower_invariant_proportion_bound: float,
    upper_invariant_proportion_bound: float,
    function_evaluation_count: int,
    converged: bool,
    emit_boundary_warnings: bool,
) -> ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    site_likelihoods: list[DiscreteGammaInvariantMixtureSiteLikelihood] = []
    total_log_likelihood = 0.0
    for pattern in compressed_patterns.patterns:
        category_log_likelihoods = (
            _empirical_protein_discrete_gamma_category_log_likelihoods(
                tree,
                pattern.states,
                taxon_order=compressed_patterns.taxon_order,
                validated_rate_matrix=validated_rate_matrix,
                categories=categories,
                model_name="empirical protein matrix +G+I",
                root_prior=root_prior,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
            )
        )
        variable_component_log_likelihood = _discrete_gamma_mixture_log_likelihood(
            categories,
            category_log_likelihoods,
        )
        variable_component_likelihood = _linear_likelihood_from_log_value(
            variable_component_log_likelihood
        )
        invariant_component_likelihood = invariant_component_site_likelihood(
            pattern.states,
            root_prior=root_prior,
            model_name="empirical protein matrix +G+I",
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        )
        pattern_log_likelihood = _invariant_mixture_log_likelihood(
            invariant_proportion=invariant_proportion,
            invariant_component_likelihood=invariant_component_likelihood,
            variable_component_log_likelihood=variable_component_log_likelihood,
        )
        mixture_likelihood = _linear_likelihood_from_log_value(pattern_log_likelihood)
        category_likelihoods = _linear_likelihoods_from_log_values(
            category_log_likelihoods
        )
        for site_position in pattern.site_positions:
            site_likelihoods.append(
                DiscreteGammaInvariantMixtureSiteLikelihood(
                    pattern_id=pattern.pattern_id,
                    pattern_weight=pattern.weight,
                    site_position=site_position,
                    site_states=pattern.states,
                    category_likelihoods=list(category_likelihoods),
                    invariant_component_likelihood=invariant_component_likelihood,
                    variable_component_likelihood=variable_component_likelihood,
                    mixture_likelihood=mixture_likelihood,
                    log_likelihood=pattern_log_likelihood,
                )
            )
        total_log_likelihood += pattern.weight * pattern_log_likelihood
    validate_site_log_likelihood_reconstruction(
        site_likelihoods,
        expected_total_log_likelihood=total_log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name="empirical protein matrix +G+I likelihood",
    )
    if emit_boundary_warnings:
        (
            hit_lower_invariant_proportion_boundary,
            hit_upper_invariant_proportion_boundary,
            boundary_warnings,
        ) = invariant_proportion_boundary_warnings(
            invariant_proportion=invariant_proportion,
            lower_invariant_proportion_bound=lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        )
    else:
        hit_lower_invariant_proportion_boundary = False
        hit_upper_invariant_proportion_boundary = False
        boundary_warnings = []
    return ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        state_count=len(PROTEIN_STATE_ORDER),
        matrix_label=matrix_label,
        root_prior_source=root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        alpha=alpha,
        category_count=len(categories),
        category_rates=list(categories),
        initial_invariant_proportion=initial_invariant_proportion,
        invariant_proportion=invariant_proportion,
        initial_log_likelihood=total_log_likelihood,
        log_likelihood=total_log_likelihood,
        function_evaluation_count=function_evaluation_count,
        converged=converged,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        hit_lower_invariant_proportion_boundary=hit_lower_invariant_proportion_boundary,
        hit_upper_invariant_proportion_boundary=hit_upper_invariant_proportion_boundary,
        boundary_warnings=boundary_warnings,
        site_likelihoods=site_likelihoods,
    )


def _evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    root_prior: numpy.ndarray,
    root_prior_source: str,
    matrix_label: str,
    gap_policy: str,
    missing_policy: str,
    invariant_proportion: float,
    initial_invariant_proportion: float,
    lower_invariant_proportion_bound: float,
    upper_invariant_proportion_bound: float,
    transition_by_node_id: dict[str | None, numpy.ndarray],
    function_evaluation_count: int,
    converged: bool,
) -> ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    site_likelihoods: list[InvariantMixtureSiteLikelihood] = []
    total_log_likelihood = 0.0
    for pattern in compressed_patterns.patterns:
        variable_component_log_likelihood = _empirical_protein_site_log_likelihood(
            tree,
            pattern.states,
            taxon_order=compressed_patterns.taxon_order,
            model_name="empirical protein matrix +I",
            root_prior=root_prior,
            transition_by_node_id=transition_by_node_id,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        )
        variable_component_likelihood = _linear_likelihood_from_log_value(
            variable_component_log_likelihood
        )
        invariant_component_likelihood = invariant_component_site_likelihood(
            pattern.states,
            root_prior=root_prior,
            model_name="empirical protein matrix +I",
            gap_policy=gap_policy,
            missing_policy=missing_policy,
        )
        pattern_log_likelihood = _invariant_mixture_log_likelihood(
            invariant_proportion=invariant_proportion,
            invariant_component_likelihood=invariant_component_likelihood,
            variable_component_log_likelihood=variable_component_log_likelihood,
        )
        mixture_likelihood = _linear_likelihood_from_log_value(pattern_log_likelihood)
        for site_position in pattern.site_positions:
            site_likelihoods.append(
                InvariantMixtureSiteLikelihood(
                    pattern_id=pattern.pattern_id,
                    pattern_weight=pattern.weight,
                    site_position=site_position,
                    site_states=pattern.states,
                    invariant_component_likelihood=invariant_component_likelihood,
                    variable_component_likelihood=variable_component_likelihood,
                    mixture_likelihood=mixture_likelihood,
                    log_likelihood=pattern_log_likelihood,
                )
            )
        total_log_likelihood += pattern.weight * pattern_log_likelihood
    validate_site_log_likelihood_reconstruction(
        site_likelihoods,
        expected_total_log_likelihood=total_log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name="empirical protein matrix +I likelihood",
    )
    return ProteinEmpiricalInvariantMixtureTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        state_count=len(PROTEIN_STATE_ORDER),
        matrix_label=matrix_label,
        root_prior_source=root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        initial_invariant_proportion=initial_invariant_proportion,
        invariant_proportion=invariant_proportion,
        initial_log_likelihood=total_log_likelihood,
        log_likelihood=total_log_likelihood,
        function_evaluation_count=function_evaluation_count,
        converged=converged,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        site_likelihoods=site_likelihoods,
    )

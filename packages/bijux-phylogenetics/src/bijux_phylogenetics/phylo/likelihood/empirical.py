from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    solve_ctmc_stationary_distribution,
)
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
    FiniteStateDiscreteGammaInvariantTreeLikelihoodReport,
    FiniteStateInvariantMixtureTreeLikelihoodReport,
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
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    compile_finite_state_discrete_gamma_contexts,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns,
    evaluate_finite_state_tree_likelihood_with_invariant_mixture,
    evaluate_finite_state_tree_likelihood_with_invariant_mixture_from_patterns,
    evaluate_finite_state_site_log_likelihood_rows_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    evaluate_fixed_topology_protein_site_log_likelihood,
    normalize_protein_likelihood_records,
    normalize_unambiguous_protein_records,
    protein_leaf_likelihood_vector,
    resolve_protein_root_prior,
    validate_empirical_protein_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.protein_empirical_models import (
    empirical_protein_rate_matrix,
    validate_empirical_protein_model_name,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_empirical_protein_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.sites import (
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


def _stationary_distribution_for_empirical_rate_matrix(
    rate_matrix: numpy.ndarray,
) -> numpy.ndarray:
    report = solve_ctmc_stationary_distribution(
        rate_matrix,
        state_labels=PROTEIN_STATE_ORDER,
    )
    return numpy.asarray(report.probabilities, dtype=float)


def _resolve_empirical_protein_root_prior(
    records: list[AlignmentRecord],
    *,
    owner_name: str,
    validated_rate_matrix: numpy.ndarray,
    default_policy: str = "equal",
    root_prior_policy: str | None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None,
    fixed_root_state: str | None,
) -> tuple[numpy.ndarray, str]:
    resolved = resolve_protein_root_prior(
        records,
        owner_name=owner_name,
        default_policy=default_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=_stationary_distribution_for_empirical_rate_matrix(
            validated_rate_matrix
        ),
    )
    return resolved.root_prior, resolved.root_prior_source


def evaluate_empirical_protein_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix",
    )
    validated_root_prior, root_prior_source = _resolve_empirical_protein_root_prior(
        normalized_records,
        owner_name="empirical protein matrix",
        validated_rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return _evaluate_empirical_protein_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        rate_matrix=validated_rate_matrix,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        matrix_label=matrix_label,
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def evaluate_empirical_protein_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def evaluate_named_empirical_protein_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    """Evaluate one fixed-topology protein likelihood from one named empirical model."""
    validated_model_name = validate_empirical_protein_model_name(model_name)
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name=f"{validated_model_name} empirical protein matrix",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        empirical_protein_rate_matrix(validated_model_name),
        model_name=f"{validated_model_name} empirical protein matrix",
    )
    validated_root_prior, root_prior_source = _resolve_empirical_protein_root_prior(
        normalized_records,
        owner_name=f"{validated_model_name} empirical protein matrix",
        validated_rate_matrix=validated_rate_matrix,
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return _evaluate_empirical_protein_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        rate_matrix=validated_rate_matrix,
        root_prior=validated_root_prior,
        root_prior_source=root_prior_source,
        matrix_label=validated_model_name,
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )


def evaluate_named_empirical_protein_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    """Evaluate one fixed-topology protein likelihood from paths and one named empirical model."""
    return evaluate_named_empirical_protein_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein branch optimization",
    )
    validated_root_prior, root_prior_source = _resolve_empirical_protein_root_prior(
        normalized_records,
        owner_name="empirical protein branch optimization",
        validated_rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
        owner_name="empirical protein branch-length fit",
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
        parameter_diagnostics=search_result.parameter_diagnostics,
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
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +G",
    )
    validated_root_prior, root_prior_source = _resolve_empirical_protein_root_prior(
        normalized_records,
        owner_name="empirical protein matrix +G",
        validated_rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    return (
        _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns(
            tree,
            compressed_patterns,
            rate_matrix=validated_rate_matrix,
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
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    surface = build_empirical_protein_finite_state_surface(
        normalized_records,
        rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    categories = build_discrete_gamma_rate_categories(
        alpha=alpha,
        category_count=category_count,
    )
    contexts = compile_finite_state_discrete_gamma_contexts(
        tree,
        surface,
        categories=categories,
    )

    def evaluate_candidate(
        invariant_proportion: float,
    ) -> tuple[ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport, float]:
        finite_state_report = (
            evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
                contexts,
                compressed_patterns,
                surface=surface,
                alpha=alpha,
                category_rates=categories,
                invariant_proportion=invariant_proportion,
                lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
                upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
                emit_boundary_warnings=False,
            )
        )
        report = _protein_empirical_discrete_gamma_invariant_report_from_finite_state(
            finite_state_report,
            alpha=alpha,
            matrix_label=matrix_label,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            initial_invariant_proportion=validated_initial_invariant_proportion,
            lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
            function_evaluation_count=0,
            converged=False,
            initial_log_likelihood=finite_state_report.log_likelihood,
        )
        return report, report.log_likelihood

    _initial_report, initial_log_likelihood = evaluate_candidate(
        validated_initial_invariant_proportion
    )
    search_result = run_bounded_likelihood_search(
        lower_bound=validated_lower_invariant_proportion_bound,
        upper_bound=validated_upper_invariant_proportion_bound,
        evaluate=evaluate_candidate,
        parameter_name="invariant_proportion",
        owner_name="empirical protein gamma-invariant fit",
    )
    optimized_finite_state_report = (
        evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_patterns(
            contexts,
            compressed_patterns,
            surface=surface,
            alpha=alpha,
            category_rates=categories,
            invariant_proportion=search_result.parameter_value,
            lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
            emit_boundary_warnings=True,
        )
    )
    optimized_report = _protein_empirical_discrete_gamma_invariant_report_from_finite_state(
        optimized_finite_state_report,
        alpha=alpha,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        initial_invariant_proportion=validated_initial_invariant_proportion,
        lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        initial_log_likelihood=initial_log_likelihood,
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
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
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
        boundary_warning_rows=list(optimized_report.boundary_warning_rows),
        optimizer_diagnostics=search_result.diagnostics,
        site_likelihoods=optimized_report.site_likelihoods,
    )


def optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    category_count: int = 4,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    validated_invariant_proportion = validate_invariant_proportion(
        invariant_proportion,
        model_name="empirical protein matrix +G+I",
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +G+I",
    )
    surface = build_empirical_protein_finite_state_surface(
        normalized_records,
        rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    finite_state_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=alpha,
        category_count=category_count,
        invariant_proportion=validated_invariant_proportion,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        emit_boundary_warnings=emit_boundary_warnings,
    )
    return _protein_empirical_discrete_gamma_invariant_report_from_finite_state(
        finite_state_report,
        alpha=alpha,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        initial_invariant_proportion=validated_invariant_proportion,
        lower_invariant_proportion_bound=(
            validated_invariant_proportion
            if lower_invariant_proportion_bound is None
            else lower_invariant_proportion_bound
        ),
        upper_invariant_proportion_bound=(
            validated_invariant_proportion
            if upper_invariant_proportion_bound is None
            else upper_invariant_proportion_bound
        ),
        function_evaluation_count=1,
        converged=True,
        initial_log_likelihood=finite_state_report.log_likelihood,
    )


def evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        emit_boundary_warnings=emit_boundary_warnings,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def optimize_empirical_protein_tree_likelihood_with_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    surface = build_empirical_protein_finite_state_surface(
        normalized_records,
        rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    def evaluate_candidate(
        invariant_proportion: float,
    ) -> tuple[ProteinEmpiricalInvariantMixtureTreeLikelihoodReport, float]:
        finite_state_report = (
            evaluate_finite_state_tree_likelihood_with_invariant_mixture_from_patterns(
                context,
                compressed_patterns,
                invariant_proportion=invariant_proportion,
            )
        )
        report = _protein_empirical_invariant_mixture_report_from_finite_state(
            finite_state_report,
            matrix_label=matrix_label,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            initial_invariant_proportion=validated_initial_invariant_proportion,
            lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
            function_evaluation_count=0,
            converged=False,
            initial_log_likelihood=finite_state_report.log_likelihood,
        )
        return report, report.log_likelihood

    initial_report, initial_log_likelihood = evaluate_candidate(
        validated_initial_invariant_proportion
    )
    search_result = run_bounded_likelihood_search(
        lower_bound=validated_lower_invariant_proportion_bound,
        upper_bound=validated_upper_invariant_proportion_bound,
        evaluate=evaluate_candidate,
        parameter_name="invariant_proportion",
        owner_name="empirical protein invariant-proportion fit",
    )
    optimized_finite_state_report = (
        evaluate_finite_state_tree_likelihood_with_invariant_mixture_from_patterns(
            context,
            compressed_patterns,
            invariant_proportion=search_result.parameter_value,
            lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
            emit_boundary_warnings=True,
        )
    )
    optimized_report = _protein_empirical_invariant_mixture_report_from_finite_state(
        optimized_finite_state_report,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        initial_invariant_proportion=validated_initial_invariant_proportion,
        lower_invariant_proportion_bound=validated_lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=validated_upper_invariant_proportion_bound,
        function_evaluation_count=search_result.function_evaluation_count,
        converged=search_result.converged,
        initial_log_likelihood=initial_log_likelihood,
    )
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
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
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
        boundary_warning_rows=list(optimized_report.boundary_warning_rows),
        optimizer_diagnostics=search_result.diagnostics,
        site_likelihoods=optimized_report.site_likelihoods,
    )


def optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
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
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
    validated_invariant_proportion = validate_invariant_proportion(
        invariant_proportion,
        model_name="empirical protein matrix +I",
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix +I",
    )
    surface = build_empirical_protein_finite_state_surface(
        normalized_records,
        rate_matrix=validated_rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    finite_state_report = evaluate_finite_state_tree_likelihood_with_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        invariant_proportion=validated_invariant_proportion,
    )
    return _protein_empirical_invariant_mixture_report_from_finite_state(
        finite_state_report,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        initial_invariant_proportion=validated_invariant_proportion,
        lower_invariant_proportion_bound=validated_invariant_proportion,
        upper_invariant_proportion_bound=validated_invariant_proportion,
        function_evaluation_count=1,
        converged=True,
        initial_log_likelihood=finite_state_report.log_likelihood,
    )


def evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    invariant_proportion: float,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
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
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )


def _protein_empirical_invariant_mixture_report_from_finite_state(
    report: FiniteStateInvariantMixtureTreeLikelihoodReport,
    *,
    matrix_label: str,
    gap_policy: str,
    missing_policy: str,
    initial_invariant_proportion: float,
    lower_invariant_proportion_bound: float,
    upper_invariant_proportion_bound: float,
    function_evaluation_count: int,
    converged: bool,
    initial_log_likelihood: float,
) -> ProteinEmpiricalInvariantMixtureTreeLikelihoodReport:
    normalization_policy = (
        CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
        if report.normalization_policy is None
        else report.normalization_policy
    )
    return ProteinEmpiricalInvariantMixtureTreeLikelihoodReport(
        taxa=report.taxa,
        site_count=report.site_count,
        pattern_count=report.pattern_count,
        compression_used=report.compression_used,
        tree_newick=report.tree_newick,
        state_count=report.state_count,
        matrix_label=matrix_label,
        root_prior_source=report.root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        normalization_policy=normalization_policy,
        initial_invariant_proportion=initial_invariant_proportion,
        invariant_proportion=report.invariant_proportion,
        initial_log_likelihood=initial_log_likelihood,
        log_likelihood=report.log_likelihood,
        function_evaluation_count=function_evaluation_count,
        converged=converged,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        hit_lower_invariant_proportion_boundary=(
            report.hit_lower_invariant_proportion_boundary
        ),
        hit_upper_invariant_proportion_boundary=(
            report.hit_upper_invariant_proportion_boundary
        ),
        boundary_warnings=list(report.boundary_warnings),
        boundary_warning_rows=list(report.boundary_warning_rows),
        optimizer_diagnostics=None,
        site_likelihoods=report.site_likelihoods,
    )


def _protein_empirical_discrete_gamma_invariant_report_from_finite_state(
    report: FiniteStateDiscreteGammaInvariantTreeLikelihoodReport,
    *,
    alpha: float,
    matrix_label: str,
    gap_policy: str,
    missing_policy: str,
    initial_invariant_proportion: float,
    lower_invariant_proportion_bound: float,
    upper_invariant_proportion_bound: float,
    function_evaluation_count: int,
    converged: bool,
    initial_log_likelihood: float,
) -> ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport:
    normalization_policy = (
        CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY
        if report.normalization_policy is None
        else report.normalization_policy
    )
    return ProteinEmpiricalDiscreteGammaInvariantTreeLikelihoodReport(
        taxa=report.taxa,
        site_count=report.site_count,
        pattern_count=report.pattern_count,
        compression_used=report.compression_used,
        tree_newick=report.tree_newick,
        state_count=report.state_count,
        matrix_label=matrix_label,
        root_prior_source=report.root_prior_source,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        normalization_policy=normalization_policy,
        alpha=alpha,
        category_count=report.category_count,
        category_rates=report.category_rates,
        initial_invariant_proportion=initial_invariant_proportion,
        invariant_proportion=report.invariant_proportion,
        initial_log_likelihood=initial_log_likelihood,
        log_likelihood=report.log_likelihood,
        function_evaluation_count=function_evaluation_count,
        converged=converged,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        hit_lower_invariant_proportion_boundary=report.hit_lower_invariant_proportion_boundary,
        hit_upper_invariant_proportion_boundary=report.hit_upper_invariant_proportion_boundary,
        boundary_warnings=list(report.boundary_warnings),
        boundary_warning_rows=list(report.boundary_warning_rows),
        optimizer_diagnostics=None,
        site_likelihoods=report.site_likelihoods,
    )


def _evaluate_empirical_protein_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray,
    root_prior_source: str,
    matrix_label: str,
    normalization_policy: str,
    gap_policy: str,
    missing_policy: str,
    ambiguity_policy: str,
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix",
    )
    context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="empirical protein matrix",
            state_labels=PROTEIN_STATE_ORDER,
            root_prior=root_prior,
            root_prior_source=root_prior_source,
            observation_policy=(
                f"gap={gap_policy};missing={missing_policy};ambiguity={ambiguity_policy}"
            ),
            normalization_policy=normalization_policy,
            transition_rate_matrix=validated_rate_matrix,
            transition_matrix_for_branch_length=lambda branch_length: (
                transition_probability_matrix(validated_rate_matrix, branch_length)
            ),
            leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
                protein_leaf_likelihood_vector(
                    {taxon_name: observed_state},
                    model_name="empirical protein matrix",
                    node_name=taxon_name,
                    gap_policy=gap_policy,
                    missing_policy=missing_policy,
                    ambiguity_policy=ambiguity_policy,
                )
            ),
        ),
    )
    site_log_likelihoods, log_likelihood = (
        evaluate_finite_state_site_log_likelihood_rows_from_patterns(
            context,
            compressed_patterns,
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
        normalization_policy=normalization_policy,
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
            normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
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
        weights=[category.weight for category in categories],
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
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
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
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
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
        optimizer_diagnostics=None,
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
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
        initial_invariant_proportion=initial_invariant_proportion,
        invariant_proportion=invariant_proportion,
        initial_log_likelihood=total_log_likelihood,
        log_likelihood=total_log_likelihood,
        function_evaluation_count=function_evaluation_count,
        converged=converged,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        optimizer_diagnostics=None,
        site_likelihoods=site_likelihoods,
    )

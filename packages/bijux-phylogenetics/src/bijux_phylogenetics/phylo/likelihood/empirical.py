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
from bijux_phylogenetics.phylo.likelihood.models import (
    DiscreteGammaSiteLikelihood,
    ProteinEmpiricalDiscreteGammaTreeLikelihoodReport,
    ProteinEmpiricalMatrixTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    UNIFORM_PROTEIN_ROOT_PRIOR,
    evaluate_fixed_topology_protein_likelihood_from_patterns,
    evaluate_fixed_topology_protein_site_log_likelihood,
    normalize_unambiguous_protein_records,
    validate_empirical_protein_rate_matrix,
    validate_protein_root_prior,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def evaluate_empirical_protein_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    rate_matrix: numpy.ndarray | list[list[float]] | tuple[tuple[float, ...], ...],
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    matrix_label: str = "empirical",
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
) -> ProteinEmpiricalMatrixTreeLikelihoodReport:
    """Evaluate one fixed-topology protein likelihood from one empirical rate matrix."""
    normalized_records = normalize_unambiguous_protein_records(
        records,
        model_name="empirical protein matrix",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
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
    compressed_patterns = compress_alignment_site_patterns_from_records(normalized_records)
    if root_prior is None:
        validated_root_prior = UNIFORM_PROTEIN_ROOT_PRIOR
        root_prior_source = "uniform"
    else:
        validated_root_prior = validate_protein_root_prior(
            root_prior,
            model_name="empirical protein matrix +G",
        )
        root_prior_source = "provided"
    return _evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_patterns(
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
    log_likelihood = evaluate_fixed_topology_protein_likelihood_from_patterns(
        tree,
        compressed_patterns,
        model_name="empirical protein matrix",
        root_prior=root_prior,
        transition_matrix_for_child=lambda child: transition_by_node_id[
            child.node_id or ""
        ],
        gap_policy=gap_policy,
        missing_policy=missing_policy,
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
        log_likelihood=log_likelihood,
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
        category_likelihoods: list[float] = []
        for category in categories:
            transition_by_node_id = {
                child.node_id: transition_probability_matrix(
                    validated_rate_matrix,
                    max(float(child.branch_length or 0.0), 0.0) * category.rate,
                )
                for _parent, child in tree.iter_edges()
            }
            site_log_likelihood = evaluate_fixed_topology_protein_site_log_likelihood(
                tree,
                pattern.states,
                taxon_order=compressed_patterns.taxon_order,
                model_name="empirical protein matrix +G",
                root_prior=root_prior,
                transition_matrix_for_child=lambda child: transition_by_node_id[
                    child.node_id or ""
                ],
                gap_policy=gap_policy,
                missing_policy=missing_policy,
            )
            category_likelihoods.append(math.exp(site_log_likelihood))
        mixture_likelihood = sum(
            category.weight * category_likelihood
            for category, category_likelihood in zip(
                categories,
                category_likelihoods,
                strict=True,
            )
        )
        pattern_log_likelihood = math.log(mixture_likelihood)
        for site_position in pattern.site_positions:
            site_likelihoods.append(
                DiscreteGammaSiteLikelihood(
                    pattern_id=pattern.pattern_id,
                    site_position=site_position,
                    category_likelihoods=list(category_likelihoods),
                    mixture_likelihood=mixture_likelihood,
                    log_likelihood=pattern_log_likelihood,
                )
            )
        total_log_likelihood += pattern.weight * pattern_log_likelihood
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

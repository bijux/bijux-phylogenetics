from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    dna_observation_state_order,
    normalize_dna_likelihood_records,
    resolve_dna_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CompiledFiniteStateLikelihoodAnalysis,
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_finite_state_tree_likelihood_with_invariant_mixture,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FiniteStateDiscreteGammaInvariantTreeLikelihoodReport,
    FiniteStateInvariantMixtureTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_models import (
    resolve_selected_nucleotide_likelihood_specification,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def build_selected_nucleotide_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    owner_name: str,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state DNA likelihood surface."""
    specification = resolve_selected_nucleotide_likelihood_specification(
        records,
        model_name=model_name,
        owner_name=owner_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    transition_evaluator = build_transition_matrix_evaluator(specification.rate_matrix)
    return FiniteStateLikelihoodSurface(
        model_name=specification.model_name,
        state_labels=tuple(
            dna_observation_state_order(
                observation_policy=specification.observation_policy
            )
        ),
        root_prior=specification.root_prior,
        root_prior_source=specification.root_prior_source,
        observation_policy=specification.observation_policy,
        normalization_policy=specification.normalization_policy,
        transition_rate_matrix=specification.rate_matrix,
        transition_matrix_for_branch_length=transition_evaluator.transition_probability_matrix,
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_dna_observation_leaf_vector(
                {taxon_name: observed_state},
                model_name=specification.model_name,
                node_name=taxon_name,
                observation_policy=specification.observation_policy,
            )
        ),
        parameter_values=specification.parameter_values,
    )


def compile_selected_nucleotide_likelihood_analysis(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    owner_name: str,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one reusable nucleotide finite-state likelihood analysis."""
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name=owner_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name=surface.model_name,
        observation_policy=surface.observation_policy,
    )
    return compile_finite_state_likelihood_analysis(
        tree,
        surface,
        records=normalized_records,
        pattern_compiler=lambda partition_records, partition_source_path: (
            compress_alignment_site_patterns_from_records(
                partition_records,
                source_path=partition_source_path,
            )
        ),
        missing_data_policy=surface.observation_policy,
        partition_map=partition_map,
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def evaluate_selected_nucleotide_tree_likelihood_with_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    owner_name: str,
    invariant_proportion: float,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FiniteStateInvariantMixtureTreeLikelihoodReport:
    """Evaluate one selected nucleotide likelihood surface under one active +I mixture."""
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name=owner_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name=surface.model_name,
        observation_policy=surface.observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    return evaluate_finite_state_tree_likelihood_with_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        invariant_proportion=invariant_proportion,
    )


def evaluate_selected_nucleotide_tree_likelihood_with_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    owner_name: str,
    invariant_proportion: float,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FiniteStateInvariantMixtureTreeLikelihoodReport:
    """Evaluate one selected nucleotide +I likelihood from tree and alignment paths."""
    return evaluate_selected_nucleotide_tree_likelihood_with_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        owner_name=owner_name,
        invariant_proportion=invariant_proportion,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    model_name: str,
    owner_name: str,
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FiniteStateDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one selected nucleotide likelihood surface under one active +G+I mixture."""
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name=owner_name,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name=surface.model_name,
        observation_policy=surface.observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    return evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=alpha,
        category_count=category_count,
        invariant_proportion=invariant_proportion,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        emit_boundary_warnings=emit_boundary_warnings,
    )


def evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    model_name: str,
    owner_name: str,
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
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
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> FiniteStateDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one selected nucleotide +G+I likelihood from tree and alignment paths."""
    return evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        model_name=model_name,
        owner_name=owner_name,
        alpha=alpha,
        invariant_proportion=invariant_proportion,
        category_count=category_count,
        lower_invariant_proportion_bound=lower_invariant_proportion_bound,
        upper_invariant_proportion_bound=upper_invariant_proportion_bound,
        emit_boundary_warnings=emit_boundary_warnings,
        observation_policy=observation_policy,
        kappa=kappa,
        base_frequencies=base_frequencies,
        exchangeabilities=exchangeabilities,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )

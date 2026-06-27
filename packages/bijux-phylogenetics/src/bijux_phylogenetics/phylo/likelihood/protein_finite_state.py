from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CompiledFiniteStateLikelihoodAnalysis,
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
)
from bijux_phylogenetics.phylo.likelihood.poisson import (
    protein_poisson_rate_matrix,
    protein_poisson_transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    PROTEIN_STATE_ORDER,
    UNIFORM_PROTEIN_ROOT_PRIOR,
    normalize_protein_likelihood_records,
    protein_leaf_likelihood_vector,
    resolve_protein_root_prior,
    validate_empirical_protein_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.protein_empirical_models import (
    empirical_protein_rate_matrix,
    validate_empirical_protein_model_name,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


def build_protein_poisson_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state protein Poisson likelihood surface."""
    normalize_protein_likelihood_records(
        records,
        model_name="protein Poisson",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    rate_matrix = protein_poisson_rate_matrix()
    resolved_root_prior = resolve_protein_root_prior(
        records,
        owner_name="protein Poisson finite-state surface",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=UNIFORM_PROTEIN_ROOT_PRIOR,
    )
    return FiniteStateLikelihoodSurface(
        model_name="protein Poisson",
        state_labels=PROTEIN_STATE_ORDER,
        root_prior=resolved_root_prior.root_prior,
        root_prior_source=resolved_root_prior.root_prior_source,
        observation_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=protein_poisson_transition_probability_matrix,
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            protein_leaf_likelihood_vector(
                {taxon_name: observed_state},
                model_name="protein Poisson",
                node_name=taxon_name,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
                ambiguity_policy=ambiguity_policy,
            )
        ),
    )


def build_empirical_protein_finite_state_surface(
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
    normalization_policy: str = CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state empirical protein likelihood surface."""
    normalize_protein_likelihood_records(
        records,
        model_name="empirical protein matrix",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    validated_rate_matrix = validate_empirical_protein_rate_matrix(
        rate_matrix,
        model_name="empirical protein matrix",
    )
    resolved_root_prior = resolve_protein_root_prior(
        records,
        owner_name="empirical protein finite-state surface",
        default_policy="equal",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=_stationary_distribution_for_empirical_matrix(
            validated_rate_matrix
        ),
    )
    return FiniteStateLikelihoodSurface(
        model_name=f"empirical protein matrix ({matrix_label})",
        state_labels=PROTEIN_STATE_ORDER,
        root_prior=resolved_root_prior.root_prior,
        root_prior_source=resolved_root_prior.root_prior_source,
        observation_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
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
    )


def build_named_empirical_protein_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    model_name: str,
    root_prior_policy: str | None = None,
    root_prior: numpy.ndarray | list[float] | tuple[float, ...] | None = None,
    fixed_root_state: str | None = None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state surface from one named empirical protein model."""
    validated_model_name = validate_empirical_protein_model_name(model_name)
    return build_empirical_protein_finite_state_surface(
        records,
        rate_matrix=empirical_protein_rate_matrix(validated_model_name),
        root_prior_policy=(
            "stationary" if root_prior_policy is None else root_prior_policy
        ),
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=validated_model_name,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    )


def compile_protein_poisson_likelihood_analysis(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one reusable protein Poisson finite-state likelihood analysis."""
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="protein Poisson",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    surface = build_protein_poisson_finite_state_surface(
        normalized_records,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
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
        missing_data_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        partition_map=partition_map,
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def compile_empirical_protein_likelihood_analysis(
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
    normalization_policy: str = CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one reusable empirical protein finite-state likelihood analysis."""
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="empirical protein matrix",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    surface = build_empirical_protein_finite_state_surface(
        normalized_records,
        rate_matrix=rate_matrix,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=matrix_label,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
        normalization_policy=normalization_policy,
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
        missing_data_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        partition_map=partition_map,
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def compile_named_empirical_protein_likelihood_analysis(
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
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one reusable named empirical protein finite-state likelihood analysis."""
    validated_model_name = validate_empirical_protein_model_name(model_name)
    return compile_empirical_protein_likelihood_analysis(
        tree,
        records,
        rate_matrix=empirical_protein_rate_matrix(validated_model_name),
        root_prior_policy=(
            "stationary" if root_prior_policy is None else root_prior_policy
        ),
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        matrix_label=validated_model_name,
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        partition_map=partition_map,
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def _protein_observation_policy_label(
    *,
    gap_policy: str,
    missing_policy: str,
    ambiguity_policy: str,
) -> str:
    return f"gap={gap_policy};missing={missing_policy};ambiguity={ambiguity_policy}"


def _stationary_distribution_for_empirical_matrix(
    rate_matrix: numpy.ndarray,
) -> numpy.ndarray:
    from bijux_phylogenetics.phylo.likelihood.ctmc import (
        solve_ctmc_stationary_distribution,
    )

    report = solve_ctmc_stationary_distribution(
        rate_matrix,
        state_labels=PROTEIN_STATE_ORDER,
    )
    return numpy.asarray(report.probabilities, dtype=float)

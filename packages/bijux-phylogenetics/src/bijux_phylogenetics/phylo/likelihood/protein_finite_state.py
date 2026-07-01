from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
from pathlib import Path

import numpy

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CompiledFiniteStateLikelihoodAnalysis,
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FiniteStateFrequencyOptimizationReport,
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
from bijux_phylogenetics.phylo.likelihood.state_frequency_optimization import (
    build_equal_input_rate_matrix,
    optimize_finite_state_stationary_frequencies,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


@dataclass(frozen=True, slots=True)
class _ProteinEqualInputLikelihoodEvaluation:
    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    state_count: int
    log_likelihood: float


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


def build_protein_equal_input_finite_state_surface(
    records: list[AlignmentRecord],
    *,
    stationary_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...],
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
    """Build one shared equal-input protein finite-state likelihood surface."""
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="protein equal-input",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    stationary_vector = _resolve_protein_stationary_vector(stationary_frequencies)
    rate_matrix = build_equal_input_rate_matrix(
        state_labels=PROTEIN_STATE_ORDER,
        stationary_frequencies=stationary_vector,
        owner_name="protein equal-input finite-state surface",
    )
    resolved_root_prior = resolve_protein_root_prior(
        normalized_records,
        owner_name="protein equal-input finite-state surface",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=stationary_vector,
    )
    return FiniteStateLikelihoodSurface(
        model_name="protein equal-input",
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
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            protein_leaf_likelihood_vector(
                {taxon_name: observed_state},
                model_name="protein equal-input",
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


def optimize_protein_equal_input_state_frequencies(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    initial_state_frequencies: dict[str, float] | None = None,
    gap_policy: str = "treat-as-missing",
    missing_policy: str = "treat-as-missing",
    ambiguity_policy: str = "reject",
    minimum_state_frequency: float = 1e-6,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> FiniteStateFrequencyOptimizationReport:
    """Optimize one fixed-topology protein equal-input stationary simplex."""
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="protein equal-input",
        gap_policy=gap_policy,
        missing_policy=missing_policy,
        ambiguity_policy=ambiguity_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if initial_state_frequencies is None:
        initial_values_by_name = _estimate_protein_state_frequencies(normalized_records)
    else:
        initial_values_by_name = dict(initial_state_frequencies)
    return optimize_finite_state_stationary_frequencies(
        model_name="protein equal-input",
        state_labels=PROTEIN_STATE_ORDER,
        initial_values_by_name=initial_values_by_name,
        evaluate=lambda candidate: _evaluate_protein_equal_input_frequency_candidate(
            tree,
            compressed_patterns,
            candidate,
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        minimum_state_frequency=minimum_state_frequency,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
        nonconvergence_warning="protein equal-input stationary-frequency search did not converge",
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


def _resolve_protein_stationary_vector(
    stationary_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...],
) -> numpy.ndarray:
    if isinstance(stationary_frequencies, dict):
        return numpy.asarray(
            [stationary_frequencies[state] for state in PROTEIN_STATE_ORDER],
            dtype=float,
        )
    return numpy.asarray(stationary_frequencies, dtype=float)


def _estimate_protein_state_frequencies(
    records: list[AlignmentRecord],
) -> dict[str, float]:
    amino_acid_counts = dict.fromkeys(PROTEIN_STATE_ORDER, 0.0)
    for record in records:
        for state in record.sequence:
            if state in amino_acid_counts:
                amino_acid_counts[state] += 1.0
    total_count = sum(amino_acid_counts.values())
    if total_count <= 0.0:
        uniform_mass = 1.0 / len(PROTEIN_STATE_ORDER)
        return dict.fromkeys(PROTEIN_STATE_ORDER, uniform_mass)
    return {
        state: count / total_count for state, count in amino_acid_counts.items()
    }


def _evaluate_protein_equal_input_frequency_candidate(
    tree: PhyloTree,
    compressed_patterns,
    candidate_state_frequencies: dict[str, float],
    *,
    gap_policy: str,
    missing_policy: str,
    ambiguity_policy: str,
) -> tuple[_ProteinEqualInputLikelihoodEvaluation, float]:
    state_order = PROTEIN_STATE_ORDER
    rate_matrix = build_equal_input_rate_matrix(
        state_labels=state_order,
        stationary_frequencies=candidate_state_frequencies,
        owner_name="protein equal-input",
    )
    stationary_vector = numpy.asarray(
        [candidate_state_frequencies[state] for state in state_order],
        dtype=float,
    )
    surface = FiniteStateLikelihoodSurface(
        model_name="protein equal-input",
        state_labels=state_order,
        root_prior=stationary_vector,
        root_prior_source="stationary",
        observation_policy=_protein_observation_policy_label(
            gap_policy=gap_policy,
            missing_policy=missing_policy,
            ambiguity_policy=ambiguity_policy,
        ),
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            protein_leaf_likelihood_vector(
                {taxon_name: observed_state},
                model_name="protein equal-input",
                node_name=taxon_name,
                gap_policy=gap_policy,
                missing_policy=missing_policy,
                ambiguity_policy=ambiguity_policy,
            )
        ),
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    payload = _ProteinEqualInputLikelihoodEvaluation(
        taxa=list(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        tree_newick=dumps_newick(tree),
        state_count=len(state_order),
        log_likelihood=log_likelihood,
    )
    return payload, log_likelihood


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

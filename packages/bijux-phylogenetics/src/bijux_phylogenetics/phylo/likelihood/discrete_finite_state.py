from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping
from pathlib import Path

import numpy

from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.phylo.likelihood.discrete_observation_policies import (
    resolve_discrete_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
    solve_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.finite_state_root_priors import (
    resolve_finite_state_root_prior,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CompiledFiniteStateLikelihoodAnalysis,
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_finite_state_tree_likelihood_with_invariant_mixture,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    FiniteStateFrequencyOptimizationReport,
    FiniteStateDiscreteGammaInvariantTreeLikelihoodReport,
    FiniteStateInvariantMixtureTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.state_frequency_optimization import (
    build_equal_input_rate_matrix,
    optimize_finite_state_stationary_frequencies,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


@dataclass(frozen=True, slots=True)
class _DiscreteEqualInputLikelihoodEvaluation:
    taxa: list[str]
    site_count: int
    pattern_count: int
    tree_newick: str
    state_count: int
    log_likelihood: float


def build_discrete_mk_finite_state_surface(
    *,
    state_order: list[str] | tuple[str, ...],
    rate_matrix: numpy.ndarray,
    records: list[AlignmentRecord] | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    observation_policy: str = "reject",
    model_name: str = "discrete Mk",
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state surface for one discrete Mk rate matrix."""
    validated_state_order = tuple(state_order)
    validated_rate_matrix = numpy.asarray(rate_matrix, dtype=float)
    resolved_root_prior = resolve_finite_state_root_prior(
        records or [],
        state_labels=validated_state_order,
        owner_name=f"{model_name} finite-state surface",
        default_policy="equal",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=_stationary_distribution_for_rate_matrix(
            validated_rate_matrix,
            state_labels=validated_state_order,
        ),
        sequence_tokenizer=_discrete_site_tokenizer,
    )
    state_index = {
        state: index for index, state in enumerate(validated_state_order)
    }
    return FiniteStateLikelihoodSurface(
        model_name=model_name,
        state_labels=validated_state_order,
        root_prior=resolved_root_prior.root_prior,
        root_prior_source=resolved_root_prior.root_prior_source,
        observation_policy=observation_policy,
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
        transition_rate_matrix=validated_rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(validated_rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_discrete_observation_leaf_vector(
                observed_state,
                model_name=model_name,
                node_name=taxon_name,
                state_order=list(validated_state_order),
                observation_policy=observation_policy,
            )
        ),
        parameter_values={"state_count": float(len(state_index))},
    )


def compile_discrete_mk_likelihood_analysis(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    state_order: list[str] | tuple[str, ...],
    rate_matrix: numpy.ndarray,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    observation_policy: str = "reject",
    model_name: str = "discrete Mk",
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one reusable discrete Mk finite-state likelihood analysis."""
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        records=records,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        observation_policy=observation_policy,
        model_name=model_name,
    )
    return compile_finite_state_likelihood_analysis(
        tree,
        surface,
        records=records,
        pattern_compiler=_compile_discrete_patterns,
        missing_data_policy=observation_policy,
        alignment_length_resolver=_resolve_discrete_alignment_length,
        partition_map=partition_map,
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def evaluate_discrete_mk_tree_likelihood_with_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    state_order: list[str] | tuple[str, ...],
    rate_matrix: numpy.ndarray,
    invariant_proportion: float,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    observation_policy: str = "reject",
    model_name: str = "discrete Mk",
) -> FiniteStateInvariantMixtureTreeLikelihoodReport:
    """Evaluate one discrete Mk likelihood surface under one active +I mixture."""
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        records=records,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        observation_policy=observation_policy,
        model_name=model_name,
    )
    compressed_patterns = _compile_discrete_patterns(records, source_path=None)
    return evaluate_finite_state_tree_likelihood_with_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        invariant_proportion=invariant_proportion,
    )


def evaluate_discrete_mk_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    state_order: list[str] | tuple[str, ...],
    rate_matrix: numpy.ndarray,
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    observation_policy: str = "reject",
    model_name: str = "discrete Mk",
) -> FiniteStateDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one discrete Mk likelihood surface under one active +G+I mixture."""
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        records=records,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        observation_policy=observation_policy,
        model_name=model_name,
    )
    compressed_patterns = _compile_discrete_patterns(records, source_path=None)
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


def optimize_discrete_equal_input_state_frequencies(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    state_order: list[str] | tuple[str, ...],
    initial_state_frequencies: dict[str, float] | None = None,
    observation_policy: str = "reject",
    model_name: str = "discrete equal-input",
    minimum_state_frequency: float = 1e-6,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> FiniteStateFrequencyOptimizationReport:
    """Optimize one fixed-topology equal-input stationary simplex for discrete states."""
    validated_state_order = tuple(state_order)
    compressed_patterns = _compile_discrete_patterns(records, source_path=None)
    if initial_state_frequencies is None:
        initial_values_by_name = _estimate_discrete_state_frequencies(
            records,
            state_order=validated_state_order,
        )
    else:
        initial_values_by_name = dict(initial_state_frequencies)
    return optimize_finite_state_stationary_frequencies(
        model_name=model_name,
        state_labels=validated_state_order,
        initial_values_by_name=initial_values_by_name,
        evaluate=lambda candidate: _evaluate_discrete_equal_input_frequency_candidate(
            tree,
            compressed_patterns,
            state_order=validated_state_order,
            candidate_state_frequencies=candidate,
            observation_policy=observation_policy,
            model_name=model_name,
        ),
        minimum_state_frequency=minimum_state_frequency,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
        nonconvergence_warning=(
            f"{model_name} stationary-frequency search did not converge"
        ),
    )


def _compile_discrete_patterns(
    records: list[AlignmentRecord],
    source_path: Path | None,
) -> CompressedAlignmentSitePatterns:
    return compress_alignment_site_patterns_from_records(
        records,
        source_path=source_path,
        site_tokenizer=_discrete_site_tokenizer,
    )


def _resolve_discrete_alignment_length(_records: list[AlignmentRecord]) -> int:
    return 1


def _discrete_site_tokenizer(sequence: str) -> tuple[str, ...]:
    return (sequence,)


def _estimate_discrete_state_frequencies(
    records: list[AlignmentRecord],
    *,
    state_order: tuple[str, ...],
) -> dict[str, float]:
    state_counts = dict.fromkeys(state_order, 0.0)
    if not records:
        uniform_mass = 1.0 / len(state_order)
        return dict.fromkeys(state_order, uniform_mass)
    for record in records:
        if record.sequence in state_counts:
            state_counts[record.sequence] += 1.0
    total_count = sum(state_counts.values())
    if total_count <= 0.0:
        uniform_mass = 1.0 / len(state_order)
        return dict.fromkeys(state_order, uniform_mass)
    return {
        state: count / total_count for state, count in state_counts.items()
    }


def _evaluate_discrete_equal_input_frequency_candidate(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    state_order: tuple[str, ...],
    candidate_state_frequencies: dict[str, float],
    observation_policy: str,
    model_name: str,
) -> tuple[_DiscreteEqualInputLikelihoodEvaluation, float]:
    rate_matrix = build_equal_input_rate_matrix(
        state_labels=state_order,
        stationary_frequencies=candidate_state_frequencies,
        owner_name=model_name,
    )
    stationary_vector = numpy.asarray(
        [candidate_state_frequencies[state] for state in state_order],
        dtype=float,
    )
    surface = FiniteStateLikelihoodSurface(
        model_name=model_name,
        state_labels=state_order,
        root_prior=stationary_vector,
        root_prior_source="stationary",
        observation_policy=observation_policy,
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_discrete_observation_leaf_vector(
                observed_state,
                model_name=model_name,
                node_name=taxon_name,
                state_order=list(state_order),
                observation_policy=observation_policy,
            )
        ),
    )
    context = compile_finite_state_likelihood_context(
        tree,
        surface,
    )
    log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    payload = _DiscreteEqualInputLikelihoodEvaluation(
        taxa=list(compressed_patterns.taxon_order),
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        tree_newick=dumps_newick(tree),
        state_count=len(state_order),
        log_likelihood=log_likelihood,
    )
    return payload, log_likelihood


def _stationary_distribution_for_rate_matrix(
    rate_matrix: numpy.ndarray,
    *,
    state_labels: tuple[str, ...],
) -> numpy.ndarray:
    report = solve_ctmc_stationary_distribution(
        rate_matrix,
        state_labels=state_labels,
    )
    return numpy.asarray(report.probabilities, dtype=float)

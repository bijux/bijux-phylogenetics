from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.codon_observation_policies import (
    resolve_codon_observation_leaf_vector,
    validate_codon_observation_policy,
    validate_codon_observation_state,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    CodonStateSpace,
    build_equal_rate_codon_ctmc_rate_matrix,
    resolve_codon_state_space,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_context,
    evaluate_finite_state_site_log_likelihood_rows_from_patterns,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)
from bijux_phylogenetics.phylo.likelihood.finite_state_root_priors import (
    ResolvedFiniteStateRootPrior,
    resolve_finite_state_root_prior,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    CodonCtmcTreeLikelihoodReport,
    FiniteStateDiscreteGammaInvariantTreeLikelihoodReport,
    FiniteStateFrequencyOptimizationReport,
    FiniteStateInvariantMixtureTreeLikelihoodReport,
    SiteLogLikelihoodRow,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    validate_site_log_likelihood_reconstruction,
)
from bijux_phylogenetics.phylo.likelihood.state_frequency_optimization import (
    optimize_finite_state_stationary_frequencies,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import (
    AlignmentTaxonMismatchError,
    InvalidAlignmentError,
)


@dataclass(frozen=True, slots=True)
class ResolvedCodonCtmcModel:
    """One normalized codon CTMC likelihood surface definition."""

    normalized_records: list[AlignmentRecord]
    state_space: CodonStateSpace
    rate_matrix: numpy.ndarray
    stationary_frequencies: numpy.ndarray
    codon_frequency_source: str
    observation_policy: str


def evaluate_codon_ctmc_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> CodonCtmcTreeLikelihoodReport:
    """Evaluate one fixed-topology codon CTMC likelihood on an aligned coding matrix."""
    resolved_model = resolve_codon_ctmc_model(
        records,
        codon_frequencies=codon_frequencies,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        resolved_model.normalized_records,
        site_tokenizer=iter_codon_tokens,
    )
    resolved_root_prior = resolve_codon_root_prior(
        resolved_model.normalized_records,
        state_space=resolved_model.state_space,
        owner_name="codon CTMC likelihood",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=resolved_model.stationary_frequencies,
    )
    log_likelihood = _evaluate_codon_ctmc_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        state_space=resolved_model.state_space,
        rate_matrix=resolved_model.rate_matrix,
        root_prior=resolved_root_prior.root_prior,
        observation_policy=resolved_model.observation_policy,
    )
    site_log_likelihoods, reconstructed_log_likelihood = (
        _evaluate_codon_ctmc_tree_likelihood_rows_from_patterns(
            tree,
            compressed_patterns,
            state_space=resolved_model.state_space,
            rate_matrix=resolved_model.rate_matrix,
            root_prior=resolved_root_prior.root_prior,
            observation_policy=resolved_model.observation_policy,
        )
    )
    validate_site_log_likelihood_reconstruction(
        site_log_likelihoods,
        expected_total_log_likelihood=reconstructed_log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name="codon CTMC likelihood",
    )
    return CodonCtmcTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=(
            compressed_patterns.pattern_count < compressed_patterns.alignment_length
        ),
        expansion_policy="expanded-site-rows",
        tree_newick=dumps_newick(tree),
        state_count=len(resolved_model.state_space.state_order),
        genetic_code_id=resolved_model.state_space.genetic_code_id,
        genetic_code_name=resolved_model.state_space.genetic_code_name,
        codon_frequency_source=resolved_model.codon_frequency_source,
        root_prior_source=resolved_root_prior.root_prior_source,
        observation_policy=resolved_model.observation_policy,
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        log_likelihood=log_likelihood,
        site_log_likelihoods=site_log_likelihoods,
    )


def evaluate_codon_ctmc_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> CodonCtmcTreeLikelihoodReport:
    """Evaluate one fixed-topology codon CTMC likelihood from tree and alignment paths."""
    return evaluate_codon_ctmc_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        codon_frequencies=codon_frequencies,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )


def optimize_codon_ctmc_state_frequencies(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    initial_codon_frequencies: dict[str, float] | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
    minimum_state_frequency: float = 1e-6,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> FiniteStateFrequencyOptimizationReport:
    """Optimize one fixed-topology codon CTMC stationary simplex."""
    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    state_space = resolve_codon_state_space(genetic_code)
    if initial_codon_frequencies is None:
        initial_values_by_name = _estimate_codon_state_frequencies(
            normalized_records,
            state_space=state_space,
        )
    else:
        initial_values_by_name = dict(initial_codon_frequencies)
    return optimize_finite_state_stationary_frequencies(
        model_name="codon CTMC",
        state_labels=state_space.state_order,
        initial_values_by_name=initial_values_by_name,
        evaluate=lambda candidate: _evaluate_codon_frequency_candidate(
            tree,
            normalized_records,
            candidate,
            genetic_code=genetic_code,
            observation_policy=observation_policy,
        ),
        minimum_state_frequency=minimum_state_frequency,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
        nonconvergence_warning="codon CTMC stationary-frequency search did not converge",
    )


def optimize_codon_ctmc_state_frequencies_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    initial_codon_frequencies: dict[str, float] | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
    minimum_state_frequency: float = 1e-6,
    improvement_tolerance: float = 1e-9,
    max_coordinate_passes: int = 12,
) -> FiniteStateFrequencyOptimizationReport:
    """Optimize one codon CTMC stationary simplex from tree and alignment paths."""
    return optimize_codon_ctmc_state_frequencies(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        initial_codon_frequencies=initial_codon_frequencies,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
        minimum_state_frequency=minimum_state_frequency,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )


def evaluate_codon_ctmc_tree_likelihood_with_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    invariant_proportion: float,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> FiniteStateInvariantMixtureTreeLikelihoodReport:
    """Evaluate one codon CTMC likelihood surface under one active +I mixture."""
    from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
        build_codon_finite_state_surface,
    )
    from bijux_phylogenetics.phylo.likelihood.finite_state import (
        evaluate_finite_state_tree_likelihood_with_invariant_mixture,
    )

    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_codon_site_patterns_from_records(normalized_records)
    surface = build_codon_finite_state_surface(
        normalized_records,
        codon_frequencies=codon_frequencies,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    return evaluate_finite_state_tree_likelihood_with_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        invariant_proportion=invariant_proportion,
    )


def evaluate_codon_ctmc_tree_likelihood_with_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    invariant_proportion: float,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> FiniteStateInvariantMixtureTreeLikelihoodReport:
    """Evaluate one codon CTMC +I likelihood from tree and alignment paths."""
    return evaluate_codon_ctmc_tree_likelihood_with_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        invariant_proportion=invariant_proportion,
        codon_frequencies=codon_frequencies,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )


def evaluate_codon_ctmc_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> FiniteStateDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one codon CTMC likelihood surface under one active +G+I mixture."""
    from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
        build_codon_finite_state_surface,
    )

    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_codon_site_patterns_from_records(normalized_records)
    surface = build_codon_finite_state_surface(
        normalized_records,
        codon_frequencies=codon_frequencies,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    return (
        evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
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
    )


def evaluate_codon_ctmc_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    alpha: float,
    invariant_proportion: float,
    category_count: int = 4,
    lower_invariant_proportion_bound: float | None = None,
    upper_invariant_proportion_bound: float | None = None,
    emit_boundary_warnings: bool = False,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> FiniteStateDiscreteGammaInvariantTreeLikelihoodReport:
    """Evaluate one codon CTMC +G+I likelihood from tree and alignment paths."""
    return (
        evaluate_codon_ctmc_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            load_tree(tree_path),
            load_fasta_alignment(alignment_path),
            alpha=alpha,
            invariant_proportion=invariant_proportion,
            category_count=category_count,
            lower_invariant_proportion_bound=lower_invariant_proportion_bound,
            upper_invariant_proportion_bound=upper_invariant_proportion_bound,
            emit_boundary_warnings=emit_boundary_warnings,
            codon_frequencies=codon_frequencies,
            root_prior_policy=root_prior_policy,
            root_prior=root_prior,
            fixed_root_state=fixed_root_state,
            genetic_code=genetic_code,
            observation_policy=observation_policy,
        )
    )


def compress_codon_site_patterns_from_records(
    records: list[AlignmentRecord],
    *,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
    source_path: Path | None = None,
) -> CompressedAlignmentSitePatterns:
    """Compress identical aligned codon sites through the shared pattern surface."""
    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    return compress_alignment_site_patterns_from_records(
        normalized_records,
        source_path=source_path,
        site_tokenizer=_iter_record_codons,
    )


def normalize_codon_likelihood_records(
    records: list[AlignmentRecord],
    *,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> list[AlignmentRecord]:
    """Normalize one aligned coding matrix into unambiguous sense codons only."""
    state_space = resolve_codon_state_space(genetic_code)
    if not records:
        raise InvalidAlignmentError(
            "codon CTMC likelihood requires at least one aligned coding sequence"
        )
    alignment_lengths = {len(record.sequence) for record in records}
    if len(alignment_lengths) != 1:
        raise InvalidAlignmentError(
            "codon CTMC likelihood requires equal-length aligned coding sequences"
        )
    alignment_length = next(iter(alignment_lengths))
    if alignment_length == 0 or alignment_length % 3 != 0:
        raise InvalidAlignmentError(
            "codon CTMC likelihood requires an aligned coding matrix whose length is divisible by three"
        )
    normalized_records: list[AlignmentRecord] = []
    for record in records:
        normalized_sequence = record.sequence.upper().replace("U", "T")
        for codon_index, codon in enumerate(
            _iter_record_codons(normalized_sequence),
            start=1,
        ):
            validate_codon_observation_state(
                codon,
                codon_index=codon_index,
                genetic_code_name=state_space.genetic_code_name,
                owner_name="codon CTMC likelihood",
                observation_policy=observation_policy,
                record_identifier=record.identifier,
                stop_codons=frozenset(state_space.stop_codons),
            )
        normalized_records.append(
            AlignmentRecord(
                identifier=record.identifier,
                sequence=normalized_sequence,
            )
        )
    return normalized_records


def resolve_codon_ctmc_model(
    records: list[AlignmentRecord],
    *,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> ResolvedCodonCtmcModel:
    """Resolve one normalized codon CTMC model definition for owned surfaces."""
    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    state_space, rate_matrix, frequencies, frequency_source = (
        build_equal_rate_codon_ctmc_rate_matrix(
            codon_frequencies,
            genetic_code=genetic_code,
        )
    )
    validated_observation_policy = validate_codon_observation_policy(
        observation_policy,
        owner_name="codon CTMC likelihood",
    )
    return ResolvedCodonCtmcModel(
        normalized_records=normalized_records,
        state_space=state_space,
        rate_matrix=rate_matrix,
        stationary_frequencies=frequencies,
        codon_frequency_source=frequency_source,
        observation_policy=validated_observation_policy,
    )


def resolve_codon_root_prior(
    records: list[AlignmentRecord],
    *,
    state_space: CodonStateSpace,
    owner_name: str,
    default_policy: str,
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
    stationary_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
) -> ResolvedFiniteStateRootPrior:
    return resolve_finite_state_root_prior(
        records,
        state_labels=state_space.state_order,
        owner_name=owner_name,
        default_policy=default_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=stationary_frequencies,
        sequence_tokenizer=iter_codon_tokens,
    )


def _estimate_codon_state_frequencies(
    records: list[AlignmentRecord],
    *,
    state_space: CodonStateSpace,
) -> dict[str, float]:
    codon_counts = dict.fromkeys(state_space.state_order, 0.0)
    for record in records:
        for codon in iter_codon_tokens(record.sequence):
            if codon in codon_counts:
                codon_counts[codon] += 1.0
    total_count = sum(codon_counts.values())
    if total_count <= 0.0:
        uniform_mass = 1.0 / len(state_space.state_order)
        return dict.fromkeys(state_space.state_order, uniform_mass)
    return {codon: count / total_count for codon, count in codon_counts.items()}


def _evaluate_codon_frequency_candidate(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    candidate_codon_frequencies: dict[str, float],
    *,
    genetic_code: int | str | None,
    observation_policy: str,
) -> tuple[CodonCtmcTreeLikelihoodReport, float]:
    report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        codon_frequencies=candidate_codon_frequencies,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    return report, report.log_likelihood


def _evaluate_codon_ctmc_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    state_space: CodonStateSpace,
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
    observation_policy: str,
) -> float:
    validate_explicit_branch_lengths(tree, model_name="codon CTMC")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="codon CTMC",
    )
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="codon CTMC",
            state_labels=tuple(state_space.state_order),
            root_prior=root_prior,
            root_prior_source="provided",
            observation_policy=observation_policy,
            normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
            transition_rate_matrix=rate_matrix,
            transition_matrix_for_branch_length=(
                transition_evaluator.transition_probability_matrix
            ),
            leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
                resolve_codon_observation_leaf_vector(
                    observed_state,
                    model_name="codon CTMC",
                    node_name=taxon_name,
                    observation_policy=observation_policy,
                    state_order=state_space.state_order,
                    state_index=state_space.state_index,
                )
            ),
        ),
    )
    log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    site_log_likelihoods, reconstructed_log_likelihood = (
        evaluate_finite_state_site_log_likelihood_rows_from_patterns(
            context,
            compressed_patterns,
        )
    )
    validate_site_log_likelihood_reconstruction(
        site_log_likelihoods,
        expected_total_log_likelihood=reconstructed_log_likelihood,
        expected_site_count=compressed_patterns.alignment_length,
        expected_pattern_count=compressed_patterns.pattern_count,
        owner_name="codon CTMC likelihood",
    )
    return log_likelihood


def _evaluate_codon_ctmc_tree_likelihood_rows_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    state_space: CodonStateSpace,
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
    observation_policy: str,
) -> tuple[list[SiteLogLikelihoodRow], float]:
    validate_explicit_branch_lengths(tree, model_name="codon CTMC")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="codon CTMC",
    )
    transition_evaluator = build_transition_matrix_evaluator(rate_matrix)
    context = compile_finite_state_likelihood_context(
        tree,
        FiniteStateLikelihoodSurface(
            model_name="codon CTMC",
            state_labels=tuple(state_space.state_order),
            root_prior=root_prior,
            root_prior_source="provided",
            observation_policy=observation_policy,
            normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
            transition_rate_matrix=rate_matrix,
            transition_matrix_for_branch_length=(
                transition_evaluator.transition_probability_matrix
            ),
            leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
                resolve_codon_observation_leaf_vector(
                    observed_state,
                    model_name="codon CTMC",
                    node_name=taxon_name,
                    observation_policy=observation_policy,
                    state_order=state_space.state_order,
                    state_index=state_space.state_index,
                )
            ),
        ),
    )
    return evaluate_finite_state_site_log_likelihood_rows_from_patterns(
        context,
        compressed_patterns,
    )


def _evaluate_codon_site_log_likelihood(
    tree: PhyloTree,
    states: tuple[str, ...],
    *,
    taxon_order: list[str],
    state_space: CodonStateSpace,
    root_prior: numpy.ndarray,
    transition_evaluator,
    observation_policy: str,
) -> float:
    states_by_taxon = dict(zip(taxon_order, states, strict=True))
    pruning_pass = postorder_conditional_likelihoods(
        tree,
        state_count=len(state_space.state_order),
        leaf_likelihood=lambda node: _codon_leaf_likelihood_vector(
            states_by_taxon,
            state_space=state_space,
            node_name=node.name,
            observation_policy=observation_policy,
        ),
        transition_matrix_for_child=lambda child: (
            transition_evaluator.transition_probability_matrix(
                float(child.branch_length or 0.0)
            )
        ),
    )
    return log_likelihood_from_root_prior(
        tree,
        pruning_pass,
        root_prior=root_prior,
    )


def _codon_leaf_likelihood_vector(
    states_by_taxon: dict[str, str],
    *,
    state_space: CodonStateSpace,
    node_name: str | None,
    observation_policy: str,
) -> numpy.ndarray:
    if node_name is None:
        raise AlignmentTaxonMismatchError(
            "codon CTMC likelihood requires named tree tips for alignment lookup"
        )
    codon = states_by_taxon[node_name]
    return resolve_codon_observation_leaf_vector(
        codon,
        model_name="codon CTMC",
        node_name=node_name,
        observation_policy=observation_policy,
        state_order=state_space.state_order,
        state_index=state_space.state_index,
    )


def iter_codon_tokens(sequence: str) -> tuple[str, ...]:
    return tuple(sequence[index : index + 3] for index in range(0, len(sequence), 3))


def _iter_record_codons(sequence: str) -> tuple[str, ...]:
    return iter_codon_tokens(sequence)

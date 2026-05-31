from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_EXCHANGEABILITY_ORDER,
    DNA_STATE_ORDER,
    normalize_dna_exchangeabilities_by_anchor,
    normalize_dna_rate_matrix,
    validate_dna_base_frequencies,
    validate_dna_exchangeabilities,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    augment_dna_rate_matrix_with_gap_state,
    dna_observation_state_order,
    estimate_empirical_dna_base_frequencies_from_records,
    estimate_empirical_gap_state_frequency,
    normalize_dna_likelihood_records,
    resolve_default_dna_root_prior_for_observation_policy,
    resolve_dna_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.dna_simplex_coordinates import (
    resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained,
)
from bijux_phylogenetics.phylo.likelihood.models import (
    GtrExchangeabilityOptimizationReport,
    GtrTreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.parameter_bounds import (
    validate_parameter_within_bounds,
    validate_positive_parameter_bounds,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_coordinate_likelihood_search,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
    transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    sum_compressed_site_pattern_log_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.validation import (
    validate_explicit_branch_lengths,
    validate_tree_taxa_against_patterns,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

_GTR_EXCHANGEABILITY_LABELS = tuple("".join(pair) for pair in DNA_EXCHANGEABILITY_ORDER)
_GTR_EXCHANGEABILITY_ANCHOR = "AC=1"
_GTR_FREE_EXCHANGEABILITY_LABELS = _GTR_EXCHANGEABILITY_LABELS[1:]


def gtr_rate_matrix(
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    *,
    base_frequencies: dict[str, float] | numpy.ndarray,
) -> numpy.ndarray:
    """Return the normalized GTR rate matrix with expected rate one."""
    stationary = validate_dna_base_frequencies(base_frequencies, model_name="GTR")
    validated_exchangeabilities = validate_dna_exchangeabilities(
        exchangeabilities,
        model_name="GTR",
    )
    off_diagonal_rates = numpy.zeros((4, 4), dtype=float)
    for pair_index, (left_state, right_state) in enumerate(DNA_EXCHANGEABILITY_ORDER):
        exchangeability = validated_exchangeabilities[pair_index]
        left_index = DNA_STATE_ORDER.index(left_state)
        right_index = DNA_STATE_ORDER.index(right_state)
        off_diagonal_rates[left_index, right_index] = (
            exchangeability * stationary[right_index]
        )
        off_diagonal_rates[right_index, left_index] = (
            exchangeability * stationary[left_index]
        )
    return normalize_dna_rate_matrix(
        off_diagonal_rates,
        stationary_frequencies=stationary,
        model_name="GTR",
    )


def gtr_transition_probability_matrix(
    branch_length: float,
    *,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    base_frequencies: dict[str, float] | numpy.ndarray,
) -> numpy.ndarray:
    """Return the native GTR transition matrix for one branch."""
    if branch_length <= 0.0:
        return numpy.eye(4, dtype=float)
    return transition_probability_matrix(
        gtr_rate_matrix(
            exchangeabilities,
            base_frequencies=base_frequencies,
        ),
        branch_length,
    )


def evaluate_gtr_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> GtrTreeLikelihoodReport:
    """Evaluate one fixed-topology GTR likelihood from aligned DNA records."""
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="GTR",
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if base_frequencies is None:
        stationary = estimate_empirical_dna_base_frequencies_from_records(
            normalized_records,
            model_name="GTR",
            observation_policy=observation_policy,
        )
        source = "estimated"
    else:
        stationary = validate_dna_base_frequencies(
            base_frequencies,
            model_name="GTR",
        )
        source = "provided"
    gap_state_frequency = (
        estimate_empirical_gap_state_frequency(
            normalized_records,
            model_name="GTR",
        )
        if observation_policy == "fifth-state"
        else None
    )
    resolved_root_prior = resolve_default_dna_root_prior_for_observation_policy(
        normalized_records,
        owner_name="GTR likelihood",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=stationary,
        observation_policy=observation_policy,
    )
    return _evaluate_gtr_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        stationary_frequencies=stationary,
        base_frequency_source=source,
        exchangeabilities=exchangeabilities,
        root_prior=resolved_root_prior.root_prior,
        observation_policy=observation_policy,
        gap_state_frequency=gap_state_frequency,
    )


def evaluate_gtr_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> GtrTreeLikelihoodReport:
    """Evaluate one fixed-topology GTR likelihood from one tree path and alignment."""
    return evaluate_gtr_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        exchangeabilities=exchangeabilities,
        base_frequencies=base_frequencies,
        observation_policy=observation_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    unconstrained_exchangeabilities: Sequence[float],
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> GtrTreeLikelihoodReport:
    """Evaluate one fixed-topology GTR likelihood from unconstrained simplex coordinates."""
    anchored_exchangeabilities = (
        resolve_anchor_normalized_dna_exchangeabilities_from_unconstrained(
            unconstrained_exchangeabilities
        )
    )
    return evaluate_gtr_tree_likelihood(
        tree,
        records,
        exchangeabilities=anchored_exchangeabilities,
        base_frequencies=base_frequencies,
        observation_policy=observation_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    unconstrained_exchangeabilities: Sequence[float],
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> GtrTreeLikelihoodReport:
    """Evaluate one fixed-topology GTR likelihood from unconstrained coordinates and file inputs."""
    return evaluate_gtr_tree_likelihood_from_unconstrained_exchangeabilities(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        unconstrained_exchangeabilities=unconstrained_exchangeabilities,
        base_frequencies=base_frequencies,
        observation_policy=observation_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def optimize_gtr_exchangeabilities(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    initial_exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    lower_exchangeability_bound: float = 0.05,
    upper_exchangeability_bound: float = 20.0,
    max_coordinate_passes: int = 24,
    improvement_tolerance: float = 1e-5,
) -> GtrExchangeabilityOptimizationReport:
    """Optimize one fixed-topology GTR exchangeability surface with AC anchored at one."""
    validated_lower_bound, validated_upper_bound = validate_positive_parameter_bounds(
        parameter_name="exchangeability",
        lower_bound=lower_exchangeability_bound,
        upper_bound=upper_exchangeability_bound,
        owner_name="GTR exchangeability optimization",
    )
    if max_coordinate_passes < 1:
        raise ValueError("GTR exchangeability optimization requires at least one pass")

    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="GTR",
        observation_policy="reject",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if base_frequencies is None:
        stationary = estimate_empirical_dna_base_frequencies_from_records(
            normalized_records,
            model_name="GTR",
            observation_policy="reject",
        )
        source = "estimated"
    else:
        stationary = validate_dna_base_frequencies(
            base_frequencies,
            model_name="GTR",
        )
        source = "provided"
    if initial_exchangeabilities is None:
        normalized_exchangeabilities = numpy.ones(
            len(DNA_EXCHANGEABILITY_ORDER), dtype=float
        )
    else:
        normalized_exchangeabilities = normalize_dna_exchangeabilities_by_anchor(
            initial_exchangeabilities,
            model_name="GTR",
        )
    initial_values = {
        label: float(normalized_exchangeabilities[index + 1])
        for index, label in enumerate(_GTR_FREE_EXCHANGEABILITY_LABELS)
    }
    for parameter_name, parameter_value in initial_values.items():
        validate_parameter_within_bounds(
            parameter_name=parameter_name,
            value=parameter_value,
            lower_bound=validated_lower_bound,
            upper_bound=validated_upper_bound,
            owner_name="GTR exchangeability optimization",
        )
    bounds_by_name = dict.fromkeys(
        _GTR_FREE_EXCHANGEABILITY_LABELS,
        (validated_lower_bound, validated_upper_bound),
    )
    working_tree = tree.copy()
    validate_explicit_branch_lengths(working_tree, model_name="GTR")
    initial_report = _evaluate_gtr_tree_likelihood_from_patterns(
        working_tree,
        compressed_patterns,
        stationary_frequencies=stationary,
        base_frequency_source=source,
        exchangeabilities=_named_exchangeabilities_to_vector(initial_values),
        root_prior=stationary,
        observation_policy="reject",
        gap_state_frequency=None,
    )

    def evaluate_candidate_exchangeabilities(
        candidate_values: dict[str, float],
    ) -> tuple[GtrTreeLikelihoodReport, float]:
        report = _evaluate_gtr_tree_likelihood_from_patterns(
            working_tree,
            compressed_patterns,
            stationary_frequencies=stationary,
            base_frequency_source=source,
            exchangeabilities=_named_exchangeabilities_to_vector(candidate_values),
            root_prior=stationary,
            observation_policy="reject",
            gap_state_frequency=None,
        )
        return report, report.log_likelihood

    search_result = run_bounded_coordinate_likelihood_search(
        initial_values=initial_values,
        bounds_by_name=bounds_by_name,
        evaluate=evaluate_candidate_exchangeabilities,
        improvement_tolerance=improvement_tolerance,
        max_coordinate_passes=max_coordinate_passes,
    )
    optimized_report = search_result.payload
    return GtrExchangeabilityOptimizationReport(
        taxa=optimized_report.taxa,
        site_count=optimized_report.site_count,
        pattern_count=optimized_report.pattern_count,
        tree_newick=optimized_report.tree_newick,
        base_frequency_source=optimized_report.base_frequency_source,
        base_frequency_a=optimized_report.base_frequency_a,
        base_frequency_c=optimized_report.base_frequency_c,
        base_frequency_g=optimized_report.base_frequency_g,
        base_frequency_t=optimized_report.base_frequency_t,
        exchangeability_anchor=optimized_report.exchangeability_anchor,
        exchangeability_ac=optimized_report.exchangeability_ac,
        exchangeability_ag=optimized_report.exchangeability_ag,
        exchangeability_at=optimized_report.exchangeability_at,
        exchangeability_cg=optimized_report.exchangeability_cg,
        exchangeability_ct=optimized_report.exchangeability_ct,
        exchangeability_gt=optimized_report.exchangeability_gt,
        parameter_count=optimized_report.parameter_count,
        initial_log_likelihood=initial_report.log_likelihood,
        optimized_log_likelihood=optimized_report.log_likelihood,
        initial_aic=initial_report.aic,
        optimized_aic=optimized_report.aic,
        function_evaluation_count=search_result.function_evaluation_count,
        optimization_pass_count=search_result.optimization_pass_count,
        converged=search_result.converged,
        lower_exchangeability_bound=validated_lower_bound,
        upper_exchangeability_bound=validated_upper_bound,
    )


def optimize_gtr_exchangeabilities_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    initial_exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None,
    lower_exchangeability_bound: float = 0.05,
    upper_exchangeability_bound: float = 20.0,
    max_coordinate_passes: int = 24,
    improvement_tolerance: float = 1e-5,
) -> GtrExchangeabilityOptimizationReport:
    """Optimize GTR exchangeabilities from one tree path and one alignment path."""
    return optimize_gtr_exchangeabilities(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        base_frequencies=base_frequencies,
        initial_exchangeabilities=initial_exchangeabilities,
        lower_exchangeability_bound=lower_exchangeability_bound,
        upper_exchangeability_bound=upper_exchangeability_bound,
        max_coordinate_passes=max_coordinate_passes,
        improvement_tolerance=improvement_tolerance,
    )


def _evaluate_gtr_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    stationary_frequencies: numpy.ndarray,
    base_frequency_source: str,
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
    ),
    root_prior: numpy.ndarray,
    observation_policy: str,
    gap_state_frequency: float | None,
) -> GtrTreeLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name="GTR")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="GTR",
    )
    validated_frequencies = validate_dna_base_frequencies(
        stationary_frequencies,
        model_name="GTR",
    )
    normalized_exchangeabilities = normalize_dna_exchangeabilities_by_anchor(
        exchangeabilities,
        model_name="GTR",
    )
    state_order = dna_observation_state_order(observation_policy=observation_policy)
    if observation_policy == "fifth-state":
        if gap_state_frequency is None:
            raise ValueError(
                "GTR fifth-state observation policy requires an explicit gap_state_frequency"
            )
        augmented_rate_matrix = augment_dna_rate_matrix_with_gap_state(
            gtr_rate_matrix(
                exchangeabilities=normalized_exchangeabilities,
                base_frequencies=validated_frequencies,
            ),
            nucleotide_frequencies=validated_frequencies,
            gap_state_frequency=gap_state_frequency,
            model_name="GTR",
        )
        transition_by_node_id = {
            child.node_id: transition_probability_matrix(
                augmented_rate_matrix,
                max(float(child.branch_length or 0.0), 0.0),
            )
            for _parent, child in tree.iter_edges()
        }
    else:
        transition_by_node_id = {
            child.node_id: gtr_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                exchangeabilities=normalized_exchangeabilities,
                base_frequencies=validated_frequencies,
            )
            for _parent, child in tree.iter_edges()
        }

    def site_log_likelihood(states: tuple[str, ...]) -> float:
        states_by_taxon = dict(
            zip(compressed_patterns.taxon_order, states, strict=True)
        )
        pruning_pass = postorder_conditional_likelihoods(
            tree,
            state_count=len(state_order),
            leaf_likelihood=lambda node: resolve_dna_observation_leaf_vector(
                states_by_taxon,
                model_name="GTR",
                node_name=node.name,
                observation_policy=observation_policy,
            ),
            transition_matrix_for_child=lambda child: transition_by_node_id[
                child.node_id or ""
            ],
        )
        return log_likelihood_from_root_prior(
            tree,
            pruning_pass,
            root_prior=root_prior,
        )

    log_likelihood = sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=site_log_likelihood,
    )
    return GtrTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=True,
        tree_newick=dumps_newick(tree),
        state_count=len(state_order),
        observation_policy=observation_policy,
        base_frequency_source=base_frequency_source,
        base_frequency_a=float(validated_frequencies[0]),
        base_frequency_c=float(validated_frequencies[1]),
        base_frequency_g=float(validated_frequencies[2]),
        base_frequency_t=float(validated_frequencies[3]),
        exchangeability_anchor=_GTR_EXCHANGEABILITY_ANCHOR,
        exchangeability_ac=float(normalized_exchangeabilities[0]),
        exchangeability_ag=float(normalized_exchangeabilities[1]),
        exchangeability_at=float(normalized_exchangeabilities[2]),
        exchangeability_cg=float(normalized_exchangeabilities[3]),
        exchangeability_ct=float(normalized_exchangeabilities[4]),
        exchangeability_gt=float(normalized_exchangeabilities[5]),
        parameter_count=8,
        log_likelihood=log_likelihood,
        aic=(-2.0 * log_likelihood) + (2.0 * 8.0),
    )


def _named_exchangeabilities_to_vector(
    named_values: dict[str, float],
) -> numpy.ndarray:
    return numpy.array(
        [
            1.0,
            float(named_values["AG"]),
            float(named_values["AT"]),
            float(named_values["CG"]),
            float(named_values["CT"]),
            float(named_values["GT"]),
        ],
        dtype=float,
    )

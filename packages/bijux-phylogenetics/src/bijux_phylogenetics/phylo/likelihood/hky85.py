from __future__ import annotations

from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    is_dna_transition,
    normalize_dna_rate_matrix,
    validate_dna_base_frequencies,
    validate_positive_kappa,
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
from bijux_phylogenetics.phylo.likelihood.models import (
    Hky85KappaOptimizationReport,
    Hky85TreeLikelihoodReport,
)
from bijux_phylogenetics.phylo.likelihood.parameter_bounds import (
    validate_parameter_within_bounds,
    validate_positive_parameter_bounds,
)
from bijux_phylogenetics.phylo.likelihood.parameter_search import (
    run_bounded_likelihood_search,
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


def hky85_rate_matrix(
    base_frequencies: dict[str, float] | numpy.ndarray,
    *,
    kappa: float,
) -> numpy.ndarray:
    """Return the normalized HKY85 rate matrix with expected rate one."""
    stationary = validate_dna_base_frequencies(base_frequencies, model_name="HKY85")
    validated_kappa = validate_positive_kappa(kappa, model_name="HKY85")
    off_diagonal_rates = numpy.zeros((4, 4), dtype=float)
    for left_index, left_state in enumerate(DNA_STATE_ORDER):
        for right_index, right_state in enumerate(DNA_STATE_ORDER):
            if left_index == right_index:
                continue
            multiplier = (
                validated_kappa if is_dna_transition(left_state, right_state) else 1.0
            )
            off_diagonal_rates[left_index, right_index] = (
                multiplier * stationary[right_index]
            )
    return normalize_dna_rate_matrix(
        off_diagonal_rates,
        stationary_frequencies=stationary,
        model_name="HKY85",
    )


def hky85_transition_probability_matrix(
    branch_length: float,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray,
    kappa: float,
) -> numpy.ndarray:
    """Return the native HKY85 transition matrix for one branch."""
    validate_positive_kappa(kappa, model_name="HKY85")
    if branch_length <= 0.0:
        return numpy.eye(4, dtype=float)
    return transition_probability_matrix(
        hky85_rate_matrix(base_frequencies, kappa=kappa),
        branch_length,
    )


def evaluate_hky85_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    kappa: float,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> Hky85TreeLikelihoodReport:
    """Evaluate one fixed-topology HKY85 likelihood from aligned DNA records."""
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="HKY85",
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if base_frequencies is None:
        stationary = estimate_empirical_dna_base_frequencies_from_records(
            normalized_records,
            model_name="HKY85",
            observation_policy=observation_policy,
        )
        source = "estimated"
    else:
        stationary = validate_dna_base_frequencies(
            base_frequencies,
            model_name="HKY85",
        )
        source = "provided"
    gap_state_frequency = (
        estimate_empirical_gap_state_frequency(
            normalized_records,
            model_name="HKY85",
        )
        if observation_policy == "fifth-state"
        else None
    )
    resolved_root_prior = resolve_default_dna_root_prior_for_observation_policy(
        normalized_records,
        owner_name="HKY85 likelihood",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=stationary,
        observation_policy=observation_policy,
    )
    return _evaluate_hky85_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        stationary_frequencies=stationary,
        base_frequency_source=source,
        kappa=kappa,
        root_prior=resolved_root_prior.root_prior,
        observation_policy=observation_policy,
        gap_state_frequency=gap_state_frequency,
    )


def evaluate_hky85_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    kappa: float,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> Hky85TreeLikelihoodReport:
    """Evaluate one fixed-topology HKY85 likelihood from one tree path and alignment."""
    return evaluate_hky85_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        kappa=kappa,
        base_frequencies=base_frequencies,
        observation_policy=observation_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def optimize_hky85_kappa(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    initial_kappa: float = 1.0,
    lower_kappa_bound: float = 0.05,
    upper_kappa_bound: float = 20.0,
) -> Hky85KappaOptimizationReport:
    """Optimize one fixed-topology HKY85 kappa on fixed branch lengths."""
    validate_positive_kappa(initial_kappa, model_name="HKY85")
    validated_lower_bound, validated_upper_bound = validate_positive_parameter_bounds(
        parameter_name="kappa",
        lower_bound=lower_kappa_bound,
        upper_bound=upper_kappa_bound,
        owner_name="HKY85 kappa optimization",
    )
    validated_initial_kappa = validate_parameter_within_bounds(
        parameter_name="kappa",
        value=initial_kappa,
        lower_bound=validated_lower_bound,
        upper_bound=validated_upper_bound,
        owner_name="HKY85 kappa optimization",
    )

    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="HKY85",
        observation_policy="reject",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if base_frequencies is None:
        stationary = estimate_empirical_dna_base_frequencies_from_records(
            normalized_records,
            model_name="HKY85",
            observation_policy="reject",
        )
        source = "estimated"
    else:
        stationary = validate_dna_base_frequencies(
            base_frequencies,
            model_name="HKY85",
        )
        source = "provided"
    working_tree = tree.copy()
    validate_explicit_branch_lengths(working_tree, model_name="HKY85")
    initial_report = _evaluate_hky85_tree_likelihood_from_patterns(
        working_tree,
        compressed_patterns,
        stationary_frequencies=stationary,
        base_frequency_source=source,
        kappa=validated_initial_kappa,
        root_prior=stationary,
        observation_policy="reject",
        gap_state_frequency=None,
    )

    def evaluate_candidate_kappa(
        candidate_kappa: float,
    ) -> tuple[Hky85TreeLikelihoodReport, float]:
        report = _evaluate_hky85_tree_likelihood_from_patterns(
            working_tree,
            compressed_patterns,
            stationary_frequencies=stationary,
            base_frequency_source=source,
            kappa=candidate_kappa,
            root_prior=stationary,
            observation_policy="reject",
            gap_state_frequency=None,
        )
        return report, report.log_likelihood

    search_result = run_bounded_likelihood_search(
        lower_bound=validated_lower_bound,
        upper_bound=validated_upper_bound,
        evaluate=evaluate_candidate_kappa,
    )
    optimized_report = search_result.payload
    return Hky85KappaOptimizationReport(
        taxa=optimized_report.taxa,
        site_count=optimized_report.site_count,
        pattern_count=optimized_report.pattern_count,
        tree_newick=optimized_report.tree_newick,
        base_frequency_source=optimized_report.base_frequency_source,
        base_frequency_a=optimized_report.base_frequency_a,
        base_frequency_c=optimized_report.base_frequency_c,
        base_frequency_g=optimized_report.base_frequency_g,
        base_frequency_t=optimized_report.base_frequency_t,
        initial_kappa=validated_initial_kappa,
        optimized_kappa=search_result.parameter_value,
        parameter_count=optimized_report.parameter_count,
        initial_log_likelihood=initial_report.log_likelihood,
        optimized_log_likelihood=optimized_report.log_likelihood,
        initial_aic=initial_report.aic,
        optimized_aic=optimized_report.aic,
        function_evaluation_count=search_result.function_evaluation_count + 1,
        converged=search_result.converged,
        lower_kappa_bound=validated_lower_bound,
        upper_kappa_bound=validated_upper_bound,
    )


def optimize_hky85_kappa_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    initial_kappa: float = 1.0,
    lower_kappa_bound: float = 0.05,
    upper_kappa_bound: float = 20.0,
) -> Hky85KappaOptimizationReport:
    """Optimize HKY85 kappa from one tree path and one alignment path."""
    return optimize_hky85_kappa(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        base_frequencies=base_frequencies,
        initial_kappa=initial_kappa,
        lower_kappa_bound=lower_kappa_bound,
        upper_kappa_bound=upper_kappa_bound,
    )


def _evaluate_hky85_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    stationary_frequencies: numpy.ndarray,
    base_frequency_source: str,
    kappa: float,
    root_prior: numpy.ndarray,
    observation_policy: str,
    gap_state_frequency: float | None,
) -> Hky85TreeLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name="HKY85")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="HKY85",
    )
    validated_frequencies = validate_dna_base_frequencies(
        stationary_frequencies,
        model_name="HKY85",
    )
    validated_kappa = validate_positive_kappa(kappa, model_name="HKY85")
    state_order = dna_observation_state_order(observation_policy=observation_policy)
    if observation_policy == "fifth-state":
        if gap_state_frequency is None:
            raise ValueError(
                "HKY85 fifth-state observation policy requires an explicit gap_state_frequency"
            )
        augmented_rate_matrix = augment_dna_rate_matrix_with_gap_state(
            hky85_rate_matrix(
                validated_frequencies,
                kappa=validated_kappa,
            ),
            nucleotide_frequencies=validated_frequencies,
            gap_state_frequency=gap_state_frequency,
            model_name="HKY85",
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
            child.node_id: hky85_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
                base_frequencies=validated_frequencies,
                kappa=validated_kappa,
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
                model_name="HKY85",
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
    return Hky85TreeLikelihoodReport(
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
        kappa=validated_kappa,
        parameter_count=4,
        log_likelihood=log_likelihood,
        aic=(-2.0 * log_likelihood) + (2.0 * 4.0),
    )

from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.dna import (
    DNA_STATE_ORDER,
    validate_dna_base_frequencies,
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
from bijux_phylogenetics.phylo.likelihood.models import F81TreeLikelihoodReport
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
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


def f81_rate_matrix(
    base_frequencies: dict[str, float] | numpy.ndarray,
) -> numpy.ndarray:
    """Return the normalized F81 rate matrix with expected rate one."""
    stationary = validate_dna_base_frequencies(base_frequencies, model_name="F81")
    variability = _f81_variability(stationary)
    if variability <= 0.0:
        raise InvalidAlignmentError(
            "F81 likelihood requires at least two nucleotides with nonzero stationary frequency"
        )
    rate_matrix = numpy.zeros((4, 4), dtype=float)
    for row_index in range(len(DNA_STATE_ORDER)):
        row_total = 0.0
        for column_index in range(len(DNA_STATE_ORDER)):
            if row_index == column_index:
                continue
            rate = stationary[column_index] / variability
            rate_matrix[row_index, column_index] = rate
            row_total += rate
        rate_matrix[row_index, row_index] = -row_total
    return rate_matrix


def f81_transition_probability_matrix(
    branch_length: float,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray,
) -> numpy.ndarray:
    """Return the native closed-form F81 transition matrix for one branch."""
    stationary = validate_dna_base_frequencies(base_frequencies, model_name="F81")
    if branch_length <= 0.0:
        return numpy.eye(4, dtype=float)
    variability = _f81_variability(stationary)
    if variability <= 0.0:
        raise InvalidAlignmentError(
            "F81 likelihood requires at least two nucleotides with nonzero stationary frequency"
        )
    decay = math.exp(-branch_length / variability)
    transition_matrix = numpy.zeros((4, 4), dtype=float)
    for row_index in range(len(DNA_STATE_ORDER)):
        for column_index in range(len(DNA_STATE_ORDER)):
            if row_index == column_index:
                transition_matrix[row_index, column_index] = stationary[
                    column_index
                ] + ((1.0 - stationary[column_index]) * decay)
            else:
                transition_matrix[row_index, column_index] = stationary[
                    column_index
                ] * (1.0 - decay)
    return transition_matrix


def evaluate_f81_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> F81TreeLikelihoodReport:
    """Evaluate one fixed-topology F81 likelihood from aligned DNA records."""
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="F81",
        observation_policy=observation_policy,
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    if base_frequencies is None:
        stationary = estimate_empirical_dna_base_frequencies_from_records(
            normalized_records,
            model_name="F81",
            observation_policy=observation_policy,
        )
        source = "estimated"
    else:
        stationary = validate_dna_base_frequencies(base_frequencies, model_name="F81")
        source = "provided"
    gap_state_frequency = (
        estimate_empirical_gap_state_frequency(
            normalized_records,
            model_name="F81",
        )
        if observation_policy == "fifth-state"
        else None
    )
    resolved_root_prior = resolve_default_dna_root_prior_for_observation_policy(
        normalized_records,
        owner_name="F81 likelihood",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=stationary,
        observation_policy=observation_policy,
    )
    return _evaluate_f81_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        stationary_frequencies=stationary,
        base_frequency_source=source,
        root_prior=resolved_root_prior.root_prior,
        observation_policy=observation_policy,
        gap_state_frequency=gap_state_frequency,
    )


def evaluate_f81_tree_likelihood_from_alignment(
    tree_path: Path,
    alignment_path: Path,
    *,
    base_frequencies: dict[str, float] | numpy.ndarray | None = None,
    observation_policy: str = "reject",
    root_prior_policy: str | None = None,
    root_prior: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    fixed_root_state: str | None = None,
) -> F81TreeLikelihoodReport:
    """Evaluate one fixed-topology F81 likelihood from one tree path and alignment."""
    return evaluate_f81_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        base_frequencies=base_frequencies,
        observation_policy=observation_policy,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
    )


def _evaluate_f81_tree_likelihood_from_patterns(
    tree: PhyloTree,
    compressed_patterns: CompressedAlignmentSitePatterns,
    *,
    stationary_frequencies: numpy.ndarray,
    base_frequency_source: str,
    root_prior: numpy.ndarray,
    observation_policy: str,
    gap_state_frequency: float | None,
) -> F81TreeLikelihoodReport:
    validate_explicit_branch_lengths(tree, model_name="F81")
    validate_tree_taxa_against_patterns(
        tree,
        compressed_patterns,
        model_name="F81",
    )
    validated_frequencies = validate_dna_base_frequencies(
        stationary_frequencies,
        model_name="F81",
    )
    state_order = dna_observation_state_order(observation_policy=observation_policy)
    if observation_policy == "fifth-state":
        if gap_state_frequency is None:
            raise ValueError(
                "F81 fifth-state observation policy requires an explicit gap_state_frequency"
            )
        augmented_rate_matrix = augment_dna_rate_matrix_with_gap_state(
            f81_rate_matrix(validated_frequencies),
            nucleotide_frequencies=validated_frequencies,
            gap_state_frequency=gap_state_frequency,
            model_name="F81",
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
            child.node_id: f81_transition_probability_matrix(
                max(float(child.branch_length or 0.0), 0.0),
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
                model_name="F81",
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
    return F81TreeLikelihoodReport(
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
        parameter_count=3,
        log_likelihood=log_likelihood,
        aic=(-2.0 * log_likelihood) + (2.0 * 3.0),
    )


def _f81_variability(stationary_frequencies: numpy.ndarray) -> float:
    return 1.0 - float(numpy.sum(stationary_frequencies * stationary_frequencies))

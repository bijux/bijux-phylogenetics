from __future__ import annotations

from collections import OrderedDict
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
from bijux_phylogenetics.phylo.likelihood.models import CodonCtmcTreeLikelihoodReport
from bijux_phylogenetics.phylo.likelihood.patterns import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.pruning import (
    build_transition_matrix_evaluator,
    log_likelihood_from_root_prior,
    postorder_conditional_likelihoods,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    sum_compressed_site_pattern_log_likelihoods,
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


def evaluate_codon_ctmc_tree_likelihood(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    codon_frequencies: dict[str, float]
    | numpy.ndarray
    | list[float]
    | tuple[float, ...]
    | None = None,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> CodonCtmcTreeLikelihoodReport:
    """Evaluate one fixed-topology codon CTMC likelihood on an aligned coding matrix."""
    compressed_patterns = compress_codon_site_patterns_from_records(
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
    log_likelihood = _evaluate_codon_ctmc_tree_likelihood_from_patterns(
        tree,
        compressed_patterns,
        state_space=state_space,
        rate_matrix=rate_matrix,
        root_prior=frequencies,
        observation_policy=observation_policy,
    )
    return CodonCtmcTreeLikelihoodReport(
        taxa=compressed_patterns.taxon_order,
        site_count=compressed_patterns.alignment_length,
        pattern_count=compressed_patterns.pattern_count,
        compression_used=(
            compressed_patterns.pattern_count < compressed_patterns.alignment_length
        ),
        tree_newick=dumps_newick(tree),
        state_count=len(state_space.state_order),
        genetic_code_id=state_space.genetic_code_id,
        genetic_code_name=state_space.genetic_code_name,
        codon_frequency_source=frequency_source,
        observation_policy=validate_codon_observation_policy(
            observation_policy,
            owner_name="codon CTMC likelihood",
        ),
        log_likelihood=log_likelihood,
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
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
) -> CodonCtmcTreeLikelihoodReport:
    """Evaluate one fixed-topology codon CTMC likelihood from tree and alignment paths."""
    return evaluate_codon_ctmc_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        codon_frequencies=codon_frequencies,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )


def compress_codon_site_patterns_from_records(
    records: list[AlignmentRecord],
    *,
    genetic_code: int | str | None = None,
    observation_policy: str = "reject",
    source_path: Path | None = None,
) -> CompressedAlignmentSitePatterns:
    """Compress identical aligned codon sites while preserving stable taxon order."""
    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    codon_count = len(normalized_records[0].sequence) // 3
    codon_sequences = [
        _iter_record_codons(record.sequence) for record in normalized_records
    ]
    grouped_positions: OrderedDict[tuple[str, ...], list[int]] = OrderedDict()
    for site_index in range(codon_count):
        states = tuple(sequence[site_index] for sequence in codon_sequences)
        grouped_positions.setdefault(states, []).append(site_index + 1)
    return CompressedAlignmentSitePatterns(
        source_path=source_path,
        taxon_order=[record.identifier for record in normalized_records],
        alignment_length=codon_count,
        pattern_count=len(grouped_positions),
        patterns=[
            AlignmentSitePattern(
                pattern_id=f"pattern-{index}",
                states=states,
                weight=len(site_positions),
                site_positions=site_positions,
            )
            for index, (states, site_positions) in enumerate(
                grouped_positions.items(),
                start=1,
            )
        ],
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
    return sum_compressed_site_pattern_log_likelihoods(
        compressed_patterns,
        site_log_likelihood=lambda states: _evaluate_codon_site_log_likelihood(
            tree,
            states,
            taxon_order=compressed_patterns.taxon_order,
            state_space=state_space,
            root_prior=root_prior,
            transition_evaluator=transition_evaluator,
            observation_policy=observation_policy,
        ),
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


def _iter_record_codons(sequence: str) -> tuple[str, ...]:
    return tuple(sequence[index : index + 3] for index in range(0, len(sequence), 3))

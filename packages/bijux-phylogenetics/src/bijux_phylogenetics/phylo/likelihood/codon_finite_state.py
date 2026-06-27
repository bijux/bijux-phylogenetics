from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    slice_partition_sequence,
)
from bijux_phylogenetics.phylo.likelihood.codon_observation_policies import (
    resolve_codon_observation_leaf_vector,
    validate_codon_observation_policy,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix,
)
from bijux_phylogenetics.phylo.likelihood.ctmc import (
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CompiledFiniteStateLikelihoodAnalysis,
    FiniteStateLikelihoodSurface,
    compile_finite_state_likelihood_analysis,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.topology.tree import PhyloTree
from bijux_phylogenetics.runtime.errors import InvalidPartitionError

from .codon import normalize_codon_likelihood_records
from .codon import compress_codon_site_patterns_from_records
from .codon import resolve_codon_root_prior


def build_codon_finite_state_surface(
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
) -> FiniteStateLikelihoodSurface:
    """Build one shared finite-state codon CTMC likelihood surface."""
    normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    state_space, rate_matrix, frequencies, _frequency_source = (
        build_equal_rate_codon_ctmc_rate_matrix(
            codon_frequencies,
            genetic_code=genetic_code,
        )
    )
    resolved_root_prior = resolve_codon_root_prior(
        records,
        state_space=state_space,
        owner_name="codon CTMC finite-state surface",
        default_policy="stationary",
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        stationary_frequencies=frequencies,
    )
    validated_observation_policy = validate_codon_observation_policy(
        observation_policy,
        owner_name="codon CTMC likelihood",
    )
    return FiniteStateLikelihoodSurface(
        model_name="codon CTMC",
        state_labels=tuple(state_space.state_order),
        root_prior=resolved_root_prior.root_prior,
        root_prior_source=resolved_root_prior.root_prior_source,
        observation_policy=validated_observation_policy,
        normalization_policy=EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_codon_observation_leaf_vector(
                observed_state,
                model_name="codon CTMC",
                node_name=taxon_name,
                observation_policy=validated_observation_policy,
                state_order=state_space.state_order,
                state_index=state_space.state_index,
            )
        ),
    )


def compile_codon_likelihood_analysis(
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
    partition_map: tuple[LocusPartition, ...] | None = None,
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one reusable codon finite-state likelihood analysis."""
    normalized_records = normalize_codon_likelihood_records(
        records,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    surface = build_codon_finite_state_surface(
        normalized_records,
        codon_frequencies=codon_frequencies,
        root_prior_policy=root_prior_policy,
        root_prior=root_prior,
        fixed_root_state=fixed_root_state,
        genetic_code=genetic_code,
        observation_policy=observation_policy,
    )
    return compile_finite_state_likelihood_analysis(
        tree,
        surface,
        records=normalized_records,
        pattern_compiler=_compile_codon_patterns,
        missing_data_policy=surface.observation_policy,
        alignment_length_resolver=_resolve_codon_alignment_length,
        partition_map=partition_map,
        partition_alignment_length=len(normalized_records[0].sequence),
        partition_record_builder=_slice_codon_partition_records,
        partition_source_site_positions_resolver=(
            _resolve_codon_partition_site_positions
        ),
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def _compile_codon_patterns(
    records: list[AlignmentRecord],
    source_path: Path | None,
) -> CompressedAlignmentSitePatterns:
    return compress_codon_site_patterns_from_records(
        records,
        source_path=source_path,
    )


def _resolve_codon_alignment_length(records: list[AlignmentRecord]) -> int:
    return len(records[0].sequence) // 3


def _slice_codon_partition_records(
    records: list[AlignmentRecord],
    partition: LocusPartition,
) -> list[AlignmentRecord]:
    _validate_codon_partition(partition)
    sliced_records = [
        AlignmentRecord(
            identifier=record.identifier,
            sequence=slice_partition_sequence(record.sequence, partition),
        )
        for record in records
    ]
    sliced_lengths = {len(record.sequence) for record in sliced_records}
    if len(sliced_lengths) != 1:
        raise InvalidPartitionError(
            "codon CTMC partitions must yield equal-length coding slices"
        )
    sliced_length = next(iter(sliced_lengths))
    if sliced_length == 0 or sliced_length % 3 != 0:
        raise InvalidPartitionError(
            "codon CTMC partitions must preserve complete codon triplets"
        )
    return sliced_records


def _resolve_codon_partition_site_positions(
    partition: LocusPartition,
) -> tuple[int, ...]:
    _validate_codon_partition(partition)
    return tuple(
        ((segment_start - 1) // 3) + 1
        for segment in partition.segments
        for segment_start in range(segment.start, segment.end + 1, 3)
    )


def _validate_codon_partition(partition: LocusPartition) -> None:
    for segment in partition.segments:
        if segment.step != 1:
            raise InvalidPartitionError(
                f"codon CTMC partition '{partition.name}' must use contiguous codon blocks without stride"
            )
        if (segment.start - 1) % 3 != 0 or segment.end % 3 != 0:
            raise InvalidPartitionError(
                f"codon CTMC partition '{partition.name}' must start on codon base 1 and end on codon base 3"
            )

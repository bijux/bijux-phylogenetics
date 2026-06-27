from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    CompiledFiniteStateLikelihoodAnalysis,
    compile_partitioned_finite_state_likelihood_analysis,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.dna_observation_policies import (
    normalize_dna_likelihood_records,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree


@dataclass(frozen=True, slots=True)
class PartitionedNucleotideModelAssignment:
    """One explicit nucleotide likelihood model assignment for one named partition."""

    partition_name: str
    model_name: str
    context_linkage_key: str | None = None
    kappa: float | None = None
    base_frequencies: dict[str, float] | numpy.ndarray | None = None
    exchangeabilities: (
        dict[tuple[str, str], float]
        | dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None
    root_prior_policy: str | None = None
    root_prior: (
        dict[str, float]
        | numpy.ndarray
        | list[float]
        | tuple[float, ...]
        | None
    ) = None
    fixed_root_state: str | None = None


def compile_partitioned_nucleotide_likelihood_analysis(
    tree: PhyloTree,
    records: list[AlignmentRecord],
    *,
    partition_map: tuple[LocusPartition, ...],
    partition_model_assignments: Sequence[PartitionedNucleotideModelAssignment],
    observation_policy: str = "reject",
    site_weight_by_position: Mapping[int, float] | None = None,
    source_path: Path | None = None,
) -> CompiledFiniteStateLikelihoodAnalysis:
    """Compile one partitioned nucleotide likelihood analysis with linked or unlinked blocks."""
    normalized_records = normalize_dna_likelihood_records(
        records,
        model_name="partitioned nucleotide",
        observation_policy=observation_policy,
    )
    assignment_by_name = _resolve_partition_model_assignments(
        partition_map=partition_map,
        partition_model_assignments=partition_model_assignments,
    )
    return compile_partitioned_finite_state_likelihood_analysis(
        tree,
        records=normalized_records,
        pattern_compiler=lambda partition_records, partition_source_path: (
            compress_alignment_site_patterns_from_records(
                partition_records,
                source_path=partition_source_path,
            )
        ),
        partition_surface_builder=lambda partition_name, _data_type, partition_records: (
            _build_partition_surface(
                partition_records=partition_records,
                assignment=assignment_by_name[partition_name],
                observation_policy=observation_policy,
            )
        ),
        missing_data_policy=observation_policy,
        partition_map=partition_map,
        partition_context_linkage_key_resolver=(
            lambda partition_name, _data_type, _surface: (
                assignment_by_name[partition_name].context_linkage_key
                or partition_name
            )
        ),
        site_weight_by_position=site_weight_by_position,
        source_path=source_path,
    )


def _resolve_partition_model_assignments(
    *,
    partition_map: tuple[LocusPartition, ...],
    partition_model_assignments: Sequence[PartitionedNucleotideModelAssignment],
) -> dict[str, PartitionedNucleotideModelAssignment]:
    partition_names = tuple(partition.name for partition in partition_map)
    assignment_by_name: dict[str, PartitionedNucleotideModelAssignment] = {}
    for assignment in partition_model_assignments:
        if assignment.partition_name in assignment_by_name:
            raise ValueError(
                f"partitioned nucleotide likelihood received duplicate model assignments for partition '{assignment.partition_name}'"
            )
        assignment_by_name[assignment.partition_name] = assignment
    assignment_names = tuple(sorted(assignment_by_name))
    expected_names = tuple(sorted(partition_names))
    if assignment_names != expected_names:
        raise ValueError(
            "partitioned nucleotide likelihood requires one model assignment for every declared partition"
        )
    return assignment_by_name


def _build_partition_surface(
    *,
    partition_records: list[AlignmentRecord],
    assignment: PartitionedNucleotideModelAssignment,
    observation_policy: str,
):
    return build_selected_nucleotide_finite_state_surface(
        partition_records,
        model_name=assignment.model_name,
        owner_name=f"partitioned nucleotide partition '{assignment.partition_name}'",
        observation_policy=observation_policy,
        kappa=assignment.kappa,
        base_frequencies=assignment.base_frequencies,
        exchangeabilities=assignment.exchangeabilities,
        root_prior_policy=assignment.root_prior_policy,
        root_prior=assignment.root_prior,
        fixed_root_state=assignment.fixed_root_state,
    )

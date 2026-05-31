from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionSubstitutionModelDefinition,
)
from bijux_phylogenetics.bayesian.state import BayesianPhylogeneticState
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    validate_locus_partitions,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_DNA_STATES = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class PosteriorAncestralSequenceDefinition:
    """Validated nucleotide posterior ancestral-sequence summary configuration."""

    records: tuple[AlignmentRecord, ...]
    posterior_probability_threshold: float
    minimum_clade_posterior_support: float
    low_confidence_state_symbol: str
    locus_partitions: tuple[LocusPartition, ...] | None = None
    partition_models: tuple[PartitionSubstitutionModelDefinition, ...] | None = None


@dataclass(frozen=True, slots=True)
class PosteriorAncestralStateProbabilityRow:
    """One clade-site-state posterior probability aggregated across sampled states."""

    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    site_position: int
    state: str
    clade_posterior_probability: float
    conditional_posterior_probability: float
    marginal_posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorAncestralSiteSummaryRow:
    """One clade-site posterior summary across all DNA states."""

    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    site_position: int
    consensus_state: str
    exported_state: str
    clade_posterior_probability: float
    max_conditional_posterior_probability: float
    max_marginal_posterior_probability: float
    low_confidence: bool
    posterior_probability_a: float
    posterior_probability_c: float
    posterior_probability_g: float
    posterior_probability_t: float


@dataclass(frozen=True, slots=True)
class PosteriorAncestralSequenceRecord:
    """One consensus ancestral sequence record for one posterior-supported clade."""

    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: list[str]
    clade_posterior_probability: float
    sequence: str


@dataclass(frozen=True, slots=True)
class PosteriorAncestralSequenceReport:
    """Posterior ancestral-sequence summary across one sampled Bayesian chain."""

    sample_count: int
    site_count: int
    taxon_count: int
    distinct_topology_count: int
    sampled_substitution_models: list[str]
    tree_uncertainty_policy: str
    posterior_probability_threshold: float
    minimum_clade_posterior_support: float
    low_confidence_state_symbol: str
    state_probability_rows: list[PosteriorAncestralStateProbabilityRow]
    site_summary_rows: list[PosteriorAncestralSiteSummaryRow]
    sequence_records: list[PosteriorAncestralSequenceRecord]
    warnings: list[str]


def build_posterior_ancestral_sequence_definition(
    *,
    records: list[AlignmentRecord] | tuple[AlignmentRecord, ...],
    posterior_probability_threshold: float = 0.5,
    minimum_clade_posterior_support: float = 0.5,
    low_confidence_state_symbol: str = "N",
    locus_partitions: list[LocusPartition] | tuple[LocusPartition, ...] | None = None,
    partition_models: (
        list[PartitionSubstitutionModelDefinition]
        | tuple[PartitionSubstitutionModelDefinition, ...]
        | None
    ) = None,
) -> PosteriorAncestralSequenceDefinition:
    """Build one validated posterior ancestral-sequence summary definition."""
    validated_records = _validate_alignment_records(records)
    validated_probability_threshold = _validate_probability_threshold(
        value=posterior_probability_threshold,
        field_name="posterior_probability_threshold",
        owner_name="posterior ancestral sequence definition",
    )
    validated_clade_support_threshold = _validate_probability_threshold(
        value=minimum_clade_posterior_support,
        field_name="minimum_clade_posterior_support",
        owner_name="posterior ancestral sequence definition",
    )
    validated_low_confidence_symbol = _validate_low_confidence_state_symbol(
        low_confidence_state_symbol
    )
    validated_locus_partitions, validated_partition_models = (
        _validate_optional_partition_surface(
            locus_partitions=locus_partitions,
            partition_models=partition_models,
            records=validated_records,
        )
    )
    return PosteriorAncestralSequenceDefinition(
        records=validated_records,
        posterior_probability_threshold=validated_probability_threshold,
        minimum_clade_posterior_support=validated_clade_support_threshold,
        low_confidence_state_symbol=validated_low_confidence_symbol,
        locus_partitions=validated_locus_partitions,
        partition_models=validated_partition_models,
    )


def summarize_nucleotide_posterior_ancestral_sequences(
    sampled_states: list[BayesianPhylogeneticState] | tuple[BayesianPhylogeneticState, ...],
    *,
    definition: PosteriorAncestralSequenceDefinition,
) -> PosteriorAncestralSequenceReport:
    """Summarize nucleotide ancestral sequence probabilities across posterior samples."""
    if not isinstance(definition, PosteriorAncestralSequenceDefinition):
        raise PhylogeneticsError(
            "posterior ancestral sequence summary requires one PosteriorAncestralSequenceDefinition",
            code="posterior_ancestral_sequence_definition_type_invalid",
        )
    validated_sampled_states = _validate_sampled_states(sampled_states)
    raise NotImplementedError


def write_posterior_ancestral_sequence_fasta(
    path: Path,
    report: PosteriorAncestralSequenceReport,
) -> Path:
    """Write posterior consensus ancestral sequences for posterior-supported clades."""
    return write_fasta_alignment(
        path,
        [
            AlignmentRecord(
                identifier=record.clade_id,
                sequence=record.sequence,
            )
            for record in report.sequence_records
        ],
    )


def write_posterior_ancestral_state_probability_table(
    path: Path,
    report: PosteriorAncestralSequenceReport,
) -> Path:
    """Write one TSV table of posterior clade-site-state probabilities."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "clade_id",
                "representative_node_id",
                "node_name",
                "descendant_taxa",
                "site_position",
                "state",
                "clade_posterior_probability",
                "conditional_posterior_probability",
                "marginal_posterior_probability",
                "supporting_sample_count",
                "total_sample_count",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for row in report.state_probability_rows:
            writer.writerow(
                {
                    "clade_id": row.clade_id,
                    "representative_node_id": row.representative_node_id,
                    "node_name": row.node_name or "",
                    "descendant_taxa": "|".join(row.descendant_taxa),
                    "site_position": row.site_position,
                    "state": row.state,
                    "clade_posterior_probability": repr(
                        row.clade_posterior_probability
                    ),
                    "conditional_posterior_probability": repr(
                        row.conditional_posterior_probability
                    ),
                    "marginal_posterior_probability": repr(
                        row.marginal_posterior_probability
                    ),
                    "supporting_sample_count": row.supporting_sample_count,
                    "total_sample_count": row.total_sample_count,
                }
            )
    return path


def _validate_alignment_records(
    records: list[AlignmentRecord] | tuple[AlignmentRecord, ...],
) -> tuple[AlignmentRecord, ...]:
    validated_records = tuple(records)
    if not validated_records:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires at least one alignment record",
            code="posterior_ancestral_sequence_records_empty",
        )
    if any(not isinstance(record, AlignmentRecord) for record in validated_records):
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires every alignment record to be one AlignmentRecord",
            code="posterior_ancestral_sequence_record_type_invalid",
        )
    sequence_lengths = {len(record.sequence) for record in validated_records}
    if len(sequence_lengths) != 1:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires one aligned sequence matrix with equal sequence lengths",
            code="posterior_ancestral_sequence_alignment_length_invalid",
        )
    return validated_records


def _validate_probability_threshold(
    *,
    value: float,
    field_name: str,
    owner_name: str,
) -> float:
    if not 0.0 < value <= 1.0:
        raise PhylogeneticsError(
            f"{owner_name} requires '{field_name}' to be greater than zero and at most one",
            code="posterior_ancestral_sequence_probability_threshold_invalid",
            details={field_name: value},
        )
    return float(format(value, ".15g"))


def _validate_low_confidence_state_symbol(low_confidence_state_symbol: str) -> str:
    normalized_symbol = low_confidence_state_symbol.strip().upper()
    if len(normalized_symbol) != 1:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires one single-character low-confidence state symbol",
            code="posterior_ancestral_sequence_low_confidence_symbol_invalid",
            details={"low_confidence_state_symbol": low_confidence_state_symbol},
        )
    return normalized_symbol


def _validate_optional_partition_surface(
    *,
    locus_partitions: list[LocusPartition] | tuple[LocusPartition, ...] | None,
    partition_models: (
        list[PartitionSubstitutionModelDefinition]
        | tuple[PartitionSubstitutionModelDefinition, ...]
        | None
    ),
    records: tuple[AlignmentRecord, ...],
) -> tuple[
    tuple[LocusPartition, ...] | None,
    tuple[PartitionSubstitutionModelDefinition, ...] | None,
]:
    if locus_partitions is None and partition_models is None:
        return None, None
    if locus_partitions is None or partition_models is None:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires locus partitions and partition models together when summarizing partitioned DNA samples",
            code="posterior_ancestral_sequence_partition_surface_incomplete",
        )
    validated_locus_partitions = tuple(locus_partitions)
    validated_partition_models = tuple(partition_models)
    if not validated_locus_partitions or not validated_partition_models:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires non-empty partition surfaces when partitioned DNA support is configured",
            code="posterior_ancestral_sequence_partition_surface_empty",
        )
    alignment_length = len(records[0].sequence)
    _assigned_site_count, unassigned_site_count = validate_locus_partitions(
        validated_locus_partitions,
        alignment_length=alignment_length,
    )
    if unassigned_site_count != 0:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires partition loci to cover every alignment site",
            code="posterior_ancestral_sequence_partition_coverage_incomplete",
            details={"unassigned_site_count": unassigned_site_count},
        )
    locus_partition_names = tuple(partition.name for partition in validated_locus_partitions)
    partition_model_names = tuple(
        partition_model.partition_name for partition_model in validated_partition_models
    )
    if locus_partition_names != partition_model_names:
        raise PhylogeneticsError(
            "posterior ancestral sequence definition requires locus partitions and partition models to use the same ordered partition names",
            code="posterior_ancestral_sequence_partition_names_mismatched",
            details={
                "locus_partition_names": list(locus_partition_names),
                "partition_model_names": list(partition_model_names),
            },
        )
    return validated_locus_partitions, validated_partition_models


def _validate_sampled_states(
    sampled_states: list[BayesianPhylogeneticState] | tuple[BayesianPhylogeneticState, ...],
) -> tuple[BayesianPhylogeneticState, ...]:
    validated_sampled_states = tuple(sampled_states)
    if not validated_sampled_states:
        raise PhylogeneticsError(
            "posterior ancestral sequence summary requires at least one sampled Bayesian phylogenetic state",
            code="posterior_ancestral_sequence_sampled_states_empty",
        )
    if any(
        not isinstance(sampled_state, BayesianPhylogeneticState)
        for sampled_state in validated_sampled_states
    ):
        raise PhylogeneticsError(
            "posterior ancestral sequence summary requires every sampled state to be one BayesianPhylogeneticState",
            code="posterior_ancestral_sequence_sampled_state_type_invalid",
        )
    return validated_sampled_states


__all__ = [
    "PosteriorAncestralSequenceDefinition",
    "PosteriorAncestralSequenceRecord",
    "PosteriorAncestralSequenceReport",
    "PosteriorAncestralSiteSummaryRow",
    "PosteriorAncestralStateProbabilityRow",
    "build_posterior_ancestral_sequence_definition",
    "summarize_nucleotide_posterior_ancestral_sequences",
    "write_posterior_ancestral_sequence_fasta",
    "write_posterior_ancestral_state_probability_table",
]

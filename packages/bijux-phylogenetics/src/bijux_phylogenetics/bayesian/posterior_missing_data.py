from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import LocusPartition
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

from .partition_model_priors import PartitionSubstitutionModelDefinition
from .posterior_ancestral_sequences import (
    _validate_alignment_records,
    _validate_low_confidence_state_symbol,
    _validate_optional_partition_surface,
    _validate_probability_threshold,
)

_DNA_STATES = ("A", "C", "G", "T")


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideDefinition:
    """Validated posterior missing-state summary configuration for DNA alignments."""

    records: tuple[AlignmentRecord, ...]
    missing_state_symbols: tuple[str, ...]
    posterior_probability_threshold: float
    low_confidence_state_symbol: str
    locus_partitions: tuple[LocusPartition, ...] | None = None
    partition_models: tuple[PartitionSubstitutionModelDefinition, ...] | None = None


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideStateProbabilityRow:
    """One taxon-site-state posterior probability for one masked nucleotide."""

    taxon: str
    site_position: int
    state: str
    posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideSiteSummaryRow:
    """One taxon-site posterior summary across all DNA states."""

    taxon: str
    site_position: int
    observed_symbol: str
    consensus_state: str
    exported_state: str
    max_posterior_probability: float
    low_confidence: bool
    posterior_probability_a: float
    posterior_probability_c: float
    posterior_probability_g: float
    posterior_probability_t: float


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideSequenceRecord:
    """One imputed nucleotide sequence for one taxon with masked tip sites."""

    identifier: str
    masked_site_count: int
    sequence: str


@dataclass(frozen=True, slots=True)
class PosteriorMissingNucleotideReport:
    """Posterior missing-state summary for masked nucleotide observations."""

    sample_count: int
    site_count: int
    taxon_count: int
    masked_site_count: int
    distinct_topology_count: int
    sampled_substitution_models: list[str]
    state_probability_rows: list[PosteriorMissingNucleotideStateProbabilityRow]
    site_summary_rows: list[PosteriorMissingNucleotideSiteSummaryRow]
    sequence_records: list[PosteriorMissingNucleotideSequenceRecord]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitStateProbabilityRow:
    """One taxon-state posterior probability for one masked discrete trait."""

    taxon: str
    state: str
    posterior_probability: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitTaxonSummaryRow:
    """One taxon-level posterior summary for one masked discrete trait."""

    taxon: str
    observed_symbol: str
    most_likely_state: str
    max_posterior_probability: float
    posterior_entropy: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingDiscreteTraitReport:
    """Posterior missing-state summary for masked discrete traits."""

    sample_count: int
    distinct_topology_count: int
    sampled_transition_models: list[str]
    state_order: list[str]
    state_probability_rows: list[PosteriorMissingDiscreteTraitStateProbabilityRow]
    taxon_summary_rows: list[PosteriorMissingDiscreteTraitTaxonSummaryRow]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class PosteriorMissingContinuousTraitTaxonSummaryRow:
    """One taxon-level posterior summary for one masked continuous trait."""

    taxon: str
    posterior_mean: float
    posterior_median: float
    posterior_hpd_95_lower: float
    posterior_hpd_95_upper: float
    mean_conditional_standard_deviation: float
    supporting_sample_count: int
    total_sample_count: int


@dataclass(frozen=True, slots=True)
class PosteriorMissingContinuousTraitReport:
    """Posterior missing-value summary for masked continuous traits."""

    sample_count: int
    distinct_topology_count: int
    sampled_trait_models: list[str]
    taxon_summary_rows: list[PosteriorMissingContinuousTraitTaxonSummaryRow]
    warnings: list[str]


def build_posterior_missing_nucleotide_definition(
    *,
    records: list[AlignmentRecord] | tuple[AlignmentRecord, ...],
    missing_state_symbols: tuple[str, ...] | list[str] = ("N", "?"),
    posterior_probability_threshold: float = 0.5,
    low_confidence_state_symbol: str = "N",
    locus_partitions: list[LocusPartition] | tuple[LocusPartition, ...] | None = None,
    partition_models: (
        list[PartitionSubstitutionModelDefinition]
        | tuple[PartitionSubstitutionModelDefinition, ...]
        | None
    ) = None,
) -> PosteriorMissingNucleotideDefinition:
    """Build one validated posterior missing-nucleotide summary definition."""
    validated_records = _validate_alignment_records(records)
    validated_missing_symbols = _validate_missing_state_symbols(missing_state_symbols)
    validated_locus_partitions, validated_partition_models = (
        _validate_optional_partition_surface(
            locus_partitions=locus_partitions,
            partition_models=partition_models,
            records=validated_records,
        )
    )
    return PosteriorMissingNucleotideDefinition(
        records=validated_records,
        missing_state_symbols=validated_missing_symbols,
        posterior_probability_threshold=_validate_probability_threshold(
            value=posterior_probability_threshold,
            field_name="posterior_probability_threshold",
            owner_name="posterior missing-nucleotide definition",
        ),
        low_confidence_state_symbol=_validate_low_confidence_state_symbol(
            low_confidence_state_symbol
        ),
        locus_partitions=validated_locus_partitions,
        partition_models=validated_partition_models,
    )


def _validate_missing_state_symbols(
    missing_state_symbols: tuple[str, ...] | list[str],
) -> tuple[str, ...]:
    validated_symbols = tuple(symbol.strip().upper() for symbol in missing_state_symbols)
    if not validated_symbols:
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires at least one missing-state symbol",
            code="posterior_missing_nucleotide_symbols_empty",
        )
    if any(len(symbol) != 1 for symbol in validated_symbols):
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires single-character missing-state symbols",
            code="posterior_missing_nucleotide_symbol_length_invalid",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    if set(validated_symbols) & set(_DNA_STATES):
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires missing-state symbols distinct from A, C, G, and T",
            code="posterior_missing_nucleotide_symbol_conflicts_with_dna_state",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    if len(set(validated_symbols)) != len(validated_symbols):
        raise PhylogeneticsError(
            "posterior missing-nucleotide definition requires unique missing-state symbols",
            code="posterior_missing_nucleotide_symbols_duplicated",
            details={"missing_state_symbols": list(validated_symbols)},
        )
    return validated_symbols


__all__ = [
    "PosteriorMissingContinuousTraitReport",
    "PosteriorMissingContinuousTraitTaxonSummaryRow",
    "PosteriorMissingDiscreteTraitReport",
    "PosteriorMissingDiscreteTraitStateProbabilityRow",
    "PosteriorMissingDiscreteTraitTaxonSummaryRow",
    "PosteriorMissingNucleotideDefinition",
    "PosteriorMissingNucleotideReport",
    "PosteriorMissingNucleotideSequenceRecord",
    "PosteriorMissingNucleotideSiteSummaryRow",
    "PosteriorMissingNucleotideStateProbabilityRow",
    "build_posterior_missing_nucleotide_definition",
]

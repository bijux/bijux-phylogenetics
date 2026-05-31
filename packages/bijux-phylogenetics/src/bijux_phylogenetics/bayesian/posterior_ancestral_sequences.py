from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.bayesian.partition_model_priors import (
    PartitionSubstitutionModelDefinition,
)
from bijux_phylogenetics.bayesian.partition_model_state import (
    resolve_partition_parameter_linkage_plan_from_model_parameters,
    resolve_partition_parameter_states_from_model_parameters,
)
from bijux_phylogenetics.bayesian.required_values import require_present
from bijux_phylogenetics.bayesian.state import BayesianPhylogeneticState
from bijux_phylogenetics.io.fasta import write_fasta_alignment
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    iterate_partition_sites,
    slice_partition_sequence,
    validate_locus_partitions,
)
from bijux_phylogenetics.phylo.likelihood.marginal_ancestral_probabilities import (
    evaluate_nucleotide_marginal_ancestral_probabilities,
)
from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

_DNA_STATES = ("A", "C", "G", "T")
_PARTITIONED_SUBSTITUTION_MODEL_NAME = "partitioned-dna"


@dataclass(frozen=True, slots=True)
class _PosteriorCladeMetadata:
    clade_id: str
    representative_node_id: str
    node_name: str | None
    descendant_taxa: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class _PosteriorSampleModelEvaluation:
    model_label: str
    clade_probability_rows: tuple[PosteriorAncestralStateProbabilityRow, ...]


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
    sampled_states: list[BayesianPhylogeneticState]
    | tuple[BayesianPhylogeneticState, ...],
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
    clade_metadata_by_id: dict[str, _PosteriorCladeMetadata] = {}
    clade_presence_count_by_id: dict[str, int] = {}
    state_probability_sum_by_key: dict[tuple[str, int, str], float] = {}
    sampled_substitution_models: set[str] = set()
    topology_ids = {
        sampled_state.tree.topology_id for sampled_state in validated_sampled_states
    }
    for sampled_state in validated_sampled_states:
        evaluation = _evaluate_posterior_sample_ancestral_probabilities(
            sampled_state=sampled_state,
            definition=definition,
        )
        sampled_substitution_models.add(evaluation.model_label)
        present_clade_ids: set[str] = set()
        for row in evaluation.clade_probability_rows:
            clade_metadata_by_id.setdefault(
                row.clade_id,
                _PosteriorCladeMetadata(
                    clade_id=row.clade_id,
                    representative_node_id=row.representative_node_id,
                    node_name=row.node_name,
                    descendant_taxa=tuple(row.descendant_taxa),
                ),
            )
            present_clade_ids.add(row.clade_id)
            accumulation_key = (row.clade_id, row.site_position, row.state)
            state_probability_sum_by_key[accumulation_key] = float(
                format(
                    state_probability_sum_by_key.get(accumulation_key, 0.0)
                    + row.conditional_posterior_probability,
                    ".15g",
                )
            )
        for clade_id in present_clade_ids:
            clade_presence_count_by_id[clade_id] = (
                clade_presence_count_by_id.get(clade_id, 0) + 1
            )
    sample_count = len(validated_sampled_states)
    state_probability_rows = _build_state_probability_rows(
        clade_metadata_by_id=clade_metadata_by_id,
        clade_presence_count_by_id=clade_presence_count_by_id,
        state_probability_sum_by_key=state_probability_sum_by_key,
        sample_count=sample_count,
    )
    site_summary_rows = _build_site_summary_rows(
        state_probability_rows=state_probability_rows,
        definition=definition,
    )
    sequence_records = _build_sequence_records(
        site_summary_rows=site_summary_rows,
        minimum_clade_posterior_support=definition.minimum_clade_posterior_support,
    )
    warnings = _build_summary_warnings(
        distinct_topology_count=len(topology_ids),
        sampled_substitution_models=sampled_substitution_models,
        clade_presence_count_by_id=clade_presence_count_by_id,
        sample_count=sample_count,
    )
    return PosteriorAncestralSequenceReport(
        sample_count=sample_count,
        site_count=len(definition.records[0].sequence),
        taxon_count=len(definition.records),
        distinct_topology_count=len(topology_ids),
        sampled_substitution_models=sorted(sampled_substitution_models),
        tree_uncertainty_policy=(
            "clade-marginal-posterior-aggregation-across-sampled-trees"
        ),
        posterior_probability_threshold=definition.posterior_probability_threshold,
        minimum_clade_posterior_support=definition.minimum_clade_posterior_support,
        low_confidence_state_symbol=definition.low_confidence_state_symbol,
        state_probability_rows=state_probability_rows,
        site_summary_rows=site_summary_rows,
        sequence_records=sequence_records,
        warnings=warnings,
    )


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
    locus_partition_names = tuple(
        partition.name for partition in validated_locus_partitions
    )
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
    sampled_states: list[BayesianPhylogeneticState]
    | tuple[BayesianPhylogeneticState, ...],
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


def _evaluate_posterior_sample_ancestral_probabilities(
    *,
    sampled_state: BayesianPhylogeneticState,
    definition: PosteriorAncestralSequenceDefinition,
) -> _PosteriorSampleModelEvaluation:
    if (
        definition.partition_models is not None
        and definition.locus_partitions is not None
    ):
        return _evaluate_partitioned_posterior_sample_ancestral_probabilities(
            sampled_state=sampled_state,
            definition=definition,
        )
    return _evaluate_unpartitioned_posterior_sample_ancestral_probabilities(
        sampled_state=sampled_state,
        definition=definition,
    )


def _evaluate_unpartitioned_posterior_sample_ancestral_probabilities(
    *,
    sampled_state: BayesianPhylogeneticState,
    definition: PosteriorAncestralSequenceDefinition,
) -> _PosteriorSampleModelEvaluation:
    model_name = sampled_state.model_parameters.categorical_parameters.get(
        "substitution-model"
    )
    if model_name in {None, _PARTITIONED_SUBSTITUTION_MODEL_NAME}:
        raise PhylogeneticsError(
            "posterior ancestral sequence summary requires a named non-partitioned nucleotide substitution model in every sampled state when no partition surface is configured",
            code="posterior_ancestral_sequence_model_name_missing",
            details={
                "observed_model_name": model_name,
            },
        )
    marginal_report = evaluate_nucleotide_marginal_ancestral_probabilities(
        sampled_state.tree.to_tree(),
        list(definition.records),
        model_name=model_name,
        kappa=sampled_state.model_parameters.scalar_parameters.get("kappa"),
        base_frequencies=sampled_state.model_parameters.vector_parameters.get(
            "base-frequencies"
        ),
        exchangeabilities=sampled_state.model_parameters.vector_parameters.get(
            "exchangeabilities"
        ),
    )
    return _PosteriorSampleModelEvaluation(
        model_label=model_name,
        clade_probability_rows=tuple(
            _normalize_probability_row(
                node_id=row.node_id,
                node_name=row.node_name,
                descendant_taxa=row.descendant_taxa,
                site_position=row.site_position,
                state=row.state,
                posterior_probability=row.posterior_probability,
            )
            for row in marginal_report.posterior_rows
        ),
    )


def _evaluate_partitioned_posterior_sample_ancestral_probabilities(
    *,
    sampled_state: BayesianPhylogeneticState,
    definition: PosteriorAncestralSequenceDefinition,
) -> _PosteriorSampleModelEvaluation:
    model_name = sampled_state.model_parameters.categorical_parameters.get(
        "substitution-model"
    )
    if model_name != _PARTITIONED_SUBSTITUTION_MODEL_NAME:
        raise PhylogeneticsError(
            "posterior ancestral sequence summary with partition surfaces requires every sampled state to preserve the partitioned DNA substitution-model label",
            code="posterior_ancestral_sequence_partitioned_model_name_invalid",
            details={
                "expected_model_name": _PARTITIONED_SUBSTITUTION_MODEL_NAME,
                "observed_model_name": model_name,
            },
        )
    partition_models = require_present(
        definition.partition_models,
        owner_name="posterior ancestral sequence partitioned evaluation",
        field_name="partition_models",
    )
    require_present(
        definition.locus_partitions,
        owner_name="posterior ancestral sequence partitioned evaluation",
        field_name="locus_partitions",
    )
    partition_names = tuple(
        partition_model.partition_name for partition_model in partition_models
    )
    linkage_plan = resolve_partition_parameter_linkage_plan_from_model_parameters(
        model_parameters=sampled_state.model_parameters,
        partition_names=partition_names,
    )
    partition_states = resolve_partition_parameter_states_from_model_parameters(
        model_parameters=sampled_state.model_parameters,
        partition_models=partition_models,
        linkage_plan=linkage_plan,
    )
    partition_state_by_name = {
        partition_state.partition_name: partition_state
        for partition_state in partition_states
    }
    rows: list[PosteriorAncestralStateProbabilityRow] = []
    model_fragments: list[str] = []
    for locus_partition, partition_model in zip(
        definition.locus_partitions,
        definition.partition_models,
        strict=True,
    ):
        partition_state = partition_state_by_name.get(partition_model.partition_name)
        if partition_state is None:
            raise PhylogeneticsError(
                "posterior ancestral sequence summary requires one realized parameter state for every configured partition",
                code="posterior_ancestral_sequence_partition_state_missing",
                details={"partition_name": partition_model.partition_name},
            )
        site_positions = _global_partition_site_positions(locus_partition)
        partition_records = [
            AlignmentRecord(
                identifier=record.identifier,
                sequence=slice_partition_sequence(record.sequence, locus_partition),
            )
            for record in definition.records
        ]
        marginal_report = evaluate_nucleotide_marginal_ancestral_probabilities(
            sampled_state.tree.to_tree(),
            partition_records,
            model_name=partition_model.base_model_name,
            kappa=partition_state.kappa,
            base_frequencies=(
                dict(partition_state.base_frequencies)
                if partition_state.base_frequencies is not None
                else None
            ),
            exchangeabilities=(
                dict(partition_state.exchangeabilities)
                if partition_state.exchangeabilities is not None
                else None
            ),
        )
        for row in marginal_report.posterior_rows:
            rows.append(
                _normalize_probability_row(
                    node_id=row.node_id,
                    node_name=row.node_name,
                    descendant_taxa=row.descendant_taxa,
                    site_position=site_positions[row.site_position - 1],
                    state=row.state,
                    posterior_probability=row.posterior_probability,
                )
            )
        model_fragments.append(
            f"{partition_model.partition_name}={partition_model.model_name}"
        )
    return _PosteriorSampleModelEvaluation(
        model_label=(
            f"{_PARTITIONED_SUBSTITUTION_MODEL_NAME}[{','.join(model_fragments)}]"
        ),
        clade_probability_rows=tuple(rows),
    )


def _normalize_probability_row(
    *,
    node_id: str,
    node_name: str | None,
    descendant_taxa: list[str],
    site_position: int,
    state: str,
    posterior_probability: float,
) -> PosteriorAncestralStateProbabilityRow:
    canonical_descendant_taxa = tuple(sorted(descendant_taxa))
    return PosteriorAncestralStateProbabilityRow(
        clade_id=canonical_clade_id(frozenset(canonical_descendant_taxa)),
        representative_node_id=node_id,
        node_name=node_name,
        descendant_taxa=list(canonical_descendant_taxa),
        site_position=site_position,
        state=state,
        clade_posterior_probability=0.0,
        conditional_posterior_probability=float(format(posterior_probability, ".15g")),
        marginal_posterior_probability=0.0,
        supporting_sample_count=0,
        total_sample_count=0,
    )


def _build_state_probability_rows(
    *,
    clade_metadata_by_id: dict[str, _PosteriorCladeMetadata],
    clade_presence_count_by_id: dict[str, int],
    state_probability_sum_by_key: dict[tuple[str, int, str], float],
    sample_count: int,
) -> list[PosteriorAncestralStateProbabilityRow]:
    rows: list[PosteriorAncestralStateProbabilityRow] = []
    for clade_id, metadata in sorted(
        clade_metadata_by_id.items(),
        key=lambda item: (len(item[1].descendant_taxa), item[1].descendant_taxa),
    ):
        supporting_sample_count = clade_presence_count_by_id[clade_id]
        clade_probability = float(
            format(supporting_sample_count / sample_count, ".15g")
        )
        for site_position in range(
            1, _infer_site_count(state_probability_sum_by_key, clade_id) + 1
        ):
            for state in _DNA_STATES:
                posterior_sum = state_probability_sum_by_key.get(
                    (clade_id, site_position, state),
                    0.0,
                )
                conditional_probability = (
                    float(format(posterior_sum / supporting_sample_count, ".15g"))
                    if supporting_sample_count > 0
                    else 0.0
                )
                marginal_probability = float(
                    format(posterior_sum / sample_count, ".15g")
                )
                rows.append(
                    PosteriorAncestralStateProbabilityRow(
                        clade_id=clade_id,
                        representative_node_id=metadata.representative_node_id,
                        node_name=metadata.node_name,
                        descendant_taxa=list(metadata.descendant_taxa),
                        site_position=site_position,
                        state=state,
                        clade_posterior_probability=clade_probability,
                        conditional_posterior_probability=conditional_probability,
                        marginal_posterior_probability=marginal_probability,
                        supporting_sample_count=supporting_sample_count,
                        total_sample_count=sample_count,
                    )
                )
    return rows


def _build_site_summary_rows(
    *,
    state_probability_rows: list[PosteriorAncestralStateProbabilityRow],
    definition: PosteriorAncestralSequenceDefinition,
) -> list[PosteriorAncestralSiteSummaryRow]:
    grouped_rows: dict[
        tuple[str, int], dict[str, PosteriorAncestralStateProbabilityRow]
    ] = {}
    for row in state_probability_rows:
        grouped_rows.setdefault((row.clade_id, row.site_position), {})[row.state] = row
    summary_rows: list[PosteriorAncestralSiteSummaryRow] = []
    for key in sorted(grouped_rows):
        state_rows = grouped_rows[key]
        ordered_state_rows = [state_rows[state] for state in _DNA_STATES]
        best_row = max(
            ordered_state_rows,
            key=lambda row: row.marginal_posterior_probability,
        )
        max_conditional_probability = max(
            row.conditional_posterior_probability for row in ordered_state_rows
        )
        max_marginal_probability = best_row.marginal_posterior_probability
        low_confidence = (
            max_marginal_probability < definition.posterior_probability_threshold
        )
        summary_rows.append(
            PosteriorAncestralSiteSummaryRow(
                clade_id=best_row.clade_id,
                representative_node_id=best_row.representative_node_id,
                node_name=best_row.node_name,
                descendant_taxa=list(best_row.descendant_taxa),
                site_position=best_row.site_position,
                consensus_state=best_row.state,
                exported_state=(
                    definition.low_confidence_state_symbol
                    if low_confidence
                    else best_row.state
                ),
                clade_posterior_probability=best_row.clade_posterior_probability,
                max_conditional_posterior_probability=float(
                    format(max_conditional_probability, ".15g")
                ),
                max_marginal_posterior_probability=float(
                    format(max_marginal_probability, ".15g")
                ),
                low_confidence=low_confidence,
                posterior_probability_a=state_rows["A"].marginal_posterior_probability,
                posterior_probability_c=state_rows["C"].marginal_posterior_probability,
                posterior_probability_g=state_rows["G"].marginal_posterior_probability,
                posterior_probability_t=state_rows["T"].marginal_posterior_probability,
            )
        )
    return summary_rows


def _build_sequence_records(
    *,
    site_summary_rows: list[PosteriorAncestralSiteSummaryRow],
    minimum_clade_posterior_support: float,
) -> list[PosteriorAncestralSequenceRecord]:
    sequence_by_clade_id: dict[str, list[str]] = {}
    row_by_clade_id: dict[str, PosteriorAncestralSiteSummaryRow] = {}
    for row in site_summary_rows:
        if row.clade_posterior_probability < minimum_clade_posterior_support:
            continue
        sequence_by_clade_id.setdefault(row.clade_id, []).append(row.exported_state)
        row_by_clade_id.setdefault(row.clade_id, row)
    return [
        PosteriorAncestralSequenceRecord(
            clade_id=clade_id,
            representative_node_id=row_by_clade_id[clade_id].representative_node_id,
            node_name=row_by_clade_id[clade_id].node_name,
            descendant_taxa=list(row_by_clade_id[clade_id].descendant_taxa),
            clade_posterior_probability=row_by_clade_id[
                clade_id
            ].clade_posterior_probability,
            sequence="".join(sequence_by_clade_id[clade_id]),
        )
        for clade_id in sorted(
            sequence_by_clade_id,
            key=lambda current_clade_id: (
                len(row_by_clade_id[current_clade_id].descendant_taxa),
                row_by_clade_id[current_clade_id].descendant_taxa,
            ),
        )
    ]


def _build_summary_warnings(
    *,
    distinct_topology_count: int,
    sampled_substitution_models: set[str],
    clade_presence_count_by_id: dict[str, int],
    sample_count: int,
) -> list[str]:
    warnings: list[str] = []
    if distinct_topology_count > 1:
        warnings.append(
            "posterior ancestral sequence summary aggregated comparable clades across multiple sampled topologies"
        )
    if len(sampled_substitution_models) > 1:
        warnings.append(
            "posterior ancestral sequence summary aggregated ancestral probabilities across multiple sampled substitution-model surfaces"
        )
    if any(count < sample_count for count in clade_presence_count_by_id.values()):
        warnings.append(
            "one or more ancestral clades are absent from some posterior samples, so marginal probabilities include topology uncertainty"
        )
    return warnings


def _global_partition_site_positions(partition: LocusPartition) -> tuple[int, ...]:
    return tuple(
        site
        for segment in partition.segments
        for site in iterate_partition_sites(segment)
    )


def _infer_site_count(
    state_probability_sum_by_key: dict[tuple[str, int, str], float],
    clade_id: str,
) -> int:
    return max(
        site_position
        for current_clade_id, site_position, _state in state_probability_sum_by_key
        if current_clade_id == clade_id
    )


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

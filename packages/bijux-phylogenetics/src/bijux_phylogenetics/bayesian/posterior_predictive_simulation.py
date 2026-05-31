from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
import random
from statistics import mean, median

import numpy

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.runtime.errors import PhylogeneticsError
from bijux_phylogenetics.simulation.contracts import (
    SimulatedContinuousTrait,
    SimulatedDiscreteTrait,
)

from .brownian_continuous_trait import BrownianContinuousTraitRunReport
from .discrete_trait_mk import DiscreteTraitMkRunReport
from .fixed_topology_dna import FixedTopologyDnaRunReport
from .fixed_topology_partitioned_dna import FixedTopologyPartitionedDnaRunReport
from .ornstein_uhlenbeck_continuous_trait import (
    OrnsteinUhlenbeckContinuousTraitRunReport,
)

POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES = ("uniform-with-replacement",)


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveSimulationDefinition:
    """Validated posterior predictive simulation configuration."""

    replicate_count: int
    sample_selection_policy: str
    seed: int


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveObservedStatisticRow:
    """One observed-data statistic retained beside posterior predictive replicates."""

    statistic_name: str
    value: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveReplicateStatisticRow:
    """One replicate statistic generated from one sampled posterior state."""

    statistic_name: str
    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    value: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveStatisticSummaryRow:
    """One observed-vs-replicate summary for one posterior predictive statistic."""

    statistic_name: str
    observed_value: float
    replicate_count: int
    replicate_mean: float
    replicate_standard_deviation: float
    replicate_minimum: float
    replicate_median: float
    replicate_maximum: float


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveAlignmentReplicate:
    """One posterior predictive alignment replicate."""

    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    records: list[AlignmentRecord]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveDiscreteTraitReplicate:
    """One posterior predictive discrete-trait replicate."""

    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    traits: list[SimulatedDiscreteTrait]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveContinuousTraitReplicate:
    """One posterior predictive continuous-trait replicate."""

    replicate_index: int
    posterior_sample_index: int
    posterior_iteration_index: int
    traits: list[SimulatedContinuousTrait]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveAlignmentSimulationReport:
    """Posterior predictive alignment simulation output for one Bayesian run."""

    definition: PosteriorPredictiveSimulationDefinition
    model_name: str
    taxon_count: int
    alignment_length: int
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow]
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow]
    statistic_summary_rows: list[PosteriorPredictiveStatisticSummaryRow]
    replicates: list[PosteriorPredictiveAlignmentReplicate]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveDiscreteTraitSimulationReport:
    """Posterior predictive discrete-trait simulation output for one Bayesian run."""

    definition: PosteriorPredictiveSimulationDefinition
    model_name: str
    taxon_count: int
    state_order: list[str]
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow]
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow]
    statistic_summary_rows: list[PosteriorPredictiveStatisticSummaryRow]
    replicates: list[PosteriorPredictiveDiscreteTraitReplicate]


@dataclass(frozen=True, slots=True)
class PosteriorPredictiveContinuousTraitSimulationReport:
    """Posterior predictive continuous-trait simulation output for one Bayesian run."""

    definition: PosteriorPredictiveSimulationDefinition
    model_name: str
    taxon_count: int
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow]
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow]
    statistic_summary_rows: list[PosteriorPredictiveStatisticSummaryRow]
    replicates: list[PosteriorPredictiveContinuousTraitReplicate]


def build_posterior_predictive_simulation_definition(
    *,
    replicate_count: int,
    sample_selection_policy: str = "uniform-with-replacement",
    seed: int = 0,
) -> PosteriorPredictiveSimulationDefinition:
    """Build one validated posterior predictive simulation configuration."""
    validated_policy = sample_selection_policy.strip()
    if validated_policy not in POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES:
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires one supported posterior-sample selection policy",
            code="posterior_predictive_simulation_policy_invalid",
            details={
                "sample_selection_policy": sample_selection_policy,
                "supported_policies": list(
                    POSTERIOR_PREDICTIVE_SAMPLE_SELECTION_POLICIES
                ),
            },
        )
    if not isinstance(replicate_count, int):
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires replicate_count to be one integer",
            code="posterior_predictive_simulation_replicate_count_type_invalid",
        )
    if replicate_count <= 0:
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires replicate_count to be positive",
            code="posterior_predictive_simulation_replicate_count_invalid",
            details={"replicate_count": replicate_count},
        )
    if not isinstance(seed, int):
        raise PhylogeneticsError(
            "posterior predictive simulation definition requires seed to be one integer",
            code="posterior_predictive_simulation_seed_type_invalid",
        )
    return PosteriorPredictiveSimulationDefinition(
        replicate_count=replicate_count,
        sample_selection_policy=validated_policy,
        seed=seed,
    )


def _build_posterior_predictive_statistic_summary_rows(
    *,
    observed_statistic_rows: list[PosteriorPredictiveObservedStatisticRow],
    replicate_statistic_rows: list[PosteriorPredictiveReplicateStatisticRow],
) -> list[PosteriorPredictiveStatisticSummaryRow]:
    replicate_values_by_statistic: dict[str, list[float]] = {}
    for row in replicate_statistic_rows:
        replicate_values_by_statistic.setdefault(row.statistic_name, []).append(row.value)
    summary_rows: list[PosteriorPredictiveStatisticSummaryRow] = []
    for observed_row in observed_statistic_rows:
        replicate_values = replicate_values_by_statistic.get(observed_row.statistic_name)
        if not replicate_values:
            raise PhylogeneticsError(
                "posterior predictive statistic summaries require at least one replicate value for every observed statistic",
                code="posterior_predictive_statistic_replicates_missing",
                details={"statistic_name": observed_row.statistic_name},
            )
        summary_rows.append(
            PosteriorPredictiveStatisticSummaryRow(
                statistic_name=observed_row.statistic_name,
                observed_value=observed_row.value,
                replicate_count=len(replicate_values),
                replicate_mean=_round_float(mean(replicate_values)),
                replicate_standard_deviation=_round_float(
                    _sample_standard_deviation(replicate_values)
                ),
                replicate_minimum=_round_float(min(replicate_values)),
                replicate_median=_round_float(median(replicate_values)),
                replicate_maximum=_round_float(max(replicate_values)),
            )
        )
    return summary_rows


def _choose_posterior_sample_indices(
    *,
    sample_count: int,
    definition: PosteriorPredictiveSimulationDefinition,
    rng: random.Random,
) -> list[int]:
    if sample_count <= 0:
        raise PhylogeneticsError(
            "posterior predictive simulation requires at least one posterior sample",
            code="posterior_predictive_sample_count_invalid",
            details={"sample_count": sample_count},
        )
    if definition.sample_selection_policy != "uniform-with-replacement":
        raise PhylogeneticsError(
            "posterior predictive simulation encountered one unsupported posterior-sample selection policy",
            code="posterior_predictive_simulation_policy_unsupported",
            details={"sample_selection_policy": definition.sample_selection_policy},
        )
    return [rng.randrange(sample_count) for _ in range(definition.replicate_count)]


def _validate_alignment_records(
    *,
    records: Sequence[AlignmentRecord],
    owner_name: str,
) -> list[AlignmentRecord]:
    validated_records = list(records)
    if not validated_records:
        raise PhylogeneticsError(
            f"{owner_name} requires at least one alignment record",
            code="posterior_predictive_alignment_records_empty",
        )
    if any(not isinstance(record, AlignmentRecord) for record in validated_records):
        raise PhylogeneticsError(
            f"{owner_name} requires every record to be one AlignmentRecord",
            code="posterior_predictive_alignment_record_type_invalid",
        )
    sequence_lengths = {len(record.sequence) for record in validated_records}
    if len(sequence_lengths) != 1:
        raise PhylogeneticsError(
            f"{owner_name} requires one aligned sequence matrix with equal sequence lengths",
            code="posterior_predictive_alignment_length_invalid",
        )
    return validated_records


def _sample_standard_deviation(values: Sequence[float]) -> float:
    if len(values) <= 1:
        return 0.0
    center = mean(values)
    return math.sqrt(
        math.fsum((value - center) ** 2 for value in values) / (len(values) - 1)
    )


def _round_float(value: float) -> float:
    return float(format(value, ".15g"))

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class BenchmarkObservation:
    label: str
    item_count: int
    runtime_seconds: float
    peak_memory_bytes: int


@dataclass(slots=True)
class TreeValidationBenchmarkReport:
    replicates: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class TreeComparisonBenchmarkReport:
    replicates: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class AlignmentDiagnosticsBenchmarkReport:
    replicates: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class AlignmentSiteBenchmarkReport:
    replicates: int
    sequence_count: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class TreeSetConsensusBenchmarkReport:
    replicates: int
    tip_count: int
    observations: list[BenchmarkObservation]


@dataclass(slots=True)
class LargeTreeScalingWorkflowBenchmark:
    workflow: str
    scaling_axis: str
    observations: list[BenchmarkObservation]
    notes: list[str]


@dataclass(slots=True)
class LargeTreeScalingBenchmarkReport:
    replicates: int
    tip_counts: list[int]
    workflows: list[LargeTreeScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class LargeAlignmentScalingObservation:
    label: str
    sequence_count: int
    alignment_length: int
    aligned_site_count: int
    runtime_seconds: float
    peak_memory_bytes: int


@dataclass(slots=True)
class LargeAlignmentScalingWorkflowBenchmark:
    workflow: str
    scaling_axis: str
    observations: list[LargeAlignmentScalingObservation]
    notes: list[str]


@dataclass(slots=True)
class LargeAlignmentScalingBenchmarkReport:
    replicates: int
    sequence_counts: list[int]
    alignment_lengths: list[int]
    workflows: list[LargeAlignmentScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class LargeTreeSetScalingObservation:
    label: str
    tree_count: int
    tip_count: int
    pair_count: int
    runtime_seconds: float
    peak_memory_bytes: int


@dataclass(slots=True)
class LargeTreeSetScalingWorkflowBenchmark:
    workflow: str
    scaling_axis: str
    observations: list[LargeTreeSetScalingObservation]
    notes: list[str]


@dataclass(slots=True)
class LargeTreeSetScalingBenchmarkReport:
    replicates: int
    tree_counts: list[int]
    tip_counts: list[int]
    workflows: list[LargeTreeSetScalingWorkflowBenchmark]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class WorkflowPracticalLimitEntry:
    workflow: str
    evidence_source: str
    tested_taxon_limit: int | None
    tested_site_limit: int | None
    tested_tree_limit: int | None
    tested_posterior_size: int | None
    max_runtime_seconds: float
    max_peak_memory_bytes: int
    memory_observation_kind: str | None
    notes: list[str]


@dataclass(slots=True)
class WorkflowPracticalLimitReport:
    replicates: int
    stress_tiers: list[str]
    entries: list[WorkflowPracticalLimitEntry]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class LargeDatasetStressObservation:
    workload: str
    tier: str
    timeout_seconds: float
    input_size_bytes: int
    sequence_count: int | None
    alignment_length: int | None
    tree_count: int | None
    taxon_count: int | None
    locus_count: int | None
    runtime_seconds: float
    peak_memory_bytes: int
    memory_observation_kind: str
    output_row_count: int
    notes: list[str]


@dataclass(slots=True)
class LargeDatasetStressSuiteReport:
    tier: str
    observations: list[LargeDatasetStressObservation]
    limitations: list[str]


@dataclass(frozen=True, slots=True)
class _StressObservationPayload:
    workload: str
    input_size_bytes: int
    sequence_count: int | None
    alignment_length: int | None
    tree_count: int | None
    taxon_count: int | None
    locus_count: int | None
    output_row_count: int
    notes: list[str]


@dataclass(frozen=True, slots=True)
class _StressTierConfig:
    tier: str
    timeout_seconds: float
    alignment_sequence_count: int
    alignment_length: int
    supermatrix_taxon_count: int
    supermatrix_locus_lengths: tuple[int, ...]
    tree_set_tree_count: int
    tree_set_tip_count: int
    comparative_taxon_count: int
    table_tip_count: int

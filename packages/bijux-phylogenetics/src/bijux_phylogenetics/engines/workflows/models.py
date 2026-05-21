from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.engines.artifacts.iqtree import IqtreeModelSelectionSummary
from bijux_phylogenetics.engines.artifacts.support import (
    BootstrapSupportSummaryReport,
    FastTreeSupportSummaryReport,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from bijux_phylogenetics.engines.common import EngineRunReport
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    CodingSequenceExclusion,
)
from bijux_phylogenetics.phylo.alignment.partitions import PartitionSummaryReport


@dataclass(slots=True)
class EngineWorkflowReport:
    workflow: str
    engine_name: str
    input_paths: list[Path]
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object] = field(default_factory=dict)
    selected_model: str | None = None
    log_likelihood: float | None = None
    iqtree_summary: IqtreeWorkflowSummary | None = None
    model_selection_summary: IqtreeModelSelectionSummary | None = None
    bootstrap_support_summary: BootstrapSupportSummaryReport | None = None
    fasttree_support_summary: FastTreeSupportSummaryReport | None = None
    sh_alrt_support_summary: ShAlrtSupportSummaryReport | None = None
    weak_backbone_report: WeakBackboneReport | None = None
    trimming_summary: AlignmentTrimmingSummary | None = None
    resumed: bool = False
    notes: list[str] = field(default_factory=list)


@dataclass(slots=True)
class IqtreeSupportValue:
    node: str
    descendant_taxa: list[str]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(slots=True)
class IqtreeWorkflowSummary:
    selected_model: str | None
    log_likelihood: float | None
    support_value_count: int
    minimum_support: float | None
    maximum_support: float | None
    support_values: list[IqtreeSupportValue] = field(default_factory=list)


@dataclass(slots=True)
class AlignmentTrimmingSummary:
    mode: str
    gap_threshold: float | None
    input_alignment_length: int
    trimmed_alignment_length: int
    retained_site_count: int
    removed_site_count: int
    retained_site_fraction: float
    removed_site_fraction: float
    input_gap_fraction: float
    trimmed_gap_fraction: float
    input_gap_percentage: float
    trimmed_gap_percentage: float


@dataclass(slots=True)
class ExternalTreeComparisonReport:
    fast_tree_path: Path
    ml_tree_path: Path
    comparison_report: ComparisonReportBuildResult


@dataclass(slots=True)
class CodonAwareAlignmentWorkflowReport:
    workflow: str
    engine_name: str
    input_path: Path
    output_paths: dict[str, Path]
    run: EngineRunReport
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    sequence_type: AlignmentAlphabet
    genetic_code_id: int
    genetic_code_name: str
    input_sequence_count: int
    accepted_sequence_count: int
    invalid_codon_sequence_count: int
    excluded_sequences: list[CodingSequenceExclusion]
    terminal_stop_sequence_count: int
    notes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    resumed: bool = False


@dataclass(frozen=True, slots=True)
class PreparedIqtreePartitions:
    command_args: list[str]
    summary: PartitionSummaryReport
    output_paths: dict[str, Path]
    notes: list[str]
    mixed_data_types: bool

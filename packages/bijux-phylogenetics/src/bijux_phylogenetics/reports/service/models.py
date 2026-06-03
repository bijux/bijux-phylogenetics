from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.dataset import (
    DatasetAuditReport,
    DatasetCrosswalkReport,
    DatasetExclusionTable,
    DatasetReadinessSummary,
)
from bijux_phylogenetics.datasets.study_inputs import (
    MetadataJoinRow,
    TraitMissingValueReport,
)
from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
)
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAmbiguousColumnReport,
    AlignmentForensicReport,
    AlignmentLinkageReport,
    AlignmentLowInformationReport,
    AlignmentQualityReport,
    AlignmentSummary,
    CodingAlignmentDiagnostics,
    DuplicateSequencePolicyReport,
    SequenceIdentityMatrix,
    SequenceQualityRankingReport,
)
from bijux_phylogenetics.phylo.taxa import (
    TaxonAuditReport,
    TaxonStabilityReport,
    TaxonWorkflowLossReport,
)
from bijux_phylogenetics.trees import (
    TreeSetProcessingSummary,
    TreeSetWorkflowBudgetReport,
)
from bijux_phylogenetics.validation import (
    CoreWorkflowValidationReport,
    LevelOneReleaseGateReport,
    ProductionScaleReadinessReport,
    ReleaseTruthReport,
)


@dataclass(slots=True)
class TableLinkageReport:
    tree_path: Path
    table_path: Path
    tree_taxa: int
    table_rows: int
    linked_taxa: int
    missing_from_table: list[str]
    extra_table_entries: list[str]
    index_column: str
    annotated_taxa: list[str]
    joined_rows: list[MetadataJoinRow]


@dataclass(frozen=True, slots=True)
class ReportInputLedgerEntry:
    path: Path
    role: str
    checksum: str
    taxa_count: int
    usage: list[str]


@dataclass(slots=True)
class ReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    validation: TreeValidationReport
    inspection: TreeInspectionReport
    forensic: TreeForensicReport
    metadata_linkage: TableLinkageReport | None
    traits_linkage: TableLinkageReport | None
    trait_missing_values: TraitMissingValueReport | None
    alignment: AlignmentSummary | None
    alignment_quality: AlignmentQualityReport | None
    alignment_forensic: AlignmentForensicReport | None
    alignment_low_information: AlignmentLowInformationReport | None
    alignment_duplicate_policy: DuplicateSequencePolicyReport | None
    alignment_ambiguous_columns: AlignmentAmbiguousColumnReport | None
    alignment_sequence_ranking: SequenceQualityRankingReport | None
    alignment_coding: CodingAlignmentDiagnostics | None
    alignment_identity_matrix: SequenceIdentityMatrix | None
    alignment_linkage: AlignmentLinkageReport | None
    dataset_readiness: DatasetReadinessSummary | None
    dataset_audit: DatasetAuditReport | None
    input_ledger: list[ReportInputLedgerEntry]
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class DistanceReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    source_path: Path
    source_kind: str
    method_limitations: list[str]
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TreeUncertaintyReportBuildResult:
    output_path: Path
    artifact_root: Path
    artifact_manifest_path: Path
    report_kind: str
    title: str
    source_path: Path
    tree_count: int
    rooted_topology_count: int
    processing: TreeSetProcessingSummary
    budget_report: TreeSetWorkflowBudgetReport
    methods_summary_path: Path
    methods_summary_warning_count: int
    limitations: list[str]
    linked_artifact_count: int
    html_size_bytes: int
    linked_artifact_bytes: int
    total_output_bytes: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TreeSetComparisonReportBuildResult:
    output_path: Path
    report_kind: str
    title: str
    left_path: Path
    right_path: Path
    shared_rooted_topology_count: int
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class AlignmentReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    alignment: AlignmentSummary
    alignment_quality: AlignmentQualityReport
    alignment_forensic: AlignmentForensicReport
    alignment_low_information: AlignmentLowInformationReport
    alignment_duplicate_policy: DuplicateSequencePolicyReport
    alignment_ambiguous_columns: AlignmentAmbiguousColumnReport
    alignment_sequence_ranking: SequenceQualityRankingReport
    alignment_coding: CodingAlignmentDiagnostics | None
    alignment_identity_matrix: SequenceIdentityMatrix
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class TaxonReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    tree_path: Path
    taxon_audit: TaxonAuditReport
    taxon_crosswalk: DatasetCrosswalkReport | None
    taxon_exclusions: DatasetExclusionTable | None
    taxon_workflow_loss: TaxonWorkflowLossReport | None
    taxon_stability: TaxonStabilityReport | None
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class WorkflowValidationReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    validation: CoreWorkflowValidationReport
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class ReleaseGateReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    release_gate: LevelOneReleaseGateReport
    machine_manifest: dict[str, object]


@dataclass(slots=True)
class ReleaseTruthReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    release_truth: ReleaseTruthReport
    machine_manifest: dict[str, object]


@dataclass(frozen=True, slots=True)
class ProductionScaleReadinessReportBuildResult:
    output_path: Path
    machine_manifest_path: Path
    report_kind: str
    title: str
    production_scale_readiness: ProductionScaleReadinessReport
    machine_manifest: dict[str, object]


__all__ = [
    "AlignmentReportBuildResult",
    "DistanceReportBuildResult",
    "ProductionScaleReadinessReportBuildResult",
    "ReleaseGateReportBuildResult",
    "ReleaseTruthReportBuildResult",
    "ReportBuildResult",
    "ReportInputLedgerEntry",
    "TableLinkageReport",
    "TaxonReportBuildResult",
    "TreeSetComparisonReportBuildResult",
    "TreeUncertaintyReportBuildResult",
    "WorkflowValidationReportBuildResult",
]

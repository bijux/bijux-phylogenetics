from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import (
    TraitLinkageReport,
    TraitValidationReport,
)
from bijux_phylogenetics.diagnostics.validation import (
    TreeInspectionReport,
    TreeValidationReport,
)
from bijux_phylogenetics.phylo.alignment import (
    AlignmentRecord,
    CodingSequencePreparationReport,
    FastaInputValidationReport,
    FastaRepairReport,
    SequenceCompositionOutlier,
    SequenceLengthOutlier,
)


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelDataset:
    """Packaged stress panel with explicit sequence, tree, and trait defects."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    raw_alignment_path: Path
    raw_tree_path: Path
    raw_traits_path: Path
    raw_sequence_input_path: Path
    raw_coding_sequences_path: Path
    raw_trait_mismatch_path: Path
    reference_output_root: Path
    taxon_count: int
    raw_trait_row_count: int
    required_traits: tuple[str, ...]
    sequence_type: str
    source_summary: str


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelExportResult:
    """Materialized copy of the packaged stress panel."""

    output_root: Path
    readme_path: Path
    raw_alignment_path: Path
    raw_tree_path: Path
    raw_traits_path: Path
    raw_sequence_input_path: Path
    raw_coding_sequences_path: Path
    raw_trait_mismatch_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class TraitDuplicateResolution:
    """Deterministic duplicate-taxon choice for one raw trait surface."""

    taxon: str
    occurrence_count: int
    selected_row_number: int
    selected_non_missing_field_count: int
    discarded_row_numbers: list[int]
    selected_reason: str


@dataclass(slots=True)
class TraitMissingObservation:
    """One missing trait observation in the raw traits table."""

    taxon: str
    row_number: int
    trait: str
    required_for_analysis: bool
    action: str


@dataclass(slots=True)
class DataQualityRepairAction:
    """One explicit action taken to move from raw inputs to the cleaned subset."""

    action_kind: str
    affected_taxa: list[str]
    affected_nodes: list[str]
    reason: str
    result: str


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelWorkflowReport:
    """Audit and cleanup report for the packaged stress panel."""

    dataset: CatarrhineDataQualityStressPanelDataset
    raw_sequence_input_validation: FastaInputValidationReport
    raw_sequence_input_repair: FastaRepairReport
    raw_sequence_length_outliers: list[SequenceLengthOutlier]
    repaired_sequence_input_validation: FastaInputValidationReport
    coding_sequence_preparation: CodingSequencePreparationReport
    raw_trait_mismatch_linkage: TraitLinkageReport
    raw_trait_mismatch_error: str | None
    raw_alignment_validation: FastaInputValidationReport
    sequence_outliers: list[SequenceCompositionOutlier]
    raw_tree_inspection: TreeInspectionReport
    raw_tree_validation: TreeValidationReport
    trait_duplicates: list[TraitDuplicateResolution]
    missing_traits: list[TraitMissingObservation]
    cleaned_trait_validation: TraitValidationReport
    cleaned_tree_validation: TreeValidationReport
    cleaned_linkage: TraitLinkageReport
    cleaned_alignment_validation: FastaInputValidationReport
    cleaned_alignment_records: list[AlignmentRecord]
    repaired_sequence_input_path: Path
    prepared_coding_sequences_path: Path
    cleaned_tree_path: Path
    cleaned_traits_path: Path
    cleaned_alignment_path: Path
    cleaned_taxa: list[str]
    dropped_taxa: list[str]
    repair_actions: list[DataQualityRepairAction]
    repaired_branch_nodes: list[str]


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelWorkflowBundle:
    """Written workflow outputs for the packaged stress panel."""

    output_root: Path
    raw_taxon_count: int
    cleaned_taxon_count: int
    duplicate_sequence_identifier_count: int
    illegal_character_count: int
    empty_sequence_count: int
    raw_sequence_length_outlier_count: int
    duplicate_trait_taxon_count: int
    missing_trait_value_count: int
    sequence_outlier_count: int
    tree_zero_length_branch_count: int
    tree_negative_branch_count: int
    tree_long_branch_outlier_count: int
    coding_frame_error_count: int
    coding_internal_stop_count: int
    raw_trait_missing_from_traits_count: int
    raw_trait_extra_taxon_count: int
    dropped_taxon_count: int
    repaired_branch_count: int
    workflow_summary_path: Path
    raw_sequence_findings_path: Path
    raw_sequence_repair_path: Path
    repaired_sequence_input_path: Path
    repaired_sequence_validation_path: Path
    coding_sequence_exclusions_path: Path
    prepared_coding_sequences_path: Path
    raw_trait_linkage_path: Path
    trait_duplicates_path: Path
    trait_missing_values_path: Path
    sequence_outliers_path: Path
    tree_issues_path: Path
    repair_actions_path: Path
    cleaned_traits_path: Path
    cleaned_alignment_path: Path
    cleaned_tree_path: Path
    cleaned_linkage_path: Path
    cleaned_validation_path: Path


@dataclass(slots=True)
class CatarrhineDataQualityStressPanelDemoResult:
    """Dataset export plus written workflow bundle for the public demo."""

    output_root: Path
    dataset: CatarrhineDataQualityStressPanelDataset
    dataset_export: CatarrhineDataQualityStressPanelExportResult
    workflow_bundle: CatarrhineDataQualityStressPanelWorkflowBundle
    overview_path: Path

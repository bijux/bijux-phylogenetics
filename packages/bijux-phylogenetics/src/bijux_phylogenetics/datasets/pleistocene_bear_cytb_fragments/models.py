from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.engines.inference import FastaToTreeWorkflowReport
from bijux_phylogenetics.io.fasta import AlignmentRecord
from bijux_phylogenetics.phylo.alignment import (
    AlignmentQualityReport,
    AlignmentSummary,
    AlignmentTrimReport,
)


@dataclass(slots=True)
class PleistoceneBearCytbFragmentDataset:
    """Packaged ancient-DNA-style dataset for degraded sequence workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    sequences_path: Path
    reference_output_root: Path
    sequence_count: int
    sequence_type: str
    workflow_prefix: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    site_missingness_threshold: float
    sequence_missingness_threshold: float
    degraded_sequence_ids: tuple[str, ...]
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class PleistoceneBearCytbFragmentDatasetExportResult:
    """Materialized copy of the packaged ancient-DNA-style dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class PleistoceneBearMissingnessEffectRow:
    """Reviewer-facing missingness comparison for one sequence across workflow stages."""

    identifier: str
    raw_sequence_length: int
    degraded_sequence: bool
    aligned_missing_fraction: float
    engine_trimmed_missing_fraction: float
    cleaned_missing_fraction: float
    removed_by_missingness_cleanup: bool


@dataclass(slots=True)
class PleistoceneBearCytbFragmentWorkflowReport:
    """One governed degraded-sequence workflow run over the packaged bear panel."""

    dataset: PleistoceneBearCytbFragmentDataset
    workflow: FastaToTreeWorkflowReport
    aligned_summary: AlignmentSummary
    trimmed_summary: AlignmentSummary
    cleaned_summary: AlignmentSummary
    aligned_quality: AlignmentQualityReport
    trimmed_quality: AlignmentQualityReport
    cleaned_quality: AlignmentQualityReport
    missingness_cleanup: AlignmentTrimReport
    cleaned_records: list[AlignmentRecord]
    missingness_rows: list[PleistoceneBearMissingnessEffectRow]


@dataclass(slots=True)
class PleistoceneBearCytbFragmentWorkflowBundle:
    """Written degraded-sequence workflow outputs for the packaged bear panel."""

    output_root: Path
    selected_model: str
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    removed_column_count: int
    removed_sequence_count: int
    cleaned_missing_data_fraction: float
    summary_path: Path
    missingness_effects_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    cleaned_alignment_path: Path
    tree_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path


@dataclass(slots=True)
class PleistoceneBearCytbFragmentDemoResult:
    """Dataset export plus workflow outputs for the public degraded-sequence demo."""

    output_root: Path
    dataset: PleistoceneBearCytbFragmentDataset
    dataset_export: PleistoceneBearCytbFragmentDatasetExportResult
    workflow_bundle: PleistoceneBearCytbFragmentWorkflowBundle
    overview_path: Path

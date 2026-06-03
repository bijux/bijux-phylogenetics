from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.engines.inference import FastaToTreeWorkflowReport


@dataclass(slots=True)
class InfluenzaAHAReferenceDataset:
    """Packaged viral FASTA dataset for sequence-to-tree workflow review."""

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
    source_accessions: tuple[str, ...]
    source_summary: str


@dataclass(slots=True)
class InfluenzaAHAReferenceDatasetExportResult:
    """Materialized copy of the packaged viral dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class InfluenzaAHAReferenceWorkflowReport:
    """One governed FASTA-to-tree workflow run over the packaged viral dataset."""

    dataset: InfluenzaAHAReferenceDataset
    workflow: FastaToTreeWorkflowReport


@dataclass(slots=True)
class InfluenzaAHAReferenceWorkflowBundle:
    """Written sequence-to-tree outputs for the packaged viral dataset."""

    output_root: Path
    selected_model: str
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    summary_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path


@dataclass(slots=True)
class InfluenzaAHAReferenceDemoResult:
    """Dataset export plus workflow outputs for the public viral demo."""

    output_root: Path
    dataset: InfluenzaAHAReferenceDataset
    dataset_export: InfluenzaAHAReferenceDatasetExportResult
    workflow_bundle: InfluenzaAHAReferenceWorkflowBundle
    overview_path: Path

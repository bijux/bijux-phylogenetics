from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.engines.inference import FastaToTreeWorkflowReport


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkDataset:
    """Packaged protein FASTA benchmark for end-to-end phylogeny review."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    sequences_path: Path
    reference_output_root: Path
    sequence_count: int
    sequence_type: str
    workflow_prefix: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    source_reference: str
    source_transformation: str
    minimum_sequence_length: int
    maximum_sequence_length: int
    source_summary: str


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkExportResult:
    """Materialized copy of the packaged protein benchmark dataset."""

    output_root: Path
    readme_path: Path
    sequences_path: Path
    expected_output_root: Path


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkWorkflowReport:
    """One governed protein FASTA-to-tree workflow run over the packaged panel."""

    dataset: GnathostomeOrthologProteinBenchmarkDataset
    workflow: FastaToTreeWorkflowReport


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkWorkflowBundle:
    """Written alignment, inference, support, and assumption outputs."""

    output_root: Path
    selected_model: str
    sequence_count: int
    alignment_length: int
    trimmed_alignment_length: int
    minimum_support: float | None
    maximum_support: float | None
    median_support: float | None
    weakly_supported_clade_count: int
    summary_path: Path
    assumptions_path: Path
    alignment_path: Path
    trimmed_alignment_path: Path
    tree_path: Path
    model_table_path: Path
    support_table_path: Path
    log_path: Path
    manifest_path: Path
    engine_artifact_root: Path


@dataclass(slots=True)
class GnathostomeOrthologProteinBenchmarkDemoResult:
    """Dataset export plus rerun workflow outputs for the public protein demo."""

    output_root: Path
    dataset: GnathostomeOrthologProteinBenchmarkDataset
    dataset_export: GnathostomeOrthologProteinBenchmarkExportResult
    workflow_bundle: GnathostomeOrthologProteinBenchmarkWorkflowBundle
    overview_path: Path

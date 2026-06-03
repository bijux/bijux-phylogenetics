from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivityVariant:
    """One declared method combination in the packaged sensitivity matrix."""

    variant_id: str
    label: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


@dataclass(slots=True)
class RabiesMethodSensitivityPanelDataset:
    """Packaged rabies dataset for method-sensitivity workflow review."""

    dataset_id: str
    label: str
    dataset_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    reference_output_root: Path
    taxon_count: int
    sequence_type: str
    workflow_prefix: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    parallel_workers: int
    source_accessions: tuple[str, ...]
    variants: tuple[RabiesMethodSensitivityVariant, ...]
    source_summary: str


@dataclass(slots=True)
class RabiesMethodSensitivityPanelExportResult:
    """Materialized copy of the packaged rabies method-sensitivity dataset."""

    output_root: Path
    readme_path: Path
    config_path: Path
    sequences_path: Path
    metadata_path: Path
    expected_output_root: Path

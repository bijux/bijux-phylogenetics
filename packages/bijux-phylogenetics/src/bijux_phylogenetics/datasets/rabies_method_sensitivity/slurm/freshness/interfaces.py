from __future__ import annotations

from pathlib import Path
from typing import Protocol


class VariantLike(Protocol):
    variant_id: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


class DatasetLike(Protocol):
    dataset_id: str
    workflow_prefix: str
    sequence_type: str
    outgroup_taxa: tuple[str, ...]
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    sequences_path: Path
    metadata_path: Path
    variants: tuple[VariantLike, ...]

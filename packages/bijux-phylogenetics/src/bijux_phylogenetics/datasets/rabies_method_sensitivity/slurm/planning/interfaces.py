from __future__ import annotations

from pathlib import Path
from typing import Protocol


class VariantConfigLike(Protocol):
    variant_id: str
    alignment_mode: str
    trimming_mode: str


class VariantRunLike(Protocol):
    config: VariantConfigLike
    alignment_length: int
    trimmed_alignment_length: int


class TaskRecordLike(Protocol):
    variant_id: str
    output_root: Path


class DatasetLike(Protocol):
    dataset_id: str
    workflow_prefix: str
    taxon_count: int


class WorkflowReportLike(Protocol):
    dataset: DatasetLike
    task_records: tuple[TaskRecordLike, ...]
    variant_runs: tuple[VariantRunLike, ...]
    iqtree_threads: int
    bootstrap_replicates: int

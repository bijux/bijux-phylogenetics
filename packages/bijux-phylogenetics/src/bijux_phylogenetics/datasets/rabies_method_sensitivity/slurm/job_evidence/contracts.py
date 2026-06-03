from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Protocol


class TaskRecordLike(Protocol):
    variant_id: str
    label: str
    execution_mode: str
    status: str
    log_path: Path
    error_code: str | None
    error_message: str | None


class VariantConfigLike(Protocol):
    variant_id: str
    label: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float


class EngineVersionLike(Protocol):
    text: str


class EngineRunLike(Protocol):
    version: EngineVersionLike
    command: list[str]
    warning_lines: list[str]
    runtime_seconds: float
    started_at_utc: str
    ended_at_utc: str


class TrimmingSummaryLike(Protocol):
    retained_site_fraction: float
    removed_site_fraction: float


class EngineWorkflowLike(Protocol):
    manifest_path: Path
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    run: EngineRunLike
    trimming_summary: TrimmingSummaryLike | None


class ConclusionSummaryLike(Protocol):
    stable_clade_count: int
    unstable_clade_count: int
    engine_specific_clade_count: int
    serious_conflict_count: int


class InferenceComparisonLike(Protocol):
    manifest_path: Path
    step_manifests: dict[str, Path]
    commands: dict[str, list[str]]
    engine_versions: dict[str, str]
    runtime_seconds: float
    selected_model: str
    conclusion_summary: ConclusionSummaryLike
    warnings: list[str]


class RootingLike(Protocol):
    requested_taxa: tuple[str, ...]
    matched_taxa: tuple[str, ...]
    outgroup_monophyletic: bool | None
    rooted_outgroup_taxa: tuple[str, ...]
    warnings: list[str]


class RootedComparisonLike(Protocol):
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    same_taxa_different_rooting: bool


class VariantRunLike(Protocol):
    config: VariantConfigLike
    alignment_workflow: EngineWorkflowLike
    trimming_workflow: EngineWorkflowLike
    inference_comparison: InferenceComparisonLike
    fasttree_rooting: RootingLike
    iqtree_rooting: RootingLike
    rooted_engine_comparison: RootedComparisonLike
    alignment_length: int
    trimmed_alignment_length: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmJobEvidenceRow:
    """One independently reviewable provenance package for one Slurm job."""

    partition_id: str
    array_index: int
    variant_id: str
    label: str
    execution_mode: str
    status: str
    script_path: str
    output_root: str
    evidence_directory: str
    evidence_json_path: str
    evidence_html_path: str
    task_log_copy_path: str
    alignment_manifest_path: str
    trimming_manifest_path: str
    inference_manifest_path: str
    model_selection_manifest_path: str
    iqtree_support_manifest_path: str
    fasttree_manifest_path: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float
    selected_model: str
    alignment_length: int
    trimmed_alignment_length: int
    total_runtime_seconds: float
    alignment_runtime_seconds: float
    trimming_runtime_seconds: float
    inference_runtime_seconds: float
    model_selection_runtime_seconds: float
    iqtree_support_runtime_seconds: float
    fasttree_runtime_seconds: float
    rooted_engine_rf_distance: int
    rooted_engine_normalized_rf: float
    rooted_same_taxa_different_rooting: bool
    serious_conflict_count: int
    stable_clade_count: int
    unstable_clade_count: int
    engine_specific_clade_count: int
    warning_count: int
    output_file_count: int
    output_byte_count: int
    artifact_file_count: int


@dataclass(frozen=True, slots=True)
class RabiesMethodSensitivitySlurmJobEvidenceReport:
    """One reviewer-facing per-job evidence surface over the full workflow."""

    dataset_id: str
    workflow_prefix: str
    execution_mode: str
    parallel_workers: int
    bundle_root: Path
    evidence_root: Path
    index_path: Path
    summary_path: Path
    job_count: int
    completed_job_count: int
    failed_job_count: int
    total_runtime_seconds: float
    total_output_file_count: int
    total_output_byte_count: int
    total_artifact_file_count: int
    jobs: tuple[RabiesMethodSensitivitySlurmJobEvidenceRow, ...]

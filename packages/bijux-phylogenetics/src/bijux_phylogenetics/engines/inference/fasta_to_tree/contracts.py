from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.evidence.provenance.method_tiers import MethodTierAssessment
from bijux_phylogenetics.phylo.alignment import (
    AlignmentAlphabet,
    FastaInputValidationReport,
    FastaRepairReport,
)

from ...artifacts.support import BootstrapSupportSummaryReport
from ...validation import ModelSelectionValidationReport
from ...workflows.models import EngineWorkflowReport


@dataclass(frozen=True, slots=True)
class FastaToTreeModelRow:
    """One reviewer-facing record describing the selected substitution model."""

    workflow: str
    engine_name: str
    sequence_type: AlignmentAlphabet
    selected_model: str
    report_selected_model: str | None
    artifact_selected_model: str | None
    model_consistent: bool
    alignment_path: Path
    trimmed_alignment_path: Path
    manifest_path: Path


@dataclass(frozen=True, slots=True)
class FastaToTreeSupportRow:
    """One reviewer-facing branch-support record from the final tree."""

    node: str
    descendant_taxa: tuple[str, ...]
    support: float
    support_fraction: float
    is_backbone: bool


@dataclass(frozen=True, slots=True)
class FastaToTreeStageFingerprint:
    """One deterministic fingerprint record for one workflow stage."""

    stage: str
    fingerprint: str
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    config: dict[str, object]
    engine_versions: dict[str, str]
    upstream_fingerprints: dict[str, str]
    resumed: bool


@dataclass(slots=True)
class FastaToTreeWorkflowReport:
    """End-to-end result for one raw-FASTA-to-tree workflow run."""

    workflow: str
    input_path: Path
    prepared_input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet
    selected_model: str
    alignment_mode: str
    trimming_mode: str
    trim_gap_threshold: float
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    started_at_utc: str
    ended_at_utc: str
    runtime_seconds: float
    engine_artifact_dir: Path
    manifest_path: Path
    run_manifest_path: Path
    output_paths: dict[str, Path]
    config: dict[str, object]
    commands: dict[str, list[str]]
    engine_versions: dict[str, str]
    step_manifests: dict[str, Path]
    stage_fingerprints: dict[str, FastaToTreeStageFingerprint]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    input_validation: FastaInputValidationReport
    repaired_input_validation: FastaInputValidationReport | None
    input_repair: FastaRepairReport | None
    alignment_workflow: EngineWorkflowReport
    trimming_workflow: EngineWorkflowReport
    model_selection_workflow: EngineWorkflowReport
    maximum_likelihood_workflow: EngineWorkflowReport
    bootstrap_workflow: EngineWorkflowReport
    model_validation: ModelSelectionValidationReport
    support_summary: BootstrapSupportSummaryReport
    model_rows: list[FastaToTreeModelRow]
    support_rows: list[FastaToTreeSupportRow]
    method_tier: MethodTierAssessment
    warnings: list[str]
    notes: list[str]

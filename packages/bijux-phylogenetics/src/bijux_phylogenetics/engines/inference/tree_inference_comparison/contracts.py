from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.compare.influence import TaxonInfluenceReport
from bijux_phylogenetics.compare.presentation import ComparisonReportBuildResult
from bijux_phylogenetics.engines.validation import InferenceTreeComparisonReport
from bijux_phylogenetics.phylo.alignment import AlignmentAlphabet

from ...workflows.models import EngineWorkflowReport


@dataclass(frozen=True, slots=True)
class InferenceComparisonSharedCladeRow:
    """One shared clade across FastTree and IQ-TREE support workflows."""

    split_id: str
    fasttree_support: float | None
    fasttree_support_fraction: float | None
    fasttree_support_label_kind: str
    iqtree_support: float | None
    iqtree_support_fraction: float | None
    iqtree_support_label_kind: str
    support_fraction_delta: float | None
    support_disagreement: bool


@dataclass(frozen=True, slots=True)
class InferenceComparisonConflictRow:
    """One clade-level conflict record across the two inference engines."""

    split_id: str
    conflict_kind: str
    fasttree_present: bool
    iqtree_present: bool
    fasttree_support: float | None
    fasttree_support_fraction: float | None
    iqtree_support: float | None
    iqtree_support_fraction: float | None
    detail: str


@dataclass(frozen=True, slots=True)
class InferenceComparisonWeightedConflictRow:
    """One ranked conflict or disagreement with explicit support weighting."""

    split_id: str
    comparison_status: str
    conflict_kind: str
    severity_class: str
    fasttree_support_fraction: float | None
    iqtree_support_fraction: float | None
    support_fraction_delta: float | None
    strongest_support_fraction: float | None
    support_weight: float | None
    serious_conflict: bool
    detail: str


@dataclass(frozen=True, slots=True)
class InferenceComparisonConclusionRow:
    """One reviewer-facing clade conclusion from an engine comparison."""

    split_id: str
    conclusion_class: str
    evidence_class: str
    comparison_status: str
    fasttree_present: bool
    iqtree_present: bool
    fasttree_support_fraction: float | None
    iqtree_support_fraction: float | None
    support_fraction_delta: float | None
    serious_conflict: bool
    detail: str


@dataclass(frozen=True, slots=True)
class InferenceComparisonConclusionSummary:
    """Compact summary of whether compared biological conclusions are stable."""

    shared_taxa_count: int
    robinson_foulds_distance: int
    normalized_robinson_foulds: float
    branch_score_distance: float | None
    stable_clade_count: int
    unstable_clade_count: int
    engine_specific_clade_count: int
    support_weighted_conflict_count: int
    low_support_disagreement_count: int
    moderate_support_disagreement_count: int
    high_support_conflict_count: int
    high_support_disagreement_count: int
    serious_conflict_count: int
    top_conflict_driver_taxa: list[str]


@dataclass(slots=True)
class InferenceComparisonWorkflowReport:
    """End-to-end result for one engine-comparison workflow on one alignment."""

    workflow: str
    input_path: Path
    out_dir: Path
    prefix: str
    sequence_type: AlignmentAlphabet | None
    selected_model: str
    iqtree_seed: int
    iqtree_threads: int
    bootstrap_replicates: int
    timeout_seconds: float | None
    started_at_utc: str
    ended_at_utc: str
    runtime_seconds: float
    engine_artifact_dir: Path
    manifest_path: Path
    output_paths: dict[str, Path]
    step_manifests: dict[str, Path]
    config: dict[str, object]
    commands: dict[str, list[str]]
    engine_versions: dict[str, str]
    input_checksums: dict[str, str]
    output_checksums: dict[str, str]
    model_selection_workflow: EngineWorkflowReport
    iqtree_support_workflow: EngineWorkflowReport
    fasttree_workflow: EngineWorkflowReport
    engine_comparison: InferenceTreeComparisonReport
    comparison_report: ComparisonReportBuildResult
    shared_clade_rows: list[InferenceComparisonSharedCladeRow]
    conflicting_clade_rows: list[InferenceComparisonConflictRow]
    weighted_conflict_rows: list[InferenceComparisonWeightedConflictRow]
    conclusion_rows: list[InferenceComparisonConclusionRow]
    conclusion_summary: InferenceComparisonConclusionSummary
    taxon_influence_report: TaxonInfluenceReport | None
    warnings: list[str]
    notes: list[str]

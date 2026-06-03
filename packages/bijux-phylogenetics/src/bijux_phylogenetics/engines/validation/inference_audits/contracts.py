from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class InferenceReadinessDecision:
    workflow: str
    ready: bool
    blockers: list[str]
    warnings: list[str]


@dataclass(slots=True)
class InferenceReadinessAuditReport:
    alignment_path: Path
    sequence_count: int
    alignment_length: int | None
    inferred_alphabet: str
    overall_decision: str
    recommended_workflow: str
    decisions: list[InferenceReadinessDecision]
    warnings: list[str]


@dataclass(slots=True)
class ModelSelectionValidationReport:
    manifest_path: Path
    manifest_selected_model: str | None
    manifest_selected_criterion: str | None
    report_selected_model: str | None
    report_selected_criterion: str | None
    artifact_selected_model: str | None
    candidate_model_count: int
    best_model_aic: str | None
    best_model_aicc: str | None
    best_model_bic: str | None
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class MLTreeTaxonValidationReport:
    manifest_path: Path
    expected_taxa: list[str]
    observed_taxa: list[str]
    missing_taxa: list[str]
    unexpected_taxa: list[str]
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class MetadataClusterObservation:
    group: str
    tree_taxa: list[str]
    monophyletic: bool | None
    status: str
    note: str


@dataclass(slots=True)
class MetadataClusteringReport:
    tree_path: Path
    metadata_path: Path
    taxon_column: str
    group_column: str
    group_count: int
    monophyletic_group_count: int
    split_group_count: int
    observations: list[MetadataClusterObservation]


@dataclass(slots=True)
class InferenceFailureTaxonomyReport:
    workflow: str
    failure_category: str
    failure_reason: str
    scientific_explanation: str
    likely_causes: list[str]
    actionable_fixes: list[str]
    evidence: dict[str, object]
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class BootstrapTreeSetValidationReport:
    tree_set_path: Path
    tree_count: int
    expected_taxa: list[str]
    valid: bool
    issues: list[str]


@dataclass(slots=True)
class InferenceTreeComparisonReport:
    comparison_kind: str
    left_manifest_path: Path
    right_manifest_path: Path
    left_label: str
    right_label: str
    left_tree_path: Path
    right_tree_path: Path
    left_engine_name: str
    right_engine_name: str
    left_selected_model: str | None
    right_selected_model: str | None
    topology: object
    support: object
    branch_lengths: object
    warnings: list[str]


@dataclass(slots=True)
class InferenceOutputConsistencyReport:
    manifest_path: Path
    workflow: str
    failure_category: str
    current_output_checksum_match: bool
    valid: bool
    issues: list[str]

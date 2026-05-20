from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.phylo.taxa import TaxonIdentityAudit


@dataclass(slots=True)
class TreeIntegrityIssue:
    code: str
    message: str
    severity: str
    affected_nodes: list[str]


@dataclass(slots=True)
class TreeFinding:
    code: str
    message: str
    severity: str
    affected_taxa: list[str]
    affected_nodes: list[str]


@dataclass(slots=True)
class RootStateConfidenceReport:
    classification: str
    rationale: list[str]
    suspicious_placement: bool
    suspicious_reasons: list[str]


@dataclass(slots=True)
class BranchLengthContextAssessment:
    context: str
    allowed: bool
    blocked_by: list[str]
    warnings: list[str]


@dataclass(slots=True)
class BranchLengthRepairSuggestion:
    issue_code: str
    summary: str
    blocked_analyses: list[str]
    suggested_action: str


@dataclass(slots=True)
class InternalNodeLabelConflict:
    node_id: str
    label: str
    conflict_type: str
    detail: str


@dataclass(slots=True)
class StableNodeIdentity:
    node_id: str
    descendant_taxa: list[str]


@dataclass(slots=True)
class UnsafeExternalLabel:
    raw_label: str
    normalized_label: str
    engines: list[str]
    reasons: list[str]


@dataclass(slots=True)
class TreeForensicReport:
    path: Path
    source_format: str
    syntax_valid: bool
    biologically_safe: bool
    validity_decision: str
    integrity_issues: list[TreeIntegrityIssue]
    findings: list[TreeFinding]
    root_state_confidence: RootStateConfidenceReport
    branch_length_contexts: list[BranchLengthContextAssessment]
    branch_length_repair_suggestions: list[BranchLengthRepairSuggestion]
    internal_label_conflicts: list[InternalNodeLabelConflict]
    stable_node_identities: list[StableNodeIdentity]
    unsafe_external_labels: list[UnsafeExternalLabel]
    taxon_identity_audit: TaxonIdentityAudit
    safe_for_topology_comparison: bool
    safe_for_time_tree_analysis: bool
    safe_for_comparative_methods: bool
    safe_for_visualization: bool
    safe_for_publication: bool
    warnings: list[str]


@dataclass(slots=True)
class TreeValidationReport:
    path: Path
    source_format: str
    tip_count: int
    internal_node_count: int
    rooted: bool
    syntax_valid: bool
    biologically_safe: bool
    validity_decision: str
    has_complete_branch_lengths: bool
    branch_length_status: str
    total_branch_length: float
    zero_length_branches: int
    negative_branch_lengths: int
    polytomy_count: int
    polytomy_nodes: list[str]
    missing_internal_branch_nodes: list[str]
    missing_terminal_branch_taxa: list[str]
    singleton_internal_nodes: list[str]
    missing_taxa: int
    duplicate_taxa: list[str]
    ultrametric: bool | None
    integrity_issues: list[TreeIntegrityIssue]
    warning_details: list[TreeFinding]
    root_state_confidence: RootStateConfidenceReport
    branch_length_contexts: list[BranchLengthContextAssessment]
    branch_length_repair_suggestions: list[BranchLengthRepairSuggestion]
    internal_label_conflicts: list[InternalNodeLabelConflict]
    stable_node_identities: list[StableNodeIdentity]
    unsafe_external_labels: list[UnsafeExternalLabel]
    taxon_identity_audit: TaxonIdentityAudit
    warnings: list[str]


@dataclass(slots=True)
class BranchLengthSummary:
    count: int
    minimum: float
    maximum: float
    mean: float
    median: float
    first_quartile: float
    third_quartile: float


@dataclass(slots=True)
class TreeQualityWarning:
    code: str
    message: str
    penalty: float
    affected_taxa: list[str]
    affected_nodes: list[str]


@dataclass(slots=True)
class InternalNodeChildCount:
    node: str
    child_count: int


@dataclass(slots=True)
class BranchLengthOutlier:
    node: str
    branch_length: float
    branch_type: str


@dataclass(slots=True)
class InternalLabelInterpretation:
    node: str
    node_id: str
    label: str
    interpretation: str
    numeric_value: float | None


@dataclass(slots=True)
class TreeInspectionReport:
    path: Path
    source_format: str
    tip_count: int
    node_count: int
    internal_node_count: int
    edge_count: int
    clade_count: int
    rooted: bool
    root_state_confidence: RootStateConfidenceReport
    is_binary: bool
    internal_child_counts: list[InternalNodeChildCount]
    singleton_internal_nodes: list[str]
    polytomy_count: int
    polytomy_nodes: list[str]
    has_branch_lengths: bool
    branch_length_status: str
    missing_internal_branch_nodes: list[str]
    missing_terminal_branch_taxa: list[str]
    is_ultrametric: bool | None
    total_branch_length: float
    branch_length_summary: BranchLengthSummary | None
    tree_diameter: float | None
    zero_length_branch_count: int
    min_root_to_tip: float | None
    max_root_to_tip: float | None
    max_depth: int
    mean_depth: float
    colless_imbalance_index: float | None
    normalized_colless_imbalance: float | None
    sackin_imbalance_index: int
    unusually_imbalanced: bool | None
    long_branch_taxa: list[str]
    long_branch_outliers: list[BranchLengthOutlier]
    short_branch_outliers: list[BranchLengthOutlier]
    suspicious_support_value_ranges: list[str]
    mixed_support_scales: bool
    likely_support_labels: list[InternalLabelInterpretation]
    likely_named_internal_labels: list[InternalLabelInterpretation]
    internal_label_conflicts: list[InternalNodeLabelConflict]
    stable_node_identities: list[StableNodeIdentity]
    unsafe_external_labels: list[UnsafeExternalLabel]
    taxon_identity_audit: TaxonIdentityAudit
    star_like: bool
    comb_like: bool
    tree_quality_score: float
    tree_quality_warnings: list[TreeQualityWarning]
    imbalance_summary: str
    cherry_count: int
    taxa: list[str]
    warnings: list[str]


@dataclass(slots=True)
class TreeDiagnosticReport:
    path: Path
    inspection: TreeInspectionReport
    validation: TreeValidationReport
    forensic: TreeForensicReport

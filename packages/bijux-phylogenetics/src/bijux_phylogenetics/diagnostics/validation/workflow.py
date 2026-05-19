from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.runtime.errors import (
    DuplicateTaxonError,
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    UnnamedTipError,
    UnrootedTreeError,
)

from .branch_review import (
    _branch_length_contexts,
    _branch_length_health,
    _branch_length_repair_suggestions,
    _branch_length_status,
    _missing_internal_branch_nodes,
    _missing_terminal_branch_taxa,
    _singleton_internal_nodes,
)
from .inspection import inspect_tree_path
from .models import (
    InternalNodeLabelConflict,
    RootStateConfidenceReport,
    TreeDiagnosticReport,
    TreeFinding,
    TreeForensicReport,
    TreeInspectionReport,
    TreeIntegrityIssue,
    TreeValidationReport,
    UnsafeExternalLabel,
)
from .structure import (
    _duplicate_taxa,
    _integrity_issues,
    _load_tree,
    _polytomy_nodes,
    _ultrametric,
)


def _findings_from_reports(
    integrity_issues: list[TreeIntegrityIssue],
    inspection: TreeInspectionReport,
    duplicate_taxa: list[str],
    missing_taxa: int,
    root_state: RootStateConfidenceReport,
    internal_label_conflicts: list[InternalNodeLabelConflict],
    unsafe_labels: list[UnsafeExternalLabel],
) -> list[TreeFinding]:
    findings: list[TreeFinding] = [
        TreeFinding(
            code=issue.code,
            message=issue.message,
            severity=issue.severity,
            affected_taxa=[],
            affected_nodes=issue.affected_nodes,
        )
        for issue in integrity_issues
    ]
    if duplicate_taxa:
        findings.append(
            TreeFinding(
                code="duplicate_taxa",
                message="tree contains duplicate tip labels and is biologically unsafe until taxa are disambiguated",
                severity="blocker",
                affected_taxa=duplicate_taxa,
                affected_nodes=[],
            )
        )
    if missing_taxa:
        findings.append(
            TreeFinding(
                code="unnamed_tips",
                message="tree contains unnamed tips and cannot be safely reconciled to external data",
                severity="blocker",
                affected_taxa=[],
                affected_nodes=[],
            )
        )
    for warning in inspection.tree_quality_warnings:
        severity = "warning"
        if warning.code in {
            "negative_branch_lengths",
            "missing_branch_lengths",
            "partial_branch_lengths",
        }:
            severity = "blocker"
        findings.append(
            TreeFinding(
                code=warning.code,
                message=warning.message,
                severity=severity,
                affected_taxa=warning.affected_taxa,
                affected_nodes=warning.affected_nodes,
            )
        )
    if root_state.suspicious_placement:
        findings.append(
            TreeFinding(
                code="suspicious_root_placement",
                message="root placement appears biologically suspicious",
                severity="warning",
                affected_taxa=[],
                affected_nodes=["<root>"],
            )
        )
    for conflict in internal_label_conflicts:
        findings.append(
            TreeFinding(
                code=conflict.conflict_type,
                message=conflict.detail,
                severity="warning",
                affected_taxa=[],
                affected_nodes=[]
                if conflict.node_id == "<tree>"
                else [conflict.node_id],
            )
        )
    if unsafe_labels:
        findings.append(
            TreeFinding(
                code="unsafe_external_labels",
                message="one or more taxon labels are unsafe across common phylogenetics engines or shell workflows",
                severity="warning",
                affected_taxa=[row.raw_label for row in unsafe_labels],
                affected_nodes=[],
            )
        )
    if (
        inspection.taxon_identity_audit.whitespace_variants
        or inspection.taxon_identity_audit.underscore_space_collisions
        or inspection.taxon_identity_audit.case_collisions
        or inspection.taxon_identity_audit.suspicious_near_duplicates
    ):
        findings.append(
            TreeFinding(
                code="taxon_identity_conflicts",
                message="one or more taxon labels are suspiciously similar and may not represent distinct biological identities",
                severity="warning",
                affected_taxa=sorted(
                    {
                        pair.left_label
                        for pair in (
                            inspection.taxon_identity_audit.whitespace_variants
                            + inspection.taxon_identity_audit.underscore_space_collisions
                            + inspection.taxon_identity_audit.case_collisions
                            + inspection.taxon_identity_audit.suspicious_near_duplicates
                        )
                    }
                    | {
                        pair.right_label
                        for pair in (
                            inspection.taxon_identity_audit.whitespace_variants
                            + inspection.taxon_identity_audit.underscore_space_collisions
                            + inspection.taxon_identity_audit.case_collisions
                            + inspection.taxon_identity_audit.suspicious_near_duplicates
                        )
                    }
                ),
                affected_nodes=[],
            )
        )
    severity_order = {"fatal": 0, "blocker": 1, "warning": 2, "info": 3}
    return sorted(
        findings,
        key=lambda row: (severity_order.get(row.severity, 99), row.code, row.message),
    )


def _validity_decision(findings: list[TreeFinding]) -> tuple[bool, bool, str]:
    has_fatal = any(finding.severity == "fatal" for finding in findings)
    has_blocker = any(finding.severity == "blocker" for finding in findings)
    has_warning = any(finding.severity == "warning" for finding in findings)
    syntax_valid = not has_fatal
    biologically_safe = syntax_valid and not has_blocker
    if not syntax_valid or has_blocker:
        return syntax_valid, biologically_safe, "invalid"
    if has_warning:
        return syntax_valid, biologically_safe, "valid_with_warnings"
    return syntax_valid, biologically_safe, "valid"


def validate_tree_path(
    path: Path,
    *,
    source_format: str | None = None,
    allow_duplicates: bool = False,
    strict: bool = False,
    allow_negative_branch_lengths: bool = False,
    require_rooted: bool = False,
    require_ultrametric: bool = False,
) -> TreeValidationReport:
    """Validate a tree file and produce a diagnostic report."""
    tree = _load_tree(path, source_format=source_format)
    inspection = inspect_tree_path(path, source_format=source_format)
    rooted = inspection.rooted
    has_complete, zero_count, negative_count = _branch_length_health(tree)
    branch_length_status = _branch_length_status(tree)
    missing_internal_branch_nodes = _missing_internal_branch_nodes(tree)
    missing_terminal_branch_taxa = _missing_terminal_branch_taxa(tree)
    singleton_internal_nodes = _singleton_internal_nodes(tree)
    missing_taxa, duplicate_taxa = _duplicate_taxa(tree)
    integrity_issues = _integrity_issues(tree)
    if duplicate_taxa and not allow_duplicates:
        raise DuplicateTaxonError(
            f"duplicate tip labels found: {', '.join(duplicate_taxa)}"
        )
    if missing_taxa and strict:
        raise UnnamedTipError(f"tree contains {missing_taxa} unnamed tip labels")
    if negative_count and not allow_negative_branch_lengths:
        raise InvalidBranchLengthError(
            f"tree contains {negative_count} negative branch lengths"
        )
    ultrametric = _ultrametric(tree)
    if require_rooted and not rooted:
        raise UnrootedTreeError(f"tree is not rooted: {path}")
    if require_ultrametric and ultrametric is not True:
        raise NonUltrametricTreeError(
            f"tree is not ultrametric within default validation tolerance: {path}"
        )
    polytomy_nodes = _polytomy_nodes(tree)
    warnings: list[str] = []
    if duplicate_taxa:
        warnings.append("tree contains duplicate tip labels")
    if missing_taxa:
        warnings.append("tree contains unnamed tips")
    if negative_count:
        warnings.append("tree contains negative branch lengths")
    if zero_count:
        warnings.append("tree contains zero-length branches")
    if missing_internal_branch_nodes:
        warnings.append("tree contains internal branches without lengths")
    if missing_terminal_branch_taxa:
        warnings.append("tree contains terminal branches without lengths")
    if singleton_internal_nodes:
        warnings.append("tree contains singleton internal nodes")
    if polytomy_nodes:
        warnings.append("tree contains one or more polytomies")
    if inspection.root_state_confidence.suspicious_placement:
        warnings.append("tree root placement appears biologically suspicious")
    if inspection.unsafe_external_labels:
        warnings.append(
            "tree contains taxon labels unsafe across common external engines"
        )
    if (
        inspection.taxon_identity_audit.whitespace_variants
        or inspection.taxon_identity_audit.underscore_space_collisions
        or inspection.taxon_identity_audit.case_collisions
        or inspection.taxon_identity_audit.suspicious_near_duplicates
    ):
        warnings.append("tree contains potentially ambiguous taxon identity variants")
    findings = _findings_from_reports(
        integrity_issues,
        inspection,
        duplicate_taxa,
        missing_taxa,
        inspection.root_state_confidence,
        inspection.internal_label_conflicts,
        inspection.unsafe_external_labels,
    )
    syntax_valid, biologically_safe, validity_decision = _validity_decision(findings)
    return TreeValidationReport(
        path=path,
        source_format=tree.source_format,
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        rooted=rooted,
        syntax_valid=syntax_valid,
        biologically_safe=biologically_safe,
        validity_decision=validity_decision,
        has_complete_branch_lengths=has_complete,
        branch_length_status=branch_length_status,
        total_branch_length=tree.total_branch_length(),
        zero_length_branches=zero_count,
        negative_branch_lengths=negative_count,
        polytomy_count=len(polytomy_nodes),
        polytomy_nodes=polytomy_nodes,
        missing_internal_branch_nodes=missing_internal_branch_nodes,
        missing_terminal_branch_taxa=missing_terminal_branch_taxa,
        singleton_internal_nodes=singleton_internal_nodes,
        missing_taxa=missing_taxa,
        duplicate_taxa=duplicate_taxa,
        ultrametric=ultrametric,
        integrity_issues=integrity_issues,
        warning_details=findings,
        root_state_confidence=inspection.root_state_confidence,
        branch_length_contexts=_branch_length_contexts(inspection),
        branch_length_repair_suggestions=_branch_length_repair_suggestions(inspection),
        internal_label_conflicts=inspection.internal_label_conflicts,
        stable_node_identities=inspection.stable_node_identities,
        unsafe_external_labels=inspection.unsafe_external_labels,
        taxon_identity_audit=inspection.taxon_identity_audit,
        warnings=warnings,
    )


def forensic_tree_path(
    path: Path, *, source_format: str | None = None
) -> TreeForensicReport:
    """Build a reviewer-facing forensic summary of whether a tree is safe for downstream use."""
    inspection = inspect_tree_path(path, source_format=source_format)
    validation = validate_tree_path(path, source_format=source_format)
    context_lookup = {
        context.context: context for context in validation.branch_length_contexts
    }
    safe_for_topology_comparison = (
        validation.syntax_valid
        and not validation.duplicate_taxa
        and validation.missing_taxa == 0
    )
    safe_for_time_tree_analysis = (
        context_lookup["time_tree"].allowed and validation.biologically_safe
    )
    safe_for_comparative_methods = (
        context_lookup["comparative_methods"].allowed and validation.biologically_safe
    )
    safe_for_visualization = validation.syntax_valid
    safe_for_publication = (
        validation.biologically_safe and not inspection.internal_label_conflicts
    )
    warnings = list(dict.fromkeys([*validation.warnings, *inspection.warnings]))
    return TreeForensicReport(
        path=path,
        source_format=validation.source_format,
        syntax_valid=validation.syntax_valid,
        biologically_safe=validation.biologically_safe,
        validity_decision=validation.validity_decision,
        integrity_issues=validation.integrity_issues,
        findings=validation.warning_details,
        root_state_confidence=validation.root_state_confidence,
        branch_length_contexts=validation.branch_length_contexts,
        branch_length_repair_suggestions=validation.branch_length_repair_suggestions,
        internal_label_conflicts=validation.internal_label_conflicts,
        stable_node_identities=validation.stable_node_identities,
        unsafe_external_labels=validation.unsafe_external_labels,
        taxon_identity_audit=validation.taxon_identity_audit,
        safe_for_topology_comparison=safe_for_topology_comparison,
        safe_for_time_tree_analysis=safe_for_time_tree_analysis,
        safe_for_comparative_methods=safe_for_comparative_methods,
        safe_for_visualization=safe_for_visualization,
        safe_for_publication=safe_for_publication,
        warnings=warnings,
    )


def diagnose_tree_path(
    path: Path, *, source_format: str | None = None
) -> TreeDiagnosticReport:
    """Return a combined inspection and validation report for one tree."""
    return TreeDiagnosticReport(
        path=path,
        inspection=inspect_tree_path(path, source_format=source_format),
        validation=validate_tree_path(path, source_format=source_format),
        forensic=forensic_tree_path(path, source_format=source_format),
    )

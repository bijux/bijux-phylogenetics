from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.diagnostics.validation import (
    TreeForensicReport,
    TreeInspectionReport,
    TreeValidationReport,
    inspect_tree_path,
    validate_tree_path,
)
from bijux_phylogenetics.reports.publication.tree import (
    TreeBranchStatisticsRow,
    TreeSupportRow,
    summarize_tree_branch_statistics,
    summarize_tree_support,
)
from bijux_phylogenetics.trees import (
    analyze_branch_length_distribution,
    extract_tree_clades,
)

from .columns import tree_table_columns
from .models import (
    SupplementaryTreeDiagnosticsRow,
    SupplementaryTreeDiagnosticsTableResult,
)
from .shared import stringify_list, write_dict_rows


def _serialize_tree_row(
    row: SupplementaryTreeDiagnosticsRow,
) -> dict[str, object]:
    return {
        "tree_source": row.tree_source,
        "source_format": row.source_format,
        "tip_count": row.tip_count,
        "internal_node_count": row.internal_node_count,
        "edge_count": row.edge_count,
        "clade_count": row.clade_count,
        "topology_shape": row.topology_shape,
        "is_binary": row.is_binary,
        "star_like": row.star_like,
        "comb_like": row.comb_like,
        "polytomy_count": row.polytomy_count,
        "polytomy_nodes": stringify_list(row.polytomy_nodes),
        "rooted": row.rooted,
        "root_state_classification": row.root_state_classification,
        "root_state_suspicious": row.root_state_suspicious,
        "branch_length_status": row.branch_length_status,
        "has_complete_branch_lengths": row.has_complete_branch_lengths,
        "total_branch_length": row.total_branch_length,
        "minimum_branch_length": ""
        if row.minimum_branch_length is None
        else row.minimum_branch_length,
        "maximum_branch_length": ""
        if row.maximum_branch_length is None
        else row.maximum_branch_length,
        "mean_branch_length": ""
        if row.mean_branch_length is None
        else row.mean_branch_length,
        "median_branch_length": ""
        if row.median_branch_length is None
        else row.median_branch_length,
        "positive_branch_median": ""
        if row.positive_branch_median is None
        else row.positive_branch_median,
        "missing_branch_count": row.missing_branch_count,
        "zero_length_branch_count": row.zero_length_branch_count,
        "negative_branch_count": row.negative_branch_count,
        "long_branch_outlier_count": row.long_branch_outlier_count,
        "short_branch_outlier_count": row.short_branch_outlier_count,
        "supported_branch_count": row.supported_branch_count,
        "strong_support_branch_count": row.strong_support_branch_count,
        "moderate_support_branch_count": row.moderate_support_branch_count,
        "weak_support_branch_count": row.weak_support_branch_count,
        "missing_support_branch_count": row.missing_support_branch_count,
        "support_value_range_warnings": stringify_list(
            row.support_value_range_warnings
        ),
        "ultrametric": "" if row.ultrametric is None else row.ultrametric,
        "min_root_to_tip": "" if row.min_root_to_tip is None else row.min_root_to_tip,
        "max_root_to_tip": "" if row.max_root_to_tip is None else row.max_root_to_tip,
        "tree_diameter": "" if row.tree_diameter is None else row.tree_diameter,
        "tree_quality_score": row.tree_quality_score,
        "safe_for_topology_comparison": row.safe_for_topology_comparison,
        "safe_for_time_tree_analysis": row.safe_for_time_tree_analysis,
        "safe_for_comparative_methods": row.safe_for_comparative_methods,
        "safe_for_visualization": row.safe_for_visualization,
        "safe_for_publication": row.safe_for_publication,
        "warning_count": row.warning_count,
        "warnings": stringify_list(row.warnings),
    }


def _support_counts(rows: list[TreeSupportRow]) -> dict[str, int]:
    counts = {"strong": 0, "moderate": 0, "weak": 0, "missing": 0}
    for row in rows:
        counts[row.support_class] = counts.get(row.support_class, 0) + 1
    return counts


def _topology_shape(inspection: TreeInspectionReport) -> str:
    if inspection.star_like:
        return "star"
    if inspection.comb_like:
        return "comb"
    if inspection.is_binary:
        return "binary"
    if inspection.polytomy_count:
        return "polytomy"
    return "mixed"


def _tree_warning_ledger(
    *,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
) -> list[str]:
    return sorted(
        dict.fromkeys(
            [
                *validation.warnings,
                *inspection.warnings,
                *forensic.warnings,
            ]
        )
    )


def _build_tree_forensic_review(
    *,
    tree_path: Path,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
) -> TreeForensicReport:
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
    warnings = sorted(dict.fromkeys([*validation.warnings, *inspection.warnings]))
    return TreeForensicReport(
        path=tree_path,
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


def _build_tree_row(
    *,
    tree_path: Path,
    validation: TreeValidationReport,
    inspection: TreeInspectionReport,
    forensic: TreeForensicReport,
    support_rows: list[TreeSupportRow],
    branch_stats: TreeBranchStatisticsRow,
) -> SupplementaryTreeDiagnosticsRow:
    warnings = _tree_warning_ledger(
        validation=validation,
        inspection=inspection,
        forensic=forensic,
    )
    support_counts = _support_counts(support_rows)
    supported_branch_count = sum(1 for row in support_rows if row.support is not None)
    return SupplementaryTreeDiagnosticsRow(
        tree_source=str(tree_path),
        source_format=inspection.source_format,
        tip_count=inspection.tip_count,
        internal_node_count=inspection.internal_node_count,
        edge_count=inspection.edge_count,
        clade_count=inspection.clade_count,
        topology_shape=_topology_shape(inspection),
        is_binary=inspection.is_binary,
        star_like=inspection.star_like,
        comb_like=inspection.comb_like,
        polytomy_count=inspection.polytomy_count,
        polytomy_nodes=inspection.polytomy_nodes,
        rooted=inspection.rooted,
        root_state_classification=inspection.root_state_confidence.classification,
        root_state_suspicious=inspection.root_state_confidence.suspicious_placement,
        branch_length_status=inspection.branch_length_status,
        has_complete_branch_lengths=validation.has_complete_branch_lengths,
        total_branch_length=inspection.total_branch_length,
        minimum_branch_length=branch_stats.minimum_branch_length,
        maximum_branch_length=branch_stats.maximum_branch_length,
        mean_branch_length=branch_stats.mean_branch_length,
        median_branch_length=branch_stats.median_branch_length,
        positive_branch_median=branch_stats.positive_branch_median,
        missing_branch_count=branch_stats.missing_branch_count,
        zero_length_branch_count=branch_stats.zero_length_branch_count,
        negative_branch_count=branch_stats.negative_branch_count,
        long_branch_outlier_count=branch_stats.long_outlier_count,
        short_branch_outlier_count=branch_stats.short_outlier_count,
        supported_branch_count=supported_branch_count,
        strong_support_branch_count=support_counts["strong"],
        moderate_support_branch_count=support_counts["moderate"],
        weak_support_branch_count=support_counts["weak"],
        missing_support_branch_count=support_counts["missing"],
        support_value_range_warnings=inspection.suspicious_support_value_ranges,
        ultrametric=inspection.is_ultrametric,
        min_root_to_tip=inspection.min_root_to_tip,
        max_root_to_tip=inspection.max_root_to_tip,
        tree_diameter=inspection.tree_diameter,
        tree_quality_score=inspection.tree_quality_score,
        safe_for_topology_comparison=forensic.safe_for_topology_comparison,
        safe_for_time_tree_analysis=forensic.safe_for_time_tree_analysis,
        safe_for_comparative_methods=forensic.safe_for_comparative_methods,
        safe_for_visualization=forensic.safe_for_visualization,
        safe_for_publication=forensic.safe_for_publication,
        warning_count=len(warnings),
        warnings=warnings,
    )


def _write_tree_rows(
    path: Path,
    *,
    columns: list[str],
    rows: list[SupplementaryTreeDiagnosticsRow],
) -> Path:
    return write_dict_rows(
        path,
        columns=columns,
        rows=[_serialize_tree_row(row) for row in rows],
    )


def write_supplementary_tree_diagnostics_table(
    path: Path,
    *,
    tree_path: Path,
) -> SupplementaryTreeDiagnosticsTableResult:
    """Write one supplementary tree diagnostics table with topology and warning summaries."""
    validation = validate_tree_path(
        tree_path,
        allow_duplicates=True,
        allow_negative_branch_lengths=True,
    )
    inspection = inspect_tree_path(tree_path)
    forensic = _build_tree_forensic_review(
        tree_path=tree_path,
        validation=validation,
        inspection=inspection,
    )
    clades = extract_tree_clades(tree_path)
    support_rows = summarize_tree_support(clades)
    branch_lengths = analyze_branch_length_distribution(tree_path)
    branch_stats = summarize_tree_branch_statistics(branch_lengths)
    rows = [
        _build_tree_row(
            tree_path=tree_path,
            validation=validation,
            inspection=inspection,
            forensic=forensic,
            support_rows=support_rows,
            branch_stats=branch_stats,
        )
    ]
    columns = tree_table_columns()
    _write_tree_rows(path, columns=columns, rows=rows)
    return SupplementaryTreeDiagnosticsTableResult(
        output_path=path,
        row_count=len(rows),
        columns=columns,
        rows=rows,
    )

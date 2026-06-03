"""High-level report services."""

from __future__ import annotations

from importlib import import_module

_EXPORTS = {
    "AlignmentFilteringMethodsSummaryTextResult": ".methods",
    "TreeInferenceMethodsSummaryTextResult": ".methods",
    "TreeValidationMethodsSummaryTextResult": ".methods",
    "write_alignment_filtering_methods_summary_text": ".methods",
    "write_tree_inference_methods_summary_text": ".methods",
    "write_tree_validation_methods_summary_text": ".methods",
    "AlignmentFigureAudit": ".publication",
    "AlignmentFigurePackageResult": ".publication",
    "PublicationPackageComparisonArtifactRow": ".publication",
    "PublicationPackageComparisonCheckRow": ".publication",
    "PublicationPackageComparisonResult": ".publication",
    "PublicationPackageRevalidationArtifactRow": ".publication",
    "PublicationPackageRevalidationCheckRow": ".publication",
    "PublicationPackageRevalidationResult": ".publication",
    "TreeBranchStatisticsRow": ".publication",
    "TreeReportPackageResult": ".publication",
    "TreeSupportRow": ".publication",
    "build_alignment_figure_package": ".publication",
    "build_tree_report_package": ".publication",
    "summarize_tree_branch_statistics": ".publication",
    "summarize_tree_support": ".publication",
    "write_publication_package_comparison_report": ".publication",
    "write_publication_package_revalidation_report": ".publication",
    "write_tree_branch_statistics_table": ".publication",
    "write_tree_support_table": ".publication",
    "ReviewerAuditChecklist": ".review",
    "ReviewerAuditChecklistItem": ".review",
    "ReviewerAuditChecklistWriteResult": ".review",
    "build_reviewer_audit_checklist": ".review",
    "write_reviewer_audit_checklist": ".review",
    "write_reviewer_audit_checklist_from_manifest": ".review",
    "SupplementaryAlignmentDiagnosticsRow": ".supplementary_tables",
    "SupplementaryAlignmentDiagnosticsTableResult": ".supplementary_tables",
    "SupplementaryAncestralStateRow": ".supplementary_tables",
    "SupplementaryAncestralStateTableResult": ".supplementary_tables",
    "SupplementaryBatchSummaryRow": ".supplementary_tables",
    "SupplementaryBatchSummaryTableResult": ".supplementary_tables",
    "SupplementaryCladeSupportRow": ".supplementary_tables",
    "SupplementaryCladeSupportTableResult": ".supplementary_tables",
    "SupplementaryComparativeModelRow": ".supplementary_tables",
    "SupplementaryComparativeModelTableResult": ".supplementary_tables",
    "SupplementaryDiversificationRow": ".supplementary_tables",
    "SupplementaryDiversificationTableResult": ".supplementary_tables",
    "SupplementaryModelSelectionRow": ".supplementary_tables",
    "SupplementaryModelSelectionTableResult": ".supplementary_tables",
    "SupplementaryTaxonTableResult": ".supplementary_tables",
    "SupplementaryTaxonTableRow": ".supplementary_tables",
    "SupplementaryTreeDiagnosticsRow": ".supplementary_tables",
    "SupplementaryTreeDiagnosticsTableResult": ".supplementary_tables",
    "write_supplementary_alignment_diagnostics_table": ".supplementary_tables",
    "write_supplementary_ancestral_state_table": ".supplementary_tables",
    "write_supplementary_batch_summary_table": ".supplementary_tables",
    "write_supplementary_clade_support_table": ".supplementary_tables",
    "write_supplementary_comparative_model_table": ".supplementary_tables",
    "write_supplementary_diversification_table": ".supplementary_tables",
    "write_supplementary_model_selection_table": ".supplementary_tables",
    "write_supplementary_taxon_table": ".supplementary_tables",
    "write_supplementary_tree_diagnostics_table": ".supplementary_tables",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name, __name__)
    return getattr(module, name)

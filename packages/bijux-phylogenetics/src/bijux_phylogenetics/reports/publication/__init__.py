from __future__ import annotations

from importlib import import_module

_EXPORTS = {
    "AlignmentFigureAudit": ".alignment",
    "AlignmentFigurePackageResult": ".alignment",
    "build_alignment_figure_package": ".alignment",
    "PublicationPackageComparisonArtifactRow": ".comparison",
    "PublicationPackageComparisonCheckRow": ".comparison",
    "PublicationPackageComparisonResult": ".comparison",
    "write_publication_package_comparison_report": ".comparison",
    "PublicationPackageRevalidationArtifactRow": ".revalidation",
    "PublicationPackageRevalidationCheckRow": ".revalidation",
    "PublicationPackageRevalidationResult": ".revalidation",
    "write_publication_package_revalidation_report": ".revalidation",
    "SUPPORTED_PUBLICATION_PACKAGE_KIND": ".support",
    "TreeBranchStatisticsRow": ".tree",
    "TreeReportPackageResult": ".tree",
    "TreeSupportRow": ".tree",
    "build_tree_report_package": ".tree",
    "summarize_tree_branch_statistics": ".tree",
    "summarize_tree_support": ".tree",
    "write_tree_branch_statistics_table": ".tree",
    "write_tree_support_table": ".tree",
}

__all__ = sorted(_EXPORTS)


def __getattr__(name: str):
    module_name = _EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name, __name__)
    return getattr(module, name)

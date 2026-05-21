"""Owned publication tree-figure package surface."""

from __future__ import annotations

from .audits import (
    build_collapsed_clade_summaries,
    build_surface_coverage,
    build_table_consistency,
    descendant_taxa,
    iter_collapsed_nodes,
    visible_tip_taxa,
)
from .contracts import (
    FigureAnnotationCoverage,
    FigureCaptionDraft,
    FigureCollapsedCladeSummary,
    FigureLegendAudit,
    FigureLegendEntry,
    FigureLegibilityAudit,
    FigureTableConsistencyReport,
    TreeFigureAuditReport,
    TreeFigurePackageResult,
)

__all__ = [
    "FigureAnnotationCoverage",
    "FigureCaptionDraft",
    "FigureCollapsedCladeSummary",
    "FigureLegendAudit",
    "FigureLegendEntry",
    "FigureLegibilityAudit",
    "FigureTableConsistencyReport",
    "TreeFigureAuditReport",
    "TreeFigurePackageResult",
    "build_collapsed_clade_summaries",
    "build_surface_coverage",
    "build_table_consistency",
    "descendant_taxa",
    "iter_collapsed_nodes",
    "visible_tip_taxa",
]

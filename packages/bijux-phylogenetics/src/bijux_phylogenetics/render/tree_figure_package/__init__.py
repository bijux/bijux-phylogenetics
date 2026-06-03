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
from .builder import build_tree_figure_package
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
from .legends import build_legend_audit, build_legend_entries
from .review import build_caption_draft, build_legibility_audit

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
    "build_tree_figure_package",
    "build_legend_audit",
    "build_legend_entries",
    "build_caption_draft",
    "build_collapsed_clade_summaries",
    "build_legibility_audit",
    "build_surface_coverage",
    "build_table_consistency",
    "descendant_taxa",
    "iter_collapsed_nodes",
    "visible_tip_taxa",
]

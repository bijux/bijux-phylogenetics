from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.render.tree_svg import TreeRenderResult

from ...tree_sets import TreeSetProcessingSummary, TreeSetWorkflowBudgetReport
from ..methods_text import (
    TreeSetUncertaintyMethodReport,
    TreeSetUncertaintyMethodsSummaryTextResult,
)


@dataclass(frozen=True, slots=True)
class TreeSetUncertaintyLegendEntry:
    """One explicit legend entry for a tree-set uncertainty figure package."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class TreeSetUncertaintyCaptionDraft:
    """Structured caption draft for a tree-set uncertainty figure package."""

    title: str
    lead_sentence: str
    support_sentence: str
    instability_sentence: str
    cluster_sentence: str
    limitation_sentence: str
    caption_ready: bool


@dataclass(frozen=True, slots=True)
class TreeSetUncertaintyPublicationAudit:
    """Reviewer-facing publication audit for tree-set uncertainty figures."""

    publication_ready: bool
    support_labels_validated: bool
    consensus_visible: bool
    clade_support_visible: bool
    unstable_taxa_visible: bool
    topology_clusters_visible: bool
    legend_complete: bool
    caption_ready: bool
    rendered_support_count: int
    plotted_clade_support_count: int
    plotted_unstable_taxon_count: int
    plotted_topology_cluster_count: int
    unstable_taxon_count: int
    topology_cluster_count: int
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class TreeSetUncertaintyFigurePackageResult:
    output_dir: Path
    tree_count: int
    processing: TreeSetProcessingSummary
    budget_report: TreeSetWorkflowBudgetReport
    consensus_tree_path: Path
    consensus_figure_path: Path
    clade_support_plot_path: Path
    unstable_taxa_plot_path: Path
    topology_clusters_plot_path: Path
    unstable_taxa_table_path: Path
    topology_clusters_table_path: Path
    uncertainty_conclusions_table_path: Path
    conclusion_summary_path: Path
    legend_path: Path
    caption_path: Path
    methods_summary_path: Path
    review_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    consensus_render: TreeRenderResult
    methods_report: TreeSetUncertaintyMethodReport
    methods_summary: TreeSetUncertaintyMethodsSummaryTextResult
    legend_entries: list[TreeSetUncertaintyLegendEntry]
    caption_draft: TreeSetUncertaintyCaptionDraft
    audit: TreeSetUncertaintyPublicationAudit
    machine_manifest: dict[str, object]

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.biogeography.migration import GeographicMigrationEventReport
from bijux_phylogenetics.biogeography.presentation.publication_support import (
    BiogeographyCaptionDraft,
    BiogeographyPublicationAudit,
    BiogeographyPublicationLegendEntry,
)
from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.phylogeography.geographic_map import GeographicMapReport
from bijux_phylogenetics.render.tree_svg import TreeRenderResult


@dataclass(frozen=True, slots=True)
class BiogeographyRegionCountRow:
    """One observed region count row from the analyzed taxon set."""

    region: str
    tip_taxon_count: int
    analyzed_taxon_fraction: float


@dataclass(frozen=True, slots=True)
class BiogeographyReportExclusionRow:
    """One excluded row from the full biogeography report package."""

    surface: str
    subject_id: str
    subject_kind: str
    raw_left: str
    raw_right: str
    reason: str
    note: str


@dataclass(slots=True)
class BiogeographyReportPackageResult:
    """Full biogeography review bundle for one rooted tree and one region table."""

    output_dir: Path
    report_path: Path
    tree_figure_path: Path
    map_path: Path
    legend_path: Path
    caption_path: Path
    summary_table_path: Path
    region_count_table_path: Path
    node_table_path: Path
    transition_matrix_path: Path
    event_table_path: Path
    map_marker_table_path: Path
    map_line_table_path: Path
    exclusion_table_path: Path
    manifest_path: Path
    reproducibility_manifest_path: Path
    state_report: GeographicStateModelReport
    event_report: GeographicMigrationEventReport
    map_report: GeographicMapReport
    tree_render: TreeRenderResult
    region_count_rows: list[BiogeographyRegionCountRow]
    exclusion_rows: list[BiogeographyReportExclusionRow]
    reviewer_summary: list[str]
    limitations: list[str]
    warnings: list[str]
    machine_manifest: dict[str, object]
    legend_entries: list[BiogeographyPublicationLegendEntry]
    caption_draft: BiogeographyCaptionDraft
    audit: BiogeographyPublicationAudit

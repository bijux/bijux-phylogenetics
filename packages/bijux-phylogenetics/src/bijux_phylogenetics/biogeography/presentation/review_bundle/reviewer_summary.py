from __future__ import annotations

from bijux_phylogenetics.biogeography.migration import GeographicMigrationEventReport
from bijux_phylogenetics.biogeography.presentation.publication_support import (
    BiogeographyPublicationAudit,
)
from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.phylogeography.geographic_map import GeographicMapReport

from .contracts import BiogeographyRegionCountRow


def build_reviewer_summary(
    *,
    state_report: GeographicStateModelReport,
    event_report: GeographicMigrationEventReport,
    map_report: GeographicMapReport,
    region_count_rows: list[BiogeographyRegionCountRow],
    audit: BiogeographyPublicationAudit,
) -> tuple[list[str], list[str]]:
    summary = [
        f"root region: {state_report.summary.root_region} ({state_report.summary.root_region_probability:.3f})",
        f"observed regions: {state_report.summary.observed_region_count}, analyzed taxa: {state_report.summary.analyzed_taxon_count}",
        f"changed branches: {state_report.summary.changed_branch_count}, strongly supported transitions: {state_report.summary.strongly_supported_transition_count}",
        f"migration events: {event_report.summary.event_count}, visible map lines: {map_report.summary.visible_line_count}",
        f"publication ready: {'yes' if audit.publication_ready else 'no'}",
        (
            f"most sampled region: {region_count_rows[0].region}"
            if region_count_rows
            else "most sampled region: unavailable"
        ),
    ]
    limitations: list[str] = []
    if state_report.summary.ambiguous_internal_node_count:
        limitations.append(
            "some ancestral regions remain ambiguous and should be reviewed through the node probability ledger"
        )
    if map_report.summary.excluded_record_count:
        limitations.append(
            "some map records were excluded because the region metadata or centroids could not be placed cleanly on the geographic map"
        )
    if not map_report.summary.visible_line_count:
        limitations.append(
            "the geographic map does not currently show any visible transition lines"
        )
    return summary, limitations

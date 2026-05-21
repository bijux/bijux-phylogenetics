from __future__ import annotations

from hashlib import sha256
from pathlib import Path

from bijux_phylogenetics.biogeography.migration import GeographicMigrationEventReport
from bijux_phylogenetics.biogeography.presentation.publication_support import (
    BiogeographyPublicationAudit,
)
from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport
from bijux_phylogenetics.phylogeography.geographic_map import GeographicMapReport

from .contracts import BiogeographyRegionCountRow, BiogeographyReportExclusionRow


def checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_machine_manifest(
    *,
    tree_path: Path,
    traits_path: Path,
    centroids_path: Path,
    state_report: GeographicStateModelReport,
    event_report: GeographicMigrationEventReport,
    map_report: GeographicMapReport,
    region_count_rows: list[BiogeographyRegionCountRow],
    exclusion_rows: list[BiogeographyReportExclusionRow],
    report_path: Path,
    tree_figure_path: Path,
    map_path: Path,
    legend_path: Path,
    caption_path: Path,
    summary_table_path: Path,
    region_count_table_path: Path,
    node_table_path: Path,
    transition_matrix_path: Path,
    event_table_path: Path,
    map_marker_table_path: Path,
    map_line_table_path: Path,
    exclusion_table_path: Path,
    warnings: list[str],
    audit: BiogeographyPublicationAudit,
    reproducibility_manifest_path: Path,
    reproducibility_manifest: dict[str, object],
) -> dict[str, object]:
    outputs = [
        report_path,
        tree_figure_path,
        map_path,
        legend_path,
        caption_path,
        summary_table_path,
        region_count_table_path,
        node_table_path,
        transition_matrix_path,
        event_table_path,
        map_marker_table_path,
        map_line_table_path,
        exclusion_table_path,
    ]
    return {
        "report_kind": "biogeography_report_package",
        "input_paths": [
            str(tree_path),
            str(traits_path),
            str(centroids_path),
        ],
        "input_checksums": {
            str(tree_path): checksum(tree_path),
            str(traits_path): checksum(traits_path),
            str(centroids_path): checksum(centroids_path),
        },
        "output_paths": [str(path) for path in outputs],
        "output_checksums": {str(path): checksum(path) for path in outputs},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": checksum(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "metrics": {
            "analyzed_taxon_count": state_report.summary.analyzed_taxon_count,
            "observed_region_count": state_report.summary.observed_region_count,
            "region_count_row_count": len(region_count_rows),
            "internal_node_count": state_report.summary.internal_node_count,
            "transition_rate_row_count": state_report.summary.transition_rate_row_count,
            "changed_branch_count": state_report.summary.changed_branch_count,
            "event_count": event_report.summary.event_count,
            "map_marker_count": len(map_report.marker_rows),
            "map_line_count": len(map_report.line_rows),
            "visible_map_line_count": map_report.summary.visible_line_count,
            "excluded_record_count": len(exclusion_rows),
            "root_region": state_report.summary.root_region,
            "warning_count": len(warnings),
            "publication_ready": audit.publication_ready,
            "legend_entry_count": audit.legend_entry_count,
            "rendered_internal_pie_count": audit.rendered_internal_pie_count,
            "rendered_internal_probability_label_count": audit.rendered_internal_probability_label_count,
        },
        "audit": {
            "publication_ready": audit.publication_ready,
            "legend_complete": audit.legend_complete,
            "caption_ready": audit.caption_ready,
            "node_probabilities_visible": audit.node_probabilities_visible,
            "transitions_visible": audit.transitions_visible,
            "map_state_colors_complete": audit.map_state_colors_complete,
            "rendered_internal_pie_count": audit.rendered_internal_pie_count,
            "rendered_internal_probability_label_count": audit.rendered_internal_probability_label_count,
            "expected_internal_node_count": audit.expected_internal_node_count,
            "visible_transition_count": audit.visible_transition_count,
            "state_color_count": audit.state_color_count,
            "legend_entry_count": audit.legend_entry_count,
            "reviewer_summary": audit.reviewer_summary,
            "limitations": audit.limitations,
        },
        "warnings": warnings,
    }

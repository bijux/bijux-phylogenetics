from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.biogeography.migration import (
    summarize_geographic_migration_events,
    write_geographic_migration_event_table,
)
from bijux_phylogenetics.biogeography.presentation.publication_support import (
    build_biogeography_caption_draft,
    build_biogeography_publication_audit,
    build_biogeography_publication_legend_entries,
    write_biogeography_caption,
    write_biogeography_publication_legend,
)
from bijux_phylogenetics.biogeography.state_models import (
    summarize_geographic_state_model,
    write_geographic_region_probability_table,
    write_geographic_state_summary_table,
    write_geographic_transition_rate_table,
)
from bijux_phylogenetics.comparative.discrete_evolution import (
    estimate_ancestral_geographic_states,
    render_tree_with_geographic_states,
)
from bijux_phylogenetics.phylogeography.geographic_map import (
    render_geographic_map_html,
    summarize_discrete_region_map,
    write_geographic_map_line_table,
    write_geographic_map_marker_table,
)
from bijux_phylogenetics.phylogeography.region_styles import (
    build_geographic_state_color_map,
)
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)

from .contracts import BiogeographyReportPackageResult
from .exclusions import (
    build_biogeography_report_exclusion_rows,
    write_biogeography_report_exclusion_table,
)
from .machine_manifest import build_machine_manifest
from .presentation import build_biogeography_report_html
from .region_counts import (
    summarize_biogeography_region_counts,
    write_biogeography_region_count_table,
)
from .reviewer_summary import build_reviewer_summary


def build_biogeography_report_package(
    *,
    tree_path: Path,
    traits_path: Path,
    centroids_path: Path,
    trait: str,
    out_dir: Path,
    taxon_column: str | None = None,
    model: str = "er",
    region_column: str = "region",
    latitude_column: str = "latitude",
    longitude_column: str = "longitude",
) -> BiogeographyReportPackageResult:
    """Build a full reviewer-facing biogeography report package."""
    out_dir.mkdir(parents=True, exist_ok=True)
    state_report = summarize_geographic_state_model(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    event_report = summarize_geographic_migration_events(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=model,
    )
    map_report = summarize_discrete_region_map(
        tree_path,
        traits_path,
        trait=trait,
        centroids_path=centroids_path,
        taxon_column=taxon_column,
        model=model,
        region_column=region_column,
        latitude_column=latitude_column,
        longitude_column=longitude_column,
    )
    reconstruction = estimate_ancestral_geographic_states(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=taxon_column,
        model=state_report.internal_model,
    )
    state_colors = build_geographic_state_color_map(
        [
            *[row.most_likely_region for row in state_report.node_rows],
            *[row.state_label for row in map_report.marker_rows],
        ]
    )

    tree_figure_path = out_dir / "ancestral-region-tree.svg"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    tree_render = render_tree_with_geographic_states(
        tree_path,
        reconstruction,
        out_path=tree_figure_path,
        layout="phylogram",
        state_colors=state_colors,
    )
    map_path = render_geographic_map_html(
        map_report,
        out_path=out_dir / "geographic-region-map.html",
        state_colors=state_colors,
    ).output_path

    summary_table_path = write_geographic_state_summary_table(
        out_dir / "summary.tsv",
        state_report,
    )
    excluded_taxa = {row.taxon for row in state_report.exclusion_rows}
    region_count_rows = summarize_biogeography_region_counts(
        tree_path,
        traits_path,
        trait=trait,
        taxon_column=state_report.taxon_column,
        excluded_taxa=excluded_taxa,
    )
    region_count_table_path = write_biogeography_region_count_table(
        out_dir / "region-counts.tsv",
        region_count_rows,
    )
    node_table_path = write_geographic_region_probability_table(
        out_dir / "ancestral-regions.tsv",
        state_report,
    )
    transition_matrix_path = write_geographic_transition_rate_table(
        out_dir / "transition-matrix.tsv",
        state_report,
    )
    event_table_path = write_geographic_migration_event_table(
        out_dir / "event-table.tsv",
        event_report,
    )
    map_marker_table_path = write_geographic_map_marker_table(
        out_dir / "map-markers.tsv",
        map_report,
    )
    map_line_table_path = write_geographic_map_line_table(
        out_dir / "map-lines.tsv",
        map_report,
    )
    exclusion_rows = build_biogeography_report_exclusion_rows(state_report, map_report)
    exclusion_table_path = write_biogeography_report_exclusion_table(
        out_dir / "exclusions.tsv",
        exclusion_rows,
    )
    legend_entries = build_biogeography_publication_legend_entries(
        state_report=state_report,
        map_report=map_report,
    )
    legend_path = write_biogeography_publication_legend(legend_path, legend_entries)
    audit = build_biogeography_publication_audit(
        state_report=state_report,
        map_report=map_report,
        tree_render=tree_render,
        legend_entries=legend_entries,
        exclusion_count=len(exclusion_rows),
    )
    caption_draft = build_biogeography_caption_draft(
        state_report=state_report,
        map_report=map_report,
        audit=audit,
    )
    caption_path = write_biogeography_caption(caption_path, caption_draft)

    warnings = list(
        dict.fromkeys(
            [*state_report.warnings, *event_report.warnings, *map_report.warnings]
        )
    )
    reviewer_summary, limitations = build_reviewer_summary(
        state_report=state_report,
        event_report=event_report,
        map_report=map_report,
        region_count_rows=region_count_rows,
        audit=audit,
    )
    limitations.extend(warnings)
    limitations.extend(audit.limitations)
    limitations = list(dict.fromkeys(limitations))
    report_path = out_dir / "biogeography-report.html"
    report_path.write_text(
        build_biogeography_report_html(
            state_report=state_report,
            event_report=event_report,
            tree_figure_path=tree_figure_path,
            map_path=map_path,
            region_count_rows=region_count_rows,
            reviewer_summary=reviewer_summary,
            limitations=limitations,
            audit=audit,
            caption_path=caption_path,
            legend_path=legend_path,
            legend_entries=legend_entries,
        ),
        encoding="utf-8",
    )

    manifest_path = out_dir / "biogeography-report.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="biogeography_report_package",
        input_files=[
            ("tree", tree_path),
            ("traits", traits_path),
            ("centroids", centroids_path),
        ],
        generated_figures=[
            ("ancestral_region_tree", tree_figure_path),
            ("geographic_map", map_path),
        ],
        generated_tables=[
            ("summary", summary_table_path),
            ("region_counts", region_count_table_path),
            ("ancestral_regions", node_table_path),
            ("transition_matrix", transition_matrix_path),
            ("event_table", event_table_path),
            ("map_markers", map_marker_table_path),
            ("map_lines", map_line_table_path),
            ("exclusions", exclusion_table_path),
            ("legend", legend_path),
        ],
        filters=None,
        model={
            "kind": "discrete_geography",
            "name": model,
            "trait": trait,
            "internal_model": state_report.internal_model,
        },
        settings={
            "taxon_column": taxon_column,
            "region_column": region_column,
            "latitude_column": latitude_column,
            "longitude_column": longitude_column,
        },
        linked_artifacts=[
            ("caption", caption_path),
            ("review", report_path),
        ],
    )
    machine_manifest = build_machine_manifest(
        tree_path=tree_path,
        traits_path=traits_path,
        centroids_path=centroids_path,
        state_report=state_report,
        event_report=event_report,
        map_report=map_report,
        region_count_rows=region_count_rows,
        exclusion_rows=exclusion_rows,
        report_path=report_path,
        tree_figure_path=tree_figure_path,
        map_path=map_path,
        legend_path=legend_path,
        caption_path=caption_path,
        summary_table_path=summary_table_path,
        region_count_table_path=region_count_table_path,
        node_table_path=node_table_path,
        transition_matrix_path=transition_matrix_path,
        event_table_path=event_table_path,
        map_marker_table_path=map_marker_table_path,
        map_line_table_path=map_line_table_path,
        exclusion_table_path=exclusion_table_path,
        warnings=warnings,
        audit=audit,
        reproducibility_manifest_path=reproducibility_manifest_path,
        reproducibility_manifest=reproducibility_manifest,
    )
    manifest_path.write_text(
        json.dumps(machine_manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return BiogeographyReportPackageResult(
        output_dir=out_dir,
        report_path=report_path,
        tree_figure_path=tree_figure_path,
        map_path=map_path,
        legend_path=legend_path,
        caption_path=caption_path,
        summary_table_path=summary_table_path,
        region_count_table_path=region_count_table_path,
        node_table_path=node_table_path,
        transition_matrix_path=transition_matrix_path,
        event_table_path=event_table_path,
        map_marker_table_path=map_marker_table_path,
        map_line_table_path=map_line_table_path,
        exclusion_table_path=exclusion_table_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        state_report=state_report,
        event_report=event_report,
        map_report=map_report,
        tree_render=tree_render,
        region_count_rows=region_count_rows,
        exclusion_rows=exclusion_rows,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
        warnings=warnings,
        machine_manifest=machine_manifest,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
    )

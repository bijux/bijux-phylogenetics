from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.biogeography.state_models import (
    GeographicStateModelReport,
    summarize_geographic_state_model,
    write_geographic_region_probability_table,
    write_geographic_state_summary_table,
    write_geographic_transition_rate_table,
)
from bijux_phylogenetics.biogeography.migration import (
    GeographicMigrationEventReport,
    summarize_geographic_migration_events,
    write_geographic_migration_event_table,
)
from bijux_phylogenetics.biogeography.publication import (
    BiogeographyCaptionDraft,
    BiogeographyPublicationAudit,
    BiogeographyPublicationLegendEntry,
    build_biogeography_caption_draft,
    build_biogeography_publication_audit,
    build_biogeography_publication_legend_entries,
    write_biogeography_caption,
    write_biogeography_publication_legend,
)
from bijux_phylogenetics.core.metadata import load_taxon_table, write_taxon_rows
from bijux_phylogenetics.comparative.discrete_evolution import (
    estimate_ancestral_geographic_states,
    render_tree_with_geographic_states,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylogeography.geographic_map import (
    GeographicMapReport,
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
from bijux_phylogenetics.render.svg import TreeRenderResult


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


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def summarize_biogeography_region_counts(
    tree_path: Path,
    table_path: Path,
    *,
    trait: str,
    taxon_column: str,
    excluded_taxa: set[str],
) -> list[BiogeographyRegionCountRow]:
    """Count observed regions after tree overlap and exclusion auditing."""
    tree = load_tree(tree_path)
    table = load_taxon_table(table_path, taxon_column=taxon_column)
    rows_by_taxon = {row[table.taxon_column]: row for row in table.rows}
    counts: dict[str, int] = {}
    analyzed_taxa = 0
    for taxon in tree.tip_names:
        row = rows_by_taxon.get(taxon)
        if row is None or taxon in excluded_taxa:
            continue
        region = row.get(trait, "").strip()
        if not region:
            continue
        counts[region] = counts.get(region, 0) + 1
        analyzed_taxa += 1
    if analyzed_taxa == 0:
        return []
    return [
        BiogeographyRegionCountRow(
            region=region,
            tip_taxon_count=count,
            analyzed_taxon_fraction=count / analyzed_taxa,
        )
        for region, count in sorted(counts.items())
    ]


def write_biogeography_region_count_table(
    path: Path,
    rows: list[BiogeographyRegionCountRow],
) -> Path:
    """Write one observed-region count ledger."""
    return write_taxon_rows(
        path,
        columns=[
            "region",
            "tip_taxon_count",
            "analyzed_taxon_fraction",
        ],
        rows=[
            {
                "region": row.region,
                "tip_taxon_count": str(row.tip_taxon_count),
                "analyzed_taxon_fraction": str(row.analyzed_taxon_fraction),
            }
            for row in rows
        ],
    )


def write_biogeography_report_exclusion_table(
    path: Path,
    rows: list[BiogeographyReportExclusionRow],
) -> Path:
    """Write one combined exclusion ledger for the full report package."""
    return write_taxon_rows(
        path,
        columns=[
            "surface",
            "subject_id",
            "subject_kind",
            "raw_left",
            "raw_right",
            "reason",
            "note",
        ],
        rows=[
            {
                "surface": row.surface,
                "subject_id": row.subject_id,
                "subject_kind": row.subject_kind,
                "raw_left": row.raw_left,
                "raw_right": row.raw_right,
                "reason": row.reason,
                "note": row.note,
            }
            for row in rows
        ],
    )


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
    exclusion_rows = _build_exclusion_rows(state_report, map_report)
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
    reviewer_summary, limitations = _reviewer_summary(
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
        _build_report_html(
            state_report=state_report,
            event_report=event_report,
            map_report=map_report,
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
    machine_manifest = _build_machine_manifest(
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


def _build_exclusion_rows(
    state_report: GeographicStateModelReport,
    map_report: GeographicMapReport,
) -> list[BiogeographyReportExclusionRow]:
    rows: list[BiogeographyReportExclusionRow] = [
        BiogeographyReportExclusionRow(
            surface="state_model",
            subject_id=row.taxon,
            subject_kind="taxon",
            raw_left=row.raw_region,
            raw_right=row.normalized_region or "",
            reason=row.reason,
            note=row.note,
        )
        for row in state_report.exclusion_rows
    ]
    rows.extend(
        BiogeographyReportExclusionRow(
            surface="map",
            subject_id=row.subject_id,
            subject_kind=row.subject_kind,
            raw_left=row.raw_left,
            raw_right=row.raw_right,
            reason=row.reason,
            note=row.note,
        )
        for row in map_report.exclusion_rows
    )
    deduplicated: list[BiogeographyReportExclusionRow] = []
    seen: set[tuple[str, str, str, str, str, str, str]] = set()
    for row in rows:
        key = (
            row.surface,
            row.subject_id,
            row.subject_kind,
            row.raw_left,
            row.raw_right,
            row.reason,
            row.note,
        )
        if key in seen:
            continue
        seen.add(key)
        deduplicated.append(row)
    return deduplicated


def _reviewer_summary(
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


def _build_machine_manifest(
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
            str(tree_path): _checksum(tree_path),
            str(traits_path): _checksum(traits_path),
            str(centroids_path): _checksum(centroids_path),
        },
        "output_paths": [str(path) for path in outputs],
        "output_checksums": {str(path): _checksum(path) for path in outputs},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(
            reproducibility_manifest_path
        ),
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


def _build_report_html(
    *,
    state_report: GeographicStateModelReport,
    event_report: GeographicMigrationEventReport,
    map_report: GeographicMapReport,
    tree_figure_path: Path,
    map_path: Path,
    region_count_rows: list[BiogeographyRegionCountRow],
    reviewer_summary: list[str],
    limitations: list[str],
    audit: BiogeographyPublicationAudit,
    caption_path: Path,
    legend_path: Path,
    legend_entries: list[BiogeographyPublicationLegendEntry],
) -> str:
    tree_svg = tree_figure_path.read_text(encoding="utf-8")
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Biogeography Report</title>",
            "  <style>",
            "    body { font-family: Georgia, 'Times New Roman', serif; margin: 0; background: linear-gradient(180deg, #f5f7f3 0%, #e7ede7 100%); color: #193426; }",
            "    main { max-width: 1320px; margin: 0 auto; padding: 24px; }",
            "    h1 { margin: 0 0 8px; font-size: 32px; }",
            "    h2 { margin: 0 0 10px; font-size: 22px; }",
            "    p { line-height: 1.5; }",
            "    .cards { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 14px; margin: 18px 0 24px; }",
            "    .card, .panel { background: rgba(255,255,255,0.82); border: 1px solid rgba(25,52,38,0.14); border-radius: 18px; padding: 18px; box-shadow: 0 16px 42px rgba(25,52,38,0.08); }",
            "    .card strong { display: block; font-size: 22px; margin-top: 6px; }",
            "    .label { color: #587363; font-size: 13px; text-transform: uppercase; letter-spacing: 0.04em; }",
            "    .grid { display: grid; grid-template-columns: 1.05fr 0.95fr; gap: 20px; }",
            "    .figure-shell { overflow: auto; }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; margin-top: 10px; font-size: 14px; }",
            "    th, td { border-bottom: 1px solid rgba(25,52,38,0.10); padding: 8px 10px; text-align: left; vertical-align: top; }",
            "    th { color: #365443; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #16543a; }",
            "    iframe { width: 100%; min-height: 780px; border: 1px solid rgba(25,52,38,0.14); border-radius: 14px; background: white; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Biogeography Report</h1>",
            "  <p>Complete geographic evolution review for one rooted tree, one region table, and one explicit centroid table. The bundle keeps the ancestral-region tree, transition matrix, migration-event ledger, and geographic map together.</p>",
            '  <section class="cards">',
            f'    <div class="card"><span class="label">model</span><strong>{escape(state_report.model)}</strong></div>',
            f'    <div class="card"><span class="label">root region</span><strong>{escape(state_report.summary.root_region)}</strong></div>',
            f'    <div class="card"><span class="label">transition rows</span><strong>{state_report.summary.transition_rate_row_count}</strong></div>',
            f'    <div class="card"><span class="label">event rows</span><strong>{event_report.summary.event_count}</strong></div>',
            "  </section>",
            '  <section class="panel">',
            "    <h2>Reviewer Summary</h2>",
            _list(reviewer_summary),
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Publication Audit</h2>",
            _publication_audit_table(audit),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Figure Legend and Caption</h2>",
            f'      <p>The publication legend is exported to <a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a> and the caption draft to <a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a>.</p>',
            _legend_entry_table(legend_entries),
            "    </section>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Ancestral-Region Tree</h2>",
            "      <p>The tree figure shows tip regions and internal ancestral region calls from the owned geographic-state reconstruction.</p>",
            '      <div class="figure-shell">',
            tree_svg,
            "      </div>",
            "    </section>",
            '    <section class="panel">',
            "      <h2>Geographic Map</h2>",
            f'      <p>The package includes the self-contained map artifact <a href="{escape(map_path.name)}">{escape(map_path.name)}</a>.</p>',
            f'      <iframe src="{escape(map_path.name)}" title="Geographic map review"></iframe>',
            "    </section>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Observed Region Counts</h2>",
            _region_count_table(region_count_rows),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Transition Matrix Summary</h2>",
            _transition_rate_table(state_report),
            "    </section>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel">',
            "      <h2>Ancestral Region Nodes</h2>",
            _node_probability_table(state_report),
            "    </section>",
            '    <section class="panel">',
            "      <h2>Migration Events</h2>",
            _event_table(event_report),
            "    </section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Limitations</h2>",
            _list(limitations),
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )


def _list(items: list[str]) -> str:
    return "<ul>" + "".join(f"<li>{escape(item)}</li>" for item in items) + "</ul>"


def _region_count_table(rows: list[BiogeographyRegionCountRow]) -> str:
    return _table(
        headers=["region", "tip_taxon_count", "analyzed_taxon_fraction"],
        rows=[
            [
                row.region,
                str(row.tip_taxon_count),
                f"{row.analyzed_taxon_fraction:.3f}",
            ]
            for row in rows
        ],
    )


def _transition_rate_table(report: GeographicStateModelReport) -> str:
    return _table(
        headers=[
            "source_region",
            "target_region",
            "rate",
            "lower_95_interval",
            "upper_95_interval",
        ],
        rows=[
            [
                row.source_region,
                row.target_region,
                f"{row.rate:.6f}",
                f"{row.lower_95_interval:.6f}",
                f"{row.upper_95_interval:.6f}",
            ]
            for row in report.transition_rate_rows
        ],
    )


def _node_probability_table(report: GeographicStateModelReport) -> str:
    return _table(
        headers=["node", "descendant_taxa", "most_likely_region", "confidence"],
        rows=[
            [
                row.node_name or row.node,
                ", ".join(row.descendant_taxa),
                row.most_likely_region,
                f"{row.confidence:.3f}",
            ]
            for row in report.node_rows
        ],
    )


def _event_table(report: GeographicMigrationEventReport) -> str:
    return _table(
        headers=[
            "branch_id",
            "source_region",
            "target_region",
            "support",
            "midpoint_depth",
        ],
        rows=[
            [
                row.branch_id,
                row.source_region,
                row.target_region,
                f"{row.support:.3f}",
                f"{row.midpoint_depth:.6f}",
            ]
            for row in report.event_rows
        ],
    )


def _publication_audit_table(audit: BiogeographyPublicationAudit) -> str:
    return _table(
        headers=["metric", "value"],
        rows=[
            ["publication_ready", str(audit.publication_ready).lower()],
            ["legend_complete", str(audit.legend_complete).lower()],
            ["caption_ready", str(audit.caption_ready).lower()],
            ["node_probabilities_visible", str(audit.node_probabilities_visible).lower()],
            ["transitions_visible", str(audit.transitions_visible).lower()],
            ["map_state_colors_complete", str(audit.map_state_colors_complete).lower()],
            [
                "rendered_internal_pies",
                f"{audit.rendered_internal_pie_count}/{audit.expected_internal_node_count}",
            ],
            [
                "rendered_internal_probability_labels",
                f"{audit.rendered_internal_probability_label_count}/{audit.expected_internal_node_count}",
            ],
            ["visible_transition_count", str(audit.visible_transition_count)],
            ["state_color_count", str(audit.state_color_count)],
        ],
    )


def _legend_entry_table(entries: list[BiogeographyPublicationLegendEntry]) -> str:
    return _table(
        headers=["surface", "label", "swatch", "detail"],
        rows=[
            [entry.surface, entry.label, entry.swatch, entry.detail]
            for entry in entries
        ],
    )


def _table(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{escape(header)}</th>" for header in headers)
    body = "".join(
        "<tr>" + "".join(f"<td>{escape(cell)}</td>" for cell in row) + "</tr>"
        for row in rows
    )
    return f"<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>"

from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.biogeography.migration import GeographicMigrationEventReport
from bijux_phylogenetics.biogeography.presentation.publication_support import (
    BiogeographyPublicationAudit,
    BiogeographyPublicationLegendEntry,
)
from bijux_phylogenetics.biogeography.state_models import GeographicStateModelReport

from .contracts import BiogeographyRegionCountRow


def build_biogeography_report_html(
    *,
    state_report: GeographicStateModelReport,
    event_report: GeographicMigrationEventReport,
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
            [
                "node_probabilities_visible",
                str(audit.node_probabilities_visible).lower(),
            ],
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

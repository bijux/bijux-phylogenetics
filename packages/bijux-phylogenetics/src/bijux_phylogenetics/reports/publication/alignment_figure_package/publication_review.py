from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.phylo.alignment import (
    AlignmentForensicReport,
    AlignmentSummary,
)
from bijux_phylogenetics.reports.review import ReviewerAuditChecklist

from .contracts import (
    AlignmentFigureAudit,
    AlignmentFigureCaptionDraft,
    AlignmentFigureLegendEntry,
)


def build_audit(
    *,
    summary: AlignmentSummary,
    forensic: AlignmentForensicReport,
    heatmap_row_count: int,
    heatmap_bin_count: int,
    plotted_window_count: int,
    plotted_sequence_count: int,
    legend_entries: list[AlignmentFigureLegendEntry],
) -> AlignmentFigureAudit:
    """Build the reviewer-facing publication audit for one figure package."""
    heatmap_visible = heatmap_row_count > 0 and heatmap_bin_count > 0
    site_summary_visible = plotted_window_count > 0
    sequence_panel_visible = plotted_sequence_count > 0
    legend_complete = {entry.surface for entry in legend_entries} == {
        "missingness-heatmap",
        "site-quality-summary",
        "sequence-quality-panel",
    }
    caption_ready = heatmap_visible and site_summary_visible and sequence_panel_visible
    suspicious_alignment = forensic.quality.suspicious_alignment
    publication_ready = (
        caption_ready
        and legend_complete
        and forensic.quality.quality_score >= 75.0
        and not suspicious_alignment
        and not summary.invalid_characters
    )
    reviewer_summary = [
        f"alignment quality score: {forensic.quality.quality_score}",
        f"heatmap rows and bins: {heatmap_row_count} x {heatmap_bin_count}",
        f"window summaries plotted: {plotted_window_count}",
        f"sequence rows plotted: {plotted_sequence_count}",
    ]
    limitations: list[str] = []
    if suspicious_alignment:
        limitations.extend(forensic.quality.suspicious_reasons)
    if summary.invalid_characters:
        limitations.append(
            "alignment contains invalid characters for the inferred alphabet"
        )
    if forensic.quality.quality_score < 75.0:
        limitations.append(
            "alignment quality score remains below the reviewer threshold"
        )
    if not summary.near_duplicate_scan_performed:
        limitations.append(
            "near-duplicate sequence review was skipped because the alignment exceeded the governed pairwise scan threshold"
        )
    if not heatmap_visible:
        limitations.append(
            "the package does not currently render the missingness heatmap"
        )
    if not site_summary_visible:
        limitations.append(
            "the package does not currently render the site-quality summary"
        )
    if not sequence_panel_visible:
        limitations.append(
            "the package does not currently render the sequence-quality panel"
        )
    if not limitations:
        limitations.append(
            "the current package keeps the key alignment quality figures explicit enough for publication-oriented review"
        )
    return AlignmentFigureAudit(
        publication_ready=publication_ready,
        heatmap_visible=heatmap_visible,
        site_summary_visible=site_summary_visible,
        sequence_panel_visible=sequence_panel_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        suspicious_alignment=suspicious_alignment,
        quality_score=forensic.quality.quality_score,
        heatmap_row_count=heatmap_row_count,
        heatmap_bin_count=heatmap_bin_count,
        plotted_window_count=plotted_window_count,
        plotted_sequence_count=plotted_sequence_count,
        invalid_character_count=len(summary.invalid_characters),
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def build_caption_draft(
    *,
    summary: AlignmentSummary,
    audit: AlignmentFigureAudit,
) -> AlignmentFigureCaptionDraft:
    """Build the durable caption draft for one alignment figure package."""
    return AlignmentFigureCaptionDraft(
        title="Alignment quality review across missingness, site windows, and sequence burden",
        lead_sentence=(
            f"This package summarizes one {summary.inferred_alphabet} alignment with {summary.sequence_count} sequences and {summary.alignment_length} aligned sites through three explicit reviewer figures rather than leaving quality evidence buried in tables alone."
        ),
        heatmap_sentence=(
            f"The missingness heatmap keeps {audit.heatmap_row_count} sequence rows and {audit.heatmap_bin_count} site bins visible so missing-data concentration can be reviewed directly."
        ),
        site_summary_sentence=(
            f"The site-quality summary renders {audit.plotted_window_count} sliding windows and highlights suspicious over- or under-aligned regions on the figure itself."
        ),
        sequence_panel_sentence=(
            f"The sequence-quality panel keeps {audit.plotted_sequence_count} ranked sequence burdens explicit, with lower scores reflecting missingness, ambiguity, duplicate burden, or composition risk."
        ),
        limitation_sentence=audit.limitations[0],
        caption_ready=audit.caption_ready,
    )


def write_caption(path: Path, draft: AlignmentFigureCaptionDraft) -> Path:
    """Write the figure caption artifact for one alignment review package."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.heatmap_sentence,
                draft.site_summary_sentence,
                draft.sequence_panel_sentence,
                draft.limitation_sentence,
                "",
                f"caption_ready: {'true' if draft.caption_ready else 'false'}",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return path


def build_review_html(
    *,
    heatmap_figure_path: Path,
    site_summary_figure_path: Path,
    sequence_panel_figure_path: Path,
    heatmap_table_path: Path,
    window_table_path: Path,
    ranking_table_path: Path,
    legend_path: Path,
    caption_path: Path,
    reviewer_audit_checklist_path: Path,
    audit: AlignmentFigureAudit,
    reviewer_audit_checklist: ReviewerAuditChecklist,
) -> str:
    """Render the reviewer-facing HTML package review."""
    figures = {
        "heatmap": heatmap_figure_path.read_text(encoding="utf-8"),
        "site_summary": site_summary_figure_path.read_text(encoding="utf-8"),
        "sequence_panel": sequence_panel_figure_path.read_text(encoding="utf-8"),
    }
    audit_rows = "".join(
        "<tr><th>" + escape(label) + "</th><td>" + escape(value) + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("quality_score", format(audit.quality_score, ".15g")),
            ("suspicious_alignment", str(audit.suspicious_alignment).lower()),
            ("heatmap_visible", str(audit.heatmap_visible).lower()),
            ("site_summary_visible", str(audit.site_summary_visible).lower()),
            ("sequence_panel_visible", str(audit.sequence_panel_visible).lower()),
        ]
    )
    limitation_items = "".join(f"<li>{escape(item)}</li>" for item in audit.limitations)
    checklist_rows = "".join(
        "<tr><td>"
        + escape(item.section)
        + "</td><td>"
        + escape(item.status)
        + "</td><td>"
        + escape(item.summary)
        + "</td><td>"
        + escape("; ".join(item.evidence))
        + "</td></tr>"
        for item in reviewer_audit_checklist.items
    )
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Alignment Quality Review</title>",
            "  <style>",
            "    body { margin: 0; background: linear-gradient(180deg, #eef6f4 0%, #f8fafc 100%); color: #1b1f24; font: 16px/1.5 'Iowan Old Style', 'Palatino Linotype', serif; }",
            "    main { max-width: 1220px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { font-family: 'Avenir Next', 'Segoe UI', sans-serif; }",
            "    h1 { color: #0f766e; margin-top: 0; }",
            "    .grid { display: grid; grid-template-columns: 1fr; gap: 18px; }",
            "    .panel { background: rgba(255,255,255,0.84); border: 1px solid rgba(15,118,110,0.14); border-radius: 18px; padding: 18px; box-shadow: 0 18px 42px rgba(15,118,110,0.08); }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid rgba(15,118,110,0.12); vertical-align: top; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #0f766e; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Alignment Quality Review</h1>",
            "  <p>Reviewer-facing alignment figure package with one missingness heatmap, one site-quality summary, and one sequence-quality panel, backed by explicit ledgers and publication-oriented audit fields.</p>",
            '  <section class="panel">',
            "    <h2>Publication Audit</h2>",
            f"    <table><tbody>{audit_rows}</tbody></table>",
            "    <ul>" + limitation_items + "</ul>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Reviewer Audit Checklist</h2>",
            "    <table><thead><tr><th>section</th><th>status</th><th>summary</th><th>evidence</th></tr></thead><tbody>"
            + checklist_rows
            + "</tbody></table>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel"><h2>Missingness Heatmap</h2><div class="figure-shell">'
            + figures["heatmap"]
            + "</div></section>",
            '    <section class="panel"><h2>Site-Quality Summary</h2><div class="figure-shell">'
            + figures["site_summary"]
            + "</div></section>",
            '    <section class="panel"><h2>Sequence-Quality Panel</h2><div class="figure-shell">'
            + figures["sequence_panel"]
            + "</div></section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Linked Artifacts</h2>",
            "    <ul>",
            f'      <li><a href="{escape(heatmap_table_path.name)}">{escape(heatmap_table_path.name)}</a></li>',
            f'      <li><a href="{escape(window_table_path.name)}">{escape(window_table_path.name)}</a></li>',
            f'      <li><a href="{escape(ranking_table_path.name)}">{escape(ranking_table_path.name)}</a></li>',
            f'      <li><a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a></li>',
            f'      <li><a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a></li>',
            f'      <li><a href="{escape(reviewer_audit_checklist_path.name)}">{escape(reviewer_audit_checklist_path.name)}</a></li>',
            "    </ul>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )

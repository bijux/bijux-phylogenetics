from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_forensic_report,
    summarize_alignment_windows,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.phylo.alignment import (
    AlignmentForensicReport,
    AlignmentRecord,
    AlignmentSummary,
    AlignmentWindowSummary,
    SequenceQualityRankingRow,
)
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.reports.review import (
    ReviewerAuditChecklist,
    build_reviewer_audit_checklist,
    write_reviewer_audit_checklist,
)
from .contracts import (
    AlignmentFigureAudit,
    AlignmentFigureCaptionDraft,
    AlignmentFigureLegendEntry,
    AlignmentFigurePackageResult,
)
from .artifact_outputs import (
    build_legend_entries,
    write_heatmap_table,
    write_legend_table,
    write_missingness_heatmap,
    write_ranking_table,
    write_sequence_quality_panel,
    write_site_quality_summary,
    write_window_table,
)
from .heatmap_analysis import build_heatmap_cells


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _build_audit(
    *,
    summary: AlignmentSummary,
    forensic: AlignmentForensicReport,
    heatmap_row_count: int,
    heatmap_bin_count: int,
    plotted_window_count: int,
    plotted_sequence_count: int,
    legend_entries: list[AlignmentFigureLegendEntry],
) -> AlignmentFigureAudit:
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


def _build_caption_draft(
    *,
    summary: AlignmentSummary,
    audit: AlignmentFigureAudit,
) -> AlignmentFigureCaptionDraft:
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


def _write_caption(path: Path, draft: AlignmentFigureCaptionDraft) -> Path:
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


def _build_review_html(
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


def build_alignment_figure_package(
    alignment_path: Path,
    *,
    out_dir: Path,
    maximum_site_bins: int = 120,
    window_size: int = 30,
    step_size: int = 10,
    panel_row_limit: int = 12,
) -> AlignmentFigurePackageResult:
    """Build a reviewer-facing alignment figure package with explicit quality figures."""
    out_dir.mkdir(parents=True, exist_ok=True)
    heatmap_figure_path = out_dir / "alignment-missingness-heatmap.svg"
    site_summary_figure_path = out_dir / "alignment-site-quality-summary.svg"
    sequence_panel_figure_path = out_dir / "alignment-sequence-quality-panel.svg"
    heatmap_table_path = out_dir / "alignment-missingness-heatmap.tsv"
    window_table_path = out_dir / "alignment-site-quality-windows.tsv"
    ranking_table_path = out_dir / "alignment-sequence-quality-ranking.tsv"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    review_path = out_dir / "alignment-quality-review.html"
    manifest_path = out_dir / "alignment-quality-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"
    reviewer_audit_checklist_path = out_dir / "reviewer-audit-checklist.tsv"

    summary = summarise_fasta(alignment_path)
    records = load_fasta_alignment(alignment_path)
    forensic = build_alignment_forensic_report(alignment_path)
    windows = summarize_alignment_windows(
        alignment_path,
        window_size=window_size,
        step_size=step_size,
    )
    heatmap_cells, heatmap_row_count, heatmap_bin_count = build_heatmap_cells(
        summary,
        records,
        forensic.sequence_ranking.rows,
        maximum_bins=maximum_site_bins,
    )
    write_heatmap_table(heatmap_table_path, heatmap_cells)
    write_window_table(
        window_table_path,
        windows=windows,
        over_regions=forensic.over_aligned_regions,
        under_regions=forensic.under_aligned_regions,
    )
    write_ranking_table(ranking_table_path, forensic.sequence_ranking.rows)
    heatmap_row_count, heatmap_bin_count = write_missingness_heatmap(
        heatmap_figure_path,
        cells=heatmap_cells,
        ranking_rows=forensic.sequence_ranking.rows,
        heatmap_bin_count=heatmap_bin_count,
    )
    plotted_window_count = write_site_quality_summary(
        site_summary_figure_path,
        windows=windows,
        over_regions=forensic.over_aligned_regions,
        under_regions=forensic.under_aligned_regions,
    )
    plotted_sequence_count = write_sequence_quality_panel(
        sequence_panel_figure_path,
        ranking_rows=forensic.sequence_ranking.rows,
        maximum_rows=panel_row_limit,
    )
    legend_entries = build_legend_entries()
    write_legend_table(legend_path, legend_entries)
    audit = _build_audit(
        summary=summary,
        forensic=forensic,
        heatmap_row_count=heatmap_row_count,
        heatmap_bin_count=heatmap_bin_count,
        plotted_window_count=plotted_window_count,
        plotted_sequence_count=plotted_sequence_count,
        legend_entries=legend_entries,
    )
    caption_draft = _build_caption_draft(summary=summary, audit=audit)
    _write_caption(caption_path, caption_draft)
    artifact_paths = [
        heatmap_figure_path,
        site_summary_figure_path,
        sequence_panel_figure_path,
        heatmap_table_path,
        window_table_path,
        ranking_table_path,
        legend_path,
        caption_path,
        review_path,
    ]
    existing_artifact_paths = artifact_paths[:-1]
    pre_review_manifest = {
        "report_kind": "alignment_quality_figure_package",
        "input_path": str(alignment_path),
        "input_checksum": _checksum(alignment_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {
            str(path): _checksum(path) for path in existing_artifact_paths
        },
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "settings": {
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
        },
        "metrics": {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "quality_score": forensic.quality.quality_score,
            "publication_ready": audit.publication_ready,
            "heatmap_row_count": audit.heatmap_row_count,
            "heatmap_bin_count": audit.heatmap_bin_count,
            "plotted_window_count": audit.plotted_window_count,
            "plotted_sequence_count": audit.plotted_sequence_count,
        },
        "alignment_summary": _json_ready(asdict(summary)),
        "alignment_quality": _json_ready(asdict(forensic.quality)),
        "alignment_readiness": _json_ready(asdict(forensic.readiness)),
        "alignment_low_information": _json_ready(asdict(forensic.low_information)),
        "audit": _json_ready(asdict(audit)),
    }
    reviewer_audit_checklist = build_reviewer_audit_checklist(pre_review_manifest)
    review_path.write_text(
        _build_review_html(
            heatmap_figure_path=heatmap_figure_path,
            site_summary_figure_path=site_summary_figure_path,
            sequence_panel_figure_path=sequence_panel_figure_path,
            heatmap_table_path=heatmap_table_path,
            window_table_path=window_table_path,
            ranking_table_path=ranking_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            reviewer_audit_checklist_path=reviewer_audit_checklist_path,
            audit=audit,
            reviewer_audit_checklist=reviewer_audit_checklist,
        ),
        encoding="utf-8",
    )
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="alignment_quality_figure_package",
        input_files=[("alignment", alignment_path)],
        generated_figures=[
            ("missingness_heatmap", heatmap_figure_path),
            ("site_quality_summary", site_summary_figure_path),
            ("sequence_quality_panel", sequence_panel_figure_path),
        ],
        generated_tables=[
            ("missingness_heatmap", heatmap_table_path),
            ("site_quality_windows", window_table_path),
            ("sequence_quality_ranking", ranking_table_path),
            ("legend", legend_path),
        ],
        filters=None,
        model={
            "kind": "alignment_quality",
            "name": "summary-and-forensic-review",
        },
        settings={
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
            "alignment_length": summary.alignment_length,
            "sequence_count": summary.sequence_count,
        },
        linked_artifacts=[
            ("caption", caption_path),
            ("review", review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "alignment_quality_figure_package",
        "input_path": str(alignment_path),
        "input_checksum": _checksum(alignment_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "maximum_site_bins": maximum_site_bins,
            "window_size": window_size,
            "step_size": step_size,
            "panel_row_limit": panel_row_limit,
        },
        "metrics": {
            "sequence_count": summary.sequence_count,
            "alignment_length": summary.alignment_length,
            "quality_score": forensic.quality.quality_score,
            "publication_ready": audit.publication_ready,
            "heatmap_row_count": audit.heatmap_row_count,
            "heatmap_bin_count": audit.heatmap_bin_count,
            "plotted_window_count": audit.plotted_window_count,
            "plotted_sequence_count": audit.plotted_sequence_count,
        },
        "alignment_summary": _json_ready(asdict(summary)),
        "alignment_quality": _json_ready(asdict(forensic.quality)),
        "alignment_readiness": _json_ready(asdict(forensic.readiness)),
        "alignment_low_information": _json_ready(asdict(forensic.low_information)),
        "audit": _json_ready(asdict(audit)),
    }
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest["output_paths"].append(str(reviewer_audit_checklist_path))
    machine_manifest["output_checksums"][str(reviewer_audit_checklist_path)] = (
        _checksum(reviewer_audit_checklist_path)
    )
    machine_manifest["reviewer_audit_checklist_path"] = str(
        reviewer_audit_checklist_path
    )
    machine_manifest["reviewer_audit_checklist"] = _json_ready(
        asdict(reviewer_audit_checklist)
    )
    review_path.write_text(
        _build_review_html(
            heatmap_figure_path=heatmap_figure_path,
            site_summary_figure_path=site_summary_figure_path,
            sequence_panel_figure_path=sequence_panel_figure_path,
            heatmap_table_path=heatmap_table_path,
            window_table_path=window_table_path,
            ranking_table_path=ranking_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            reviewer_audit_checklist_path=reviewer_audit_checklist_path,
            audit=audit,
            reviewer_audit_checklist=reviewer_audit_checklist,
        ),
        encoding="utf-8",
    )
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return AlignmentFigurePackageResult(
        output_dir=out_dir,
        heatmap_figure_path=heatmap_figure_path,
        site_summary_figure_path=site_summary_figure_path,
        sequence_panel_figure_path=sequence_panel_figure_path,
        heatmap_table_path=heatmap_table_path,
        window_table_path=window_table_path,
        ranking_table_path=ranking_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        summary=summary,
        forensic=forensic,
        windows=windows,
        heatmap_cells=heatmap_cells,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
        reviewer_audit_checklist=reviewer_audit_checklist,
        machine_manifest=machine_manifest,
    )

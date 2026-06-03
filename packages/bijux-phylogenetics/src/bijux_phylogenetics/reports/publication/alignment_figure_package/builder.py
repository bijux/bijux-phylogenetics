from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.fasta.quality import (
    build_alignment_forensic_report,
    summarize_alignment_windows,
)
from bijux_phylogenetics.io.fasta.records import summarise_fasta
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.reports.review import (
    build_reviewer_audit_checklist,
    write_reviewer_audit_checklist,
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
from .contracts import AlignmentFigurePackageResult
from .heatmap_analysis import build_heatmap_cells
from .manifest import (
    attach_reviewer_audit_checklist,
    build_machine_manifest,
    build_pre_review_manifest,
    write_machine_manifest,
)
from .publication_review import (
    build_audit,
    build_caption_draft,
    build_review_html,
    write_caption,
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
    audit = build_audit(
        summary=summary,
        forensic=forensic,
        heatmap_row_count=heatmap_row_count,
        heatmap_bin_count=heatmap_bin_count,
        plotted_window_count=plotted_window_count,
        plotted_sequence_count=plotted_sequence_count,
        legend_entries=legend_entries,
    )
    caption_draft = build_caption_draft(summary=summary, audit=audit)
    write_caption(caption_path, caption_draft)
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
    pre_review_manifest = build_pre_review_manifest(
        alignment_path=alignment_path,
        artifact_paths=artifact_paths,
        reproducibility_manifest_path=reproducibility_manifest_path,
        maximum_site_bins=maximum_site_bins,
        window_size=window_size,
        step_size=step_size,
        panel_row_limit=panel_row_limit,
        summary=summary,
        forensic=forensic,
        audit=audit,
    )
    reviewer_audit_checklist = build_reviewer_audit_checklist(pre_review_manifest)
    review_path.write_text(
        build_review_html(
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
    machine_manifest = build_machine_manifest(
        alignment_path=alignment_path,
        artifact_paths=artifact_paths,
        reproducibility_manifest_path=reproducibility_manifest_path,
        reproducibility_manifest=reproducibility_manifest,
        maximum_site_bins=maximum_site_bins,
        window_size=window_size,
        step_size=step_size,
        panel_row_limit=panel_row_limit,
        summary=summary,
        forensic=forensic,
        audit=audit,
    )
    reviewer_audit_checklist = write_reviewer_audit_checklist(
        reviewer_audit_checklist_path,
        machine_manifest,
    ).checklist
    machine_manifest = attach_reviewer_audit_checklist(
        machine_manifest=machine_manifest,
        reviewer_audit_checklist_path=reviewer_audit_checklist_path,
        reviewer_audit_checklist=reviewer_audit_checklist,
    )
    review_path.write_text(
        build_review_html(
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
    write_machine_manifest(
        manifest_path,
        machine_manifest,
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

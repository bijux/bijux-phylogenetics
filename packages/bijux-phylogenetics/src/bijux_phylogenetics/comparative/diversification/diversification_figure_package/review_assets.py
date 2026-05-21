from __future__ import annotations

from html import escape
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows

from ..models import (
    CladeDiversificationScanReport,
    DiversificationModelComparisonReport,
    LineageThroughTimeReport,
    SamplingFractionReport,
)
from .contracts import (
    DiversificationFigureAudit,
    DiversificationFigureCaptionDraft,
    DiversificationFigureLegendEntry,
)
from .figure_surfaces import format_figure_value


def build_legend_entries() -> list[DiversificationFigureLegendEntry]:
    return [
        DiversificationFigureLegendEntry(
            surface="lineage-through-time",
            label="lineage accumulation over time",
            swatch="#0f766e",
            detail="the lineage-through-time curve keeps the rooted time axis explicit and marks each diversification event directly on the figure",
        ),
        DiversificationFigureLegendEntry(
            surface="clade-outliers",
            label="high diversification outlier",
            swatch="#b91c1c",
            detail="positive z-scores at or above one standard deviation mark clades diversifying faster than the tree-wide baseline",
        ),
        DiversificationFigureLegendEntry(
            surface="clade-outliers",
            label="low diversification outlier",
            swatch="#1d4ed8",
            detail="negative z-scores at or below minus one standard deviation mark clades diversifying slower than the tree-wide baseline",
        ),
        DiversificationFigureLegendEntry(
            surface="clade-outliers",
            label="baseline clade",
            swatch="#94a3b8",
            detail="gray bars mark clades that stay within the baseline diversification band",
        ),
        DiversificationFigureLegendEntry(
            surface="model-comparison",
            label="better-supported model",
            swatch="#ca8a04",
            detail="the highlighted point marks the lowest AIC among the candidate diversification models",
        ),
    ]


def write_legend_table(
    path: Path, entries: list[DiversificationFigureLegendEntry]
) -> Path:
    return write_taxon_rows(
        path,
        columns=["surface", "label", "swatch", "detail"],
        rows=[
            {
                "surface": entry.surface,
                "label": entry.label,
                "swatch": entry.swatch,
                "detail": entry.detail,
            }
            for entry in entries
        ],
    )


def build_audit(
    *,
    lineage_report: LineageThroughTimeReport,
    clade_report: CladeDiversificationScanReport,
    model_report: DiversificationModelComparisonReport,
    sampling_report: SamplingFractionReport | None,
    plotted_ltt_point_count: int,
    plotted_clade_count: int,
    plotted_model_count: int,
    legend_entries: list[DiversificationFigureLegendEntry],
) -> DiversificationFigureAudit:
    lineage_curve_visible = plotted_ltt_point_count >= 2
    clade_outlier_surface_visible = plotted_clade_count > 0
    model_comparison_visible = plotted_model_count >= 2
    legend_complete = {entry.surface for entry in legend_entries} == {
        "lineage-through-time",
        "clade-outliers",
        "model-comparison",
    }
    caption_ready = (
        lineage_curve_visible
        and clade_outlier_surface_visible
        and model_comparison_visible
    )
    sampling_metadata_complete = (
        None if sampling_report is None else sampling_report.complete
    )
    highlighted_outlier_count = len(clade_report.high_diversification_clades) + len(
        clade_report.low_diversification_clades
    )
    publication_ready = (
        lineage_curve_visible
        and clade_outlier_surface_visible
        and model_comparison_visible
        and legend_complete
        and caption_ready
        and (sampling_metadata_complete is not False)
    )
    reviewer_summary = [
        f"lineage-through-time points rendered: {plotted_ltt_point_count}/{len(lineage_report.points)}",
        f"clade diversification rows rendered: {plotted_clade_count}/{len(clade_report.observations)}",
        f"highlighted clade outliers: {highlighted_outlier_count}",
        f"model comparison rows rendered: {plotted_model_count}/{len(model_report.rows)}",
    ]
    limitations: list[str] = []
    if sampling_report is None:
        limitations.append(
            "sampling-aware correction metadata was not supplied, so diversification rates follow the complete-sampling assumption"
        )
    elif not sampling_report.complete:
        limitations.append(
            "sampling metadata is incomplete or invalid, so publication readiness remains blocked for the sampling-aware diversification figures"
        )
    if not lineage_curve_visible:
        limitations.append(
            "the lineage-through-time figure does not retain enough plotted events"
        )
    if not clade_outlier_surface_visible:
        limitations.append(
            "the clade outlier figure does not retain any plotted clade rows"
        )
    if not model_comparison_visible:
        limitations.append(
            "the model-comparison figure does not retain enough candidate models"
        )
    if not legend_complete:
        limitations.append(
            "the figure legend does not cover all rendered diversification surfaces"
        )
    if not limitations:
        limitations.append(
            "the current package keeps lineage accumulation, clade-rate deviations, and model support visible enough for publication-oriented review"
        )
    return DiversificationFigureAudit(
        publication_ready=publication_ready,
        lineage_curve_visible=lineage_curve_visible,
        clade_outlier_surface_visible=clade_outlier_surface_visible,
        model_comparison_visible=model_comparison_visible,
        legend_complete=legend_complete,
        caption_ready=caption_ready,
        sampling_metadata_complete=sampling_metadata_complete,
        plotted_ltt_point_count=plotted_ltt_point_count,
        plotted_clade_count=plotted_clade_count,
        highlighted_outlier_count=highlighted_outlier_count,
        plotted_model_count=plotted_model_count,
        better_model=model_report.better_model,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )


def build_caption_draft(
    *,
    lineage_report: LineageThroughTimeReport,
    clade_report: CladeDiversificationScanReport,
    audit: DiversificationFigureAudit,
) -> DiversificationFigureCaptionDraft:
    return DiversificationFigureCaptionDraft(
        title="Diversification review across lineage accumulation, clade-rate outliers, and model support",
        lead_sentence=(
            f"This package summarizes diversification patterns for one rooted ultrametric tree with {lineage_report.tip_count} tips through three explicit publication surfaces instead of leaving the evidence fragmented across tables and generic report text."
        ),
        lineage_sentence=(
            f"The lineage-through-time figure retains {audit.plotted_ltt_point_count} plotted time points across a crown age of {format_figure_value(lineage_report.root_age)} so the tempo of lineage accumulation remains visible."
        ),
        clade_sentence=(
            f"The clade diversification panel shows {audit.plotted_clade_count} evaluated clades and highlights {audit.highlighted_outlier_count} high- or low-rate outliers relative to the tree-wide baseline."
        ),
        model_sentence=(
            f"The model-comparison panel retains {audit.plotted_model_count} candidate diversification models and highlights `{audit.better_model}` as the better-supported fit by AIC."
        ),
        limitation_sentence=audit.limitations[0],
        caption_ready=audit.caption_ready,
    )


def write_caption(path: Path, draft: DiversificationFigureCaptionDraft) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"# {draft.title}",
                "",
                draft.lead_sentence,
                draft.lineage_sentence,
                draft.clade_sentence,
                draft.model_sentence,
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
    lineage_figure_path: Path,
    clade_figure_path: Path,
    model_figure_path: Path,
    lineage_table_path: Path,
    clade_table_path: Path,
    model_table_path: Path,
    legend_path: Path,
    caption_path: Path,
    methods_summary_path: Path,
    methods_summary_text: str,
    audit: DiversificationFigureAudit,
) -> str:
    figures = {
        "lineage": lineage_figure_path.read_text(encoding="utf-8"),
        "clades": clade_figure_path.read_text(encoding="utf-8"),
        "models": model_figure_path.read_text(encoding="utf-8"),
    }
    audit_rows = "".join(
        "<tr><th>" + escape(label) + "</th><td>" + escape(value) + "</td></tr>"
        for label, value in [
            ("publication_ready", str(audit.publication_ready).lower()),
            ("better_model", audit.better_model),
            ("lineage_curve_visible", str(audit.lineage_curve_visible).lower()),
            (
                "clade_outlier_surface_visible",
                str(audit.clade_outlier_surface_visible).lower(),
            ),
            ("model_comparison_visible", str(audit.model_comparison_visible).lower()),
            (
                "sampling_metadata_complete",
                "n/a"
                if audit.sampling_metadata_complete is None
                else str(audit.sampling_metadata_complete).lower(),
            ),
        ]
    )
    limitation_items = "".join(f"<li>{escape(item)}</li>" for item in audit.limitations)
    return "\n".join(
        [
            "<!DOCTYPE html>",
            '<html lang="en">',
            "<head>",
            '  <meta charset="utf-8">',
            "  <title>Bijux Diversification Figure Review</title>",
            "  <style>",
            "    body { margin: 0; background: linear-gradient(180deg, #f8fafc 0%, #fff7ed 100%); color: #1f2937; font: 16px/1.5 'Iowan Old Style', 'Palatino Linotype', serif; }",
            "    main { max-width: 1240px; margin: 0 auto; padding: 24px; }",
            "    h1, h2 { font-family: 'Avenir Next', 'Segoe UI', sans-serif; }",
            "    h1 { margin-top: 0; color: #7c2d12; }",
            "    .grid { display: grid; grid-template-columns: 1fr; gap: 18px; }",
            "    .panel { background: rgba(255,255,255,0.9); border: 1px solid rgba(124,45,18,0.12); border-radius: 18px; padding: 18px; box-shadow: 0 18px 40px rgba(124,45,18,0.08); }",
            "    .figure-shell svg { width: 100%; height: auto; display: block; }",
            "    table { width: 100%; border-collapse: collapse; }",
            "    th, td { text-align: left; padding: 8px 10px; border-bottom: 1px solid rgba(124,45,18,0.12); vertical-align: top; }",
            "    ul { margin: 8px 0 0 18px; }",
            "    a { color: #9a3412; }",
            "  </style>",
            "</head>",
            "<body>",
            "<main>",
            "  <h1>Bijux Diversification Figure Review</h1>",
            "  <p>Reviewer-facing figure package for lineage-through-time trajectories, clade-rate outlier detection, and diversification model support.</p>",
            '  <section class="panel">',
            "    <h2>Publication Audit</h2>",
            f"    <table><tbody>{audit_rows}</tbody></table>",
            "    <ul>" + limitation_items + "</ul>",
            "  </section>",
            '  <section class="grid" style="margin-top: 20px;">',
            '    <section class="panel"><h2>Lineage-Through-Time Curve</h2><div class="figure-shell">'
            + figures["lineage"]
            + "</div></section>",
            '    <section class="panel"><h2>Clade-Rate Outliers</h2><div class="figure-shell">'
            + figures["clades"]
            + "</div></section>",
            '    <section class="panel"><h2>Model Comparison</h2><div class="figure-shell">'
            + figures["models"]
            + "</div></section>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Methods Summary</h2>",
            f"    <pre>{escape(methods_summary_text)}</pre>",
            "  </section>",
            '  <section class="panel" style="margin-top: 20px;">',
            "    <h2>Linked Artifacts</h2>",
            "    <ul>",
            f'      <li><a href="{escape(lineage_table_path.name)}">{escape(lineage_table_path.name)}</a></li>',
            f'      <li><a href="{escape(clade_table_path.name)}">{escape(clade_table_path.name)}</a></li>',
            f'      <li><a href="{escape(model_table_path.name)}">{escape(model_table_path.name)}</a></li>',
            f'      <li><a href="{escape(legend_path.name)}">{escape(legend_path.name)}</a></li>',
            f'      <li><a href="{escape(caption_path.name)}">{escape(caption_path.name)}</a></li>',
            f'      <li><a href="{escape(methods_summary_path.name)}">{escape(methods_summary_path.name)}</a></li>',
            "    </ul>",
            "  </section>",
            "</main>",
            "</body>",
            "</html>",
        ]
    )

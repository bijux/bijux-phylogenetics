from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
from html import escape
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.render.reproducibility import (
    FigureReproducibilityFilter,
    write_figure_reproducibility_manifest,
)

from ..clades import write_clade_diversification_table
from ..lineage import write_lineage_through_time_table
from ..reporting import (
    build_diversification_method_report,
    write_diversification_methods_summary_text,
    write_diversification_model_comparison_table,
)
from .contracts import (
    DiversificationFigureAudit,
    DiversificationFigureCaptionDraft,
    DiversificationFigureLegendEntry,
    DiversificationFigurePackageResult,
)
from .figure_surfaces import (
    format_figure_value,
    write_clade_outlier_svg,
    write_ltt_svg,
    write_model_comparison_svg,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def _build_legend_entries() -> list[DiversificationFigureLegendEntry]:
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


def _write_legend_table(
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


def _build_audit(
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


def _build_caption_draft(
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


def _write_caption(path: Path, draft: DiversificationFigureCaptionDraft) -> Path:
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


def _build_review_html(
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


def build_diversification_figure_package(
    tree_path: Path,
    *,
    out_dir: Path,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    sampling_column: str | None = None,
    min_tip_count: int = 2,
    model: str = "birth-death",
) -> DiversificationFigurePackageResult:
    """Build a reviewer-facing diversification figure package with explicit publication artifacts."""
    out_dir.mkdir(parents=True, exist_ok=True)
    lineage_figure_path = out_dir / "lineage-through-time.svg"
    clade_figure_path = out_dir / "clade-diversification-outliers.svg"
    model_figure_path = out_dir / "diversification-model-comparison.svg"
    lineage_table_path = out_dir / "lineage-through-time.tsv"
    clade_table_path = out_dir / "clade-diversification-outliers.tsv"
    model_table_path = out_dir / "diversification-model-comparison.tsv"
    legend_path = out_dir / "figure-legend.tsv"
    caption_path = out_dir / "figure-caption.md"
    methods_summary_path = out_dir / "diversification-methods-summary.md"
    review_path = out_dir / "diversification-figure-review.html"
    manifest_path = out_dir / "diversification-figure-package.manifest.json"
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"

    methods_report = build_diversification_method_report(
        tree_path,
        metadata_path=metadata_path,
        taxon_column=taxon_column,
        sampling_column=sampling_column,
        estimate_model=model,
        clade_model=model,
        clade_min_tip_count=min_tip_count,
    )
    lineage_report = methods_report.lineage
    clade_report = methods_report.clade_scan
    model_report = methods_report.model_comparison
    sampling_report = methods_report.sampling_report

    write_lineage_through_time_table(lineage_table_path, lineage_report)
    write_clade_diversification_table(clade_table_path, clade_report)
    write_diversification_model_comparison_table(model_table_path, model_report)
    plotted_ltt_point_count = write_ltt_svg(lineage_figure_path, lineage_report)
    plotted_clade_count = write_clade_outlier_svg(clade_figure_path, clade_report)
    plotted_model_count = write_model_comparison_svg(model_figure_path, model_report)
    legend_entries = _build_legend_entries()
    _write_legend_table(legend_path, legend_entries)
    audit = _build_audit(
        lineage_report=lineage_report,
        clade_report=clade_report,
        model_report=model_report,
        sampling_report=sampling_report,
        plotted_ltt_point_count=plotted_ltt_point_count,
        plotted_clade_count=plotted_clade_count,
        plotted_model_count=plotted_model_count,
        legend_entries=legend_entries,
    )
    caption_draft = _build_caption_draft(
        lineage_report=lineage_report,
        clade_report=clade_report,
        audit=audit,
    )
    _write_caption(caption_path, caption_draft)
    methods_summary = write_diversification_methods_summary_text(
        methods_summary_path,
        methods_report,
    )
    review_path.write_text(
        _build_review_html(
            lineage_figure_path=lineage_figure_path,
            clade_figure_path=clade_figure_path,
            model_figure_path=model_figure_path,
            lineage_table_path=lineage_table_path,
            clade_table_path=clade_table_path,
            model_table_path=model_table_path,
            legend_path=legend_path,
            caption_path=caption_path,
            methods_summary_path=methods_summary_path,
            methods_summary_text=methods_summary.text,
            audit=audit,
        ),
        encoding="utf-8",
    )
    artifact_paths = [
        lineage_figure_path,
        clade_figure_path,
        model_figure_path,
        lineage_table_path,
        clade_table_path,
        model_table_path,
        legend_path,
        caption_path,
        methods_summary_path,
        review_path,
    ]
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="diversification_figure_package",
        input_files=[
            ("tree", tree_path),
            *([("metadata", metadata_path)] if metadata_path is not None else []),
        ],
        generated_figures=[
            ("lineage_through_time", lineage_figure_path),
            ("clade_outliers", clade_figure_path),
            ("model_comparison", model_figure_path),
        ],
        generated_tables=[
            ("lineage_through_time", lineage_table_path),
            ("clade_outliers", clade_table_path),
            ("model_comparison", model_table_path),
        ],
        filters=[
            FigureReproducibilityFilter(
                name="min_tip_count",
                value=str(min_tip_count),
                detail="exclude clades smaller than the configured minimum tip count from outlier review",
            )
        ],
        model={
            "kind": "diversification",
            "name": model,
            "selected_model": model_report.better_model,
            "candidate_models": [row.model for row in model_report.rows],
        },
        settings={
            "taxon_column": taxon_column,
            "sampling_column": sampling_column,
            "metadata_path": None if metadata_path is None else str(metadata_path),
            "tip_count": lineage_report.tip_count,
        },
        linked_artifacts=[
            ("legend", legend_path),
            ("caption", caption_path),
            ("methods_summary", methods_summary_path),
            ("review", review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "diversification_figure_package",
        "input_path": str(tree_path),
        "input_checksum": _checksum(tree_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(reproducibility_manifest_path),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "metadata_path": None if metadata_path is None else str(metadata_path),
            "taxon_column": taxon_column,
            "sampling_column": sampling_column,
            "min_tip_count": min_tip_count,
            "model": model,
            "methods_summary_path": str(methods_summary_path),
        },
        "metrics": {
            "tip_count": lineage_report.tip_count,
            "root_age": lineage_report.root_age,
            "publication_ready": audit.publication_ready,
            "sampling_metadata_complete": audit.sampling_metadata_complete,
            "plotted_ltt_point_count": audit.plotted_ltt_point_count,
            "plotted_clade_count": audit.plotted_clade_count,
            "highlighted_outlier_count": audit.highlighted_outlier_count,
            "plotted_model_count": audit.plotted_model_count,
            "better_model": audit.better_model,
            "methods_summary_warning_count": methods_summary.warning_count,
        },
        "outputs": {
            "methods_summary_path": str(methods_summary_path),
        },
        "lineage_report": _json_ready(asdict(lineage_report)),
        "clade_report": _json_ready(asdict(clade_report)),
        "model_report": _json_ready(asdict(model_report)),
        "sampling_report": None
        if sampling_report is None
        else _json_ready(asdict(sampling_report)),
        "methods_summary": _json_ready(asdict(methods_summary)),
        "audit": _json_ready(asdict(audit)),
    }
    manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return DiversificationFigurePackageResult(
        output_dir=out_dir,
        lineage_figure_path=lineage_figure_path,
        clade_figure_path=clade_figure_path,
        model_figure_path=model_figure_path,
        lineage_table_path=lineage_table_path,
        clade_table_path=clade_table_path,
        model_table_path=model_table_path,
        legend_path=legend_path,
        caption_path=caption_path,
        methods_summary_path=methods_summary_path,
        review_path=review_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        lineage_report=lineage_report,
        clade_report=clade_report,
        model_report=model_report,
        sampling_report=sampling_report,
        methods_report=methods_report,
        methods_summary=methods_summary,
        legend_entries=legend_entries,
        caption_draft=caption_draft,
        audit=audit,
        machine_manifest=machine_manifest,
    )

from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

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
    write_clade_outlier_svg,
    write_ltt_svg,
    write_model_comparison_svg,
)
from .review_assets import (
    build_audit,
    build_caption_draft,
    build_legend_entries,
    build_review_html,
    write_caption,
    write_legend_table,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


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
    legend_entries = build_legend_entries()
    write_legend_table(legend_path, legend_entries)
    audit = build_audit(
        lineage_report=lineage_report,
        clade_report=clade_report,
        model_report=model_report,
        sampling_report=sampling_report,
        plotted_ltt_point_count=plotted_ltt_point_count,
        plotted_clade_count=plotted_clade_count,
        plotted_model_count=plotted_model_count,
        legend_entries=legend_entries,
    )
    caption_draft = build_caption_draft(
        lineage_report=lineage_report,
        clade_report=clade_report,
        audit=audit,
    )
    write_caption(caption_path, caption_draft)
    methods_summary = write_diversification_methods_summary_text(
        methods_summary_path,
        methods_report,
    )
    review_path.write_text(
        build_review_html(
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

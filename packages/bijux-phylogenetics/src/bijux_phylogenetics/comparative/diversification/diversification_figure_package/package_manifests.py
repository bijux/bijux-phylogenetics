from __future__ import annotations

from dataclasses import asdict
from hashlib import sha256
import json
from pathlib import Path

from bijux_phylogenetics.render.reproducibility import (
    FigureReproducibilityFilter,
    write_figure_reproducibility_manifest,
)

from ..models import (
    CladeDiversificationScanReport,
    DiversificationMethodReport,
    DiversificationMethodsSummaryTextResult,
    LineageThroughTimeReport,
    SamplingFractionReport,
)
from .contracts import (
    DiversificationFigureAudit,
    DiversificationFigurePackageResult,
)


def _checksum(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _json_ready(payload: object) -> object:
    return json.loads(json.dumps(payload, default=str))


def write_package_manifests(
    *,
    result: DiversificationFigurePackageResult,
    tree_path: Path,
    metadata_path: Path | None,
    taxon_column: str | None,
    sampling_column: str | None,
    min_tip_count: int,
    model: str,
    lineage_report: LineageThroughTimeReport,
    clade_report: CladeDiversificationScanReport,
    sampling_report: SamplingFractionReport | None,
    methods_report: DiversificationMethodReport,
    methods_summary: DiversificationMethodsSummaryTextResult,
    audit: DiversificationFigureAudit,
) -> tuple[dict[str, object], dict[str, object]]:
    artifact_paths = [
        result.lineage_figure_path,
        result.clade_figure_path,
        result.model_figure_path,
        result.lineage_table_path,
        result.clade_table_path,
        result.model_table_path,
        result.legend_path,
        result.caption_path,
        result.methods_summary_path,
        result.review_path,
    ]
    reproducibility_manifest = write_figure_reproducibility_manifest(
        result.reproducibility_manifest_path,
        report_kind="diversification_figure_package",
        input_files=[
            ("tree", tree_path),
            *([("metadata", metadata_path)] if metadata_path is not None else []),
        ],
        generated_figures=[
            ("lineage_through_time", result.lineage_figure_path),
            ("clade_outliers", result.clade_figure_path),
            ("model_comparison", result.model_figure_path),
        ],
        generated_tables=[
            ("lineage_through_time", result.lineage_table_path),
            ("clade_outliers", result.clade_table_path),
            ("model_comparison", result.model_table_path),
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
            "selected_model": result.model_report.better_model,
            "candidate_models": [row.model for row in result.model_report.rows],
        },
        settings={
            "taxon_column": taxon_column,
            "sampling_column": sampling_column,
            "metadata_path": None if metadata_path is None else str(metadata_path),
            "tip_count": lineage_report.tip_count,
        },
        linked_artifacts=[
            ("legend", result.legend_path),
            ("caption", result.caption_path),
            ("methods_summary", result.methods_summary_path),
            ("review", result.review_path),
        ],
    )
    machine_manifest = {
        "report_kind": "diversification_figure_package",
        "input_path": str(tree_path),
        "input_checksum": _checksum(tree_path),
        "output_paths": [str(path) for path in artifact_paths],
        "output_checksums": {str(path): _checksum(path) for path in artifact_paths},
        "reproducibility_manifest_path": str(result.reproducibility_manifest_path),
        "reproducibility_manifest_checksum": _checksum(
            result.reproducibility_manifest_path
        ),
        "reproducibility_manifest": reproducibility_manifest,
        "settings": {
            "metadata_path": None if metadata_path is None else str(metadata_path),
            "taxon_column": taxon_column,
            "sampling_column": sampling_column,
            "min_tip_count": min_tip_count,
            "model": model,
            "methods_summary_path": str(result.methods_summary_path),
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
            "methods_summary_path": str(result.methods_summary_path),
        },
        "lineage_report": _json_ready(asdict(lineage_report)),
        "clade_report": _json_ready(asdict(clade_report)),
        "model_report": _json_ready(asdict(result.model_report)),
        "sampling_report": None
        if sampling_report is None
        else _json_ready(asdict(sampling_report)),
        "methods_report": _json_ready(asdict(methods_report)),
        "methods_summary": _json_ready(asdict(methods_summary)),
        "audit": _json_ready(asdict(audit)),
    }
    result.manifest_path.write_text(
        json.dumps(machine_manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return reproducibility_manifest, machine_manifest

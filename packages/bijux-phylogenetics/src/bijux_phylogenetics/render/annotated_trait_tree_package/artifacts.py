from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path
from typing import Any

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.render.reproducibility import (
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.render.tree_figure_package import TreeFigurePackageResult

from .contracts import (
    AnnotatedTraitTreeCoverageRow,
    AnnotatedTraitTreePublicationAudit,
    AnnotatedTraitTreeSummaryRow,
)


def sha256(path: Path) -> str:
    """Return the stable SHA-256 digest for one package artifact."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_annotation_coverage_table(
    coverage_path: Path,
    coverage_rows: list[AnnotatedTraitTreeCoverageRow],
) -> None:
    """Write the annotation coverage review table for the package."""
    write_taxon_rows(
        coverage_path,
        columns=[
            "surface",
            "source_kind",
            "required",
            "visible_taxon_count",
            "observed_taxon_count",
            "covered_taxon_count",
            "complete",
            "missing_taxa",
            "extra_taxa",
        ],
        rows=[
            {
                "surface": row.surface,
                "source_kind": row.source_kind,
                "required": row.required,
                "visible_taxon_count": row.visible_taxon_count,
                "observed_taxon_count": row.observed_taxon_count,
                "covered_taxon_count": row.covered_taxon_count,
                "complete": row.complete,
                "missing_taxa": "|".join(row.missing_taxa),
                "extra_taxa": "|".join(row.extra_taxa),
            }
            for row in coverage_rows
        ],
    )


def write_annotation_summary_table(
    summary_path: Path,
    summary_rows: list[AnnotatedTraitTreeSummaryRow],
) -> None:
    """Write the annotation surface summary review table for the package."""
    write_taxon_rows(
        summary_path,
        columns=[
            "surface",
            "source_kind",
            "value_kind",
            "observed_taxon_count",
            "missing_taxon_count",
            "distinct_value_count",
            "minimum_numeric_value",
            "maximum_numeric_value",
            "example_values",
        ],
        rows=[
            {
                "surface": row.surface,
                "source_kind": row.source_kind,
                "value_kind": row.value_kind,
                "observed_taxon_count": row.observed_taxon_count,
                "missing_taxon_count": row.missing_taxon_count,
                "distinct_value_count": row.distinct_value_count,
                "minimum_numeric_value": row.minimum_numeric_value,
                "maximum_numeric_value": row.maximum_numeric_value,
                "example_values": "|".join(row.example_values),
            }
            for row in summary_rows
        ],
    )


def build_package_manifest(
    *,
    title: str,
    tree_path: Path,
    metadata_path: Path | None,
    traits_path: Path | None,
    taxon_column: str | None,
    label_column: str | None,
    categorical_column: str | None,
    continuous_column: str | None,
    metadata_strip_columns: list[str],
    heatmap_columns: list[str],
    layout: str,
    show_support_values: bool,
    figure_package: TreeFigurePackageResult,
    coverage_path: Path,
    summary_path: Path,
    review_path: Path,
    audit: AnnotatedTraitTreePublicationAudit,
    coverage_rows: list[AnnotatedTraitTreeCoverageRow],
    summary_rows: list[AnnotatedTraitTreeSummaryRow],
) -> dict[str, Any]:
    """Build the manifest payload for one annotated trait tree package."""
    return {
        "report_kind": "annotated_trait_tree_package",
        "title": title,
        "tree_path": str(tree_path),
        "metadata_path": None if metadata_path is None else str(metadata_path),
        "traits_path": None if traits_path is None else str(traits_path),
        "taxon_column": taxon_column,
        "label_column": label_column,
        "categorical_column": categorical_column,
        "continuous_column": continuous_column,
        "metadata_strip_columns": metadata_strip_columns,
        "heatmap_columns": heatmap_columns,
        "layout": layout,
        "show_support_values": show_support_values,
        "input_checksums": {
            str(path): sha256(path)
            for path in (tree_path, metadata_path, traits_path)
            if path is not None
        },
        "figure_package_manifest_path": str(figure_package.manifest_path),
        "figure_package_manifest_checksum": sha256(figure_package.manifest_path),
        "coverage_path": str(coverage_path),
        "coverage_checksum": sha256(coverage_path),
        "summary_path": str(summary_path),
        "summary_checksum": sha256(summary_path),
        "review_path": str(review_path),
        "audit": asdict(audit),
        "coverage_rows": [asdict(row) for row in coverage_rows],
        "summary_rows": [asdict(row) for row in summary_rows],
    }


def write_package_reproducibility_manifest(
    reproducibility_manifest_path: Path,
    *,
    title: str,
    tree_path: Path,
    metadata_path: Path | None,
    traits_path: Path | None,
    taxon_column: str | None,
    label_column: str | None,
    categorical_column: str | None,
    continuous_column: str | None,
    metadata_strip_columns: list[str],
    heatmap_columns: list[str],
    layout: str,
    show_support_values: bool,
    figure_package: TreeFigurePackageResult,
    coverage_path: Path,
    summary_path: Path,
) -> dict[str, Any]:
    """Write the reproducibility manifest for one annotated trait tree package."""
    return write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="annotated_trait_tree_package",
        input_files=[
            ("tree", tree_path),
            *([] if metadata_path is None else [("metadata", metadata_path)]),
            *([] if traits_path is None else [("traits", traits_path)]),
        ],
        generated_figures=[("annotated_trait_tree", figure_package.figure_path)],
        generated_tables=[
            ("tree_legend", figure_package.legend_path),
            ("tree_annotations", figure_package.annotations_path),
            ("annotation_coverage", coverage_path),
            ("annotation_surface_summary", summary_path),
        ],
        filters=None,
        model={
            "kind": "none",
            "name": None,
            "detail": "the annotated trait tree package overlays supplied labels and trait metadata without fitting a new statistical model",
        },
        settings={
            "title": title,
            "taxon_column": taxon_column,
            "label_column": label_column,
            "categorical_column": categorical_column,
            "continuous_column": continuous_column,
            "metadata_strip_columns": metadata_strip_columns,
            "heatmap_columns": heatmap_columns,
            "layout": layout,
            "show_support_values": show_support_values,
        },
        linked_artifacts=[
            ("tree_caption", figure_package.caption_path),
            ("tree_figure_manifest", figure_package.manifest_path),
            (
                "tree_figure_reproducibility_manifest",
                figure_package.reproducibility_manifest_path,
            ),
        ],
    )


def write_package_manifest(manifest_path: Path, manifest: dict[str, Any]) -> None:
    """Persist the annotated trait tree package manifest."""
    manifest_path.write_text(
        json.dumps(manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

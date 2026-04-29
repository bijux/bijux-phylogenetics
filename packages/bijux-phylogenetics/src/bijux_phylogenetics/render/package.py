from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.svg import AnnotationStrip, TreeRenderResult, render_tree_svg


@dataclass(slots=True)
class TreeFigurePackageResult:
    output_dir: Path
    figure_path: Path
    manifest_path: Path
    caption_path: Path
    annotations_path: Path
    render: TreeRenderResult


def build_tree_figure_package(
    tree_path: Path,
    *,
    out_dir: Path,
    title: str = "Bijux Tree Figure",
    labels: dict[str, str] | None = None,
    layout: str = "phylogram",
    show_support_values: bool = False,
    categorical_traits: dict[str, str] | None = None,
    continuous_traits: dict[str, float] | None = None,
    metadata_strips: list[AnnotationStrip] | None = None,
    heatmap_columns: list[AnnotationStrip] | None = None,
    collapsed_clades: list[str] | None = None,
) -> TreeFigurePackageResult:
    """Build a standalone figure package for publication and review workflows."""
    out_dir.mkdir(parents=True, exist_ok=True)
    figure_path = out_dir / "figure.svg"
    manifest_path = out_dir / "figure-manifest.json"
    caption_path = out_dir / "figure-caption.md"
    annotations_path = out_dir / "tip-annotations.tsv"

    render = render_tree_svg(
        tree_path,
        out_path=figure_path,
        labels=labels,
        layout=layout,
        show_support_values=show_support_values,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
        collapsed_clades=collapsed_clades,
    )

    tree = load_tree(tree_path)
    labels = labels or {}
    categorical_traits = categorical_traits or {}
    continuous_traits = continuous_traits or {}
    metadata_strips = metadata_strips or []
    heatmap_columns = heatmap_columns or []

    annotation_columns = ["taxon", "label", "categorical_trait", "continuous_trait"]
    annotation_columns.extend(strip.name for strip in metadata_strips)
    annotation_columns.extend(column.name for column in heatmap_columns)
    annotation_rows = []
    for taxon in tree.tip_names:
        row = {
            "taxon": taxon,
            "label": labels.get(taxon, taxon),
            "categorical_trait": categorical_traits.get(taxon, ""),
            "continuous_trait": str(continuous_traits[taxon]) if taxon in continuous_traits else "",
        }
        for strip in metadata_strips:
            row[strip.name] = strip.values.get(taxon, "")
        for column in heatmap_columns:
            row[column.name] = column.values.get(taxon, "")
        annotation_rows.append(row)
    write_taxon_rows(annotations_path, columns=annotation_columns, rows=annotation_rows)

    manifest = {
        "title": title,
        "tree_path": str(tree_path),
        "figure_path": str(figure_path),
        "caption_path": str(caption_path),
        "annotations_path": str(annotations_path),
        "layout": layout,
        "collapsed_clades": sorted(set(collapsed_clades or [])),
        "render": asdict(render),
    }
    manifest_path.write_text(json.dumps(manifest, default=str, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    caption = (
        f"# {title}\n\n"
        f"- Layout: `{layout}`\n"
        f"- Tips rendered: `{render.tip_count}`\n"
        f"- Support labels rendered: `{render.rendered_support_count}`\n"
        f"- Metadata strips: `{render.rendered_metadata_strip_count}`\n"
        f"- Heatmap columns: `{render.rendered_heatmap_column_count}`\n"
    )
    caption_path.write_text(caption, encoding="utf-8")

    return TreeFigurePackageResult(
        output_dir=out_dir,
        figure_path=figure_path,
        manifest_path=manifest_path,
        caption_path=caption_path,
        annotations_path=annotations_path,
        render=render,
    )

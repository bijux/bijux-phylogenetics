from __future__ import annotations

from dataclasses import asdict
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.reproducibility import (
    FigureReproducibilityFilter,
    write_figure_reproducibility_manifest,
)
from bijux_phylogenetics.render.tree_svg import (
    AnnotationStrip,
    SupportLabelRenderAudit,
    audit_support_label_rendering,
    render_tree_svg,
)

from .audits import (
    build_collapsed_clade_summaries,
    build_surface_coverage,
    build_table_consistency,
    visible_tip_taxa,
)
from .contracts import TreeFigureAuditReport, TreeFigurePackageResult
from .legends import build_legend_audit, build_legend_entries
from .review import build_caption_draft, build_legibility_audit


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
    reproducibility_manifest_path = out_dir / "figure-reproducibility.manifest.json"
    caption_path = out_dir / "figure-caption.md"
    legend_path = out_dir / "figure-legend.tsv"
    annotations_path = out_dir / "tip-annotations.tsv"
    collapsed_clades = collapsed_clades or []
    collapsed_clade_names = {name for name in collapsed_clades if name}
    support_audit = (
        audit_support_label_rendering(tree_path)
        if show_support_values
        else SupportLabelRenderAudit(
            validated=False,
            labels_by_node={},
            warnings=[],
        )
    )

    tree = load_tree(tree_path)
    labels = labels or {taxon: taxon for taxon in tree.tip_names}

    render = render_tree_svg(
        tree_path,
        out_path=figure_path,
        labels=labels,
        layout=layout,
        show_support_values=show_support_values and support_audit.validated,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
        collapsed_clades=collapsed_clades,
        validated_support_labels=support_audit.labels_by_node,
        support_validation_warnings=support_audit.warnings,
    )

    categorical_traits = categorical_traits or {}
    continuous_traits = continuous_traits or {}
    metadata_strips = metadata_strips or []
    heatmap_columns = heatmap_columns or []
    visible_taxa = visible_tip_taxa(tree.root, collapsed_clade_names)

    annotation_columns = ["taxon", "label", "categorical_trait", "continuous_trait"]
    annotation_columns.extend(strip.name for strip in metadata_strips)
    annotation_columns.extend(column.name for column in heatmap_columns)
    annotation_rows = []
    for taxon in tree.tip_names:
        row = {
            "taxon": taxon,
            "label": labels.get(taxon, taxon),
            "categorical_trait": categorical_traits.get(taxon, ""),
            "continuous_trait": str(continuous_traits[taxon])
            if taxon in continuous_traits
            else "",
        }
        for strip in metadata_strips:
            row[strip.name] = strip.values.get(taxon, "")
        for column in heatmap_columns:
            row[column.name] = column.values.get(taxon, "")
        annotation_rows.append(row)
    write_taxon_rows(annotations_path, columns=annotation_columns, rows=annotation_rows)

    annotation_coverage = [
        build_surface_coverage(
            surface="labels",
            visible_taxa=visible_taxa,
            observed_taxa=set(labels),
        )
    ]
    if categorical_traits:
        annotation_coverage.append(
            build_surface_coverage(
                surface="categorical traits",
                visible_taxa=visible_taxa,
                observed_taxa={
                    taxon for taxon, value in categorical_traits.items() if value
                },
            )
        )
    if continuous_traits:
        annotation_coverage.append(
            build_surface_coverage(
                surface="continuous traits",
                visible_taxa=visible_taxa,
                observed_taxa=set(continuous_traits),
            )
        )
    for strip in metadata_strips:
        annotation_coverage.append(
            build_surface_coverage(
                surface=f"metadata strip: {strip.name}",
                visible_taxa=visible_taxa,
                observed_taxa={taxon for taxon, value in strip.values.items() if value},
            )
        )
    for column in heatmap_columns:
        annotation_coverage.append(
            build_surface_coverage(
                surface=f"heatmap column: {column.name}",
                visible_taxa=visible_taxa,
                observed_taxa={
                    taxon for taxon, value in column.values.items() if value
                },
            )
        )

    collapsed_summaries = build_collapsed_clade_summaries(
        tree=tree,
        collapsed=collapsed_clade_names,
        metadata_strips=metadata_strips,
    )
    legend_audit = build_legend_audit(
        render=render,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
    )
    table_consistency = build_table_consistency(
        visible_taxa=tree.tip_names,
        labels=labels,
        annotation_rows=annotation_rows,
    )
    scale_bar_valid = (not render.has_scale_bar) or (
        render.scale_bar_length is not None
        and render.max_branch_distance is not None
        and 0 < render.scale_bar_length <= (render.max_branch_distance / 3)
    )
    scale_bar_note = (
        "scale bar withheld because the selected layout is not branch-length proportional"
        if not render.has_scale_bar
        else (
            f"scale bar length {render.scale_bar_length} audited against maximum branch distance {render.max_branch_distance}"
            if scale_bar_valid
            else "scale bar length did not satisfy proportional rendering safeguards"
        )
    )
    reviewer_summary = [
        f"rendered {render.visible_tip_count} visible tips from {render.tip_count} source taxa",
        "support labels were rendered only after support-scale audit"
        if show_support_values and support_audit.validated
        else "support labels were withheld or omitted when no validated support audit was available",
        "all figure annotation surfaces align with the exported annotation table"
        if table_consistency.consistent
        and all(row.aligned for row in annotation_coverage)
        else "one or more figure annotation surfaces require reviewer attention",
    ]
    limitations = []
    if collapsed_summaries:
        limitations.append(
            "collapsed clades summarize hidden descendants and should not be read as individual tip-level evidence"
        )
    for row in annotation_coverage:
        if not row.aligned:
            limitations.append(
                f"{row.surface} is missing values for one or more visible taxa"
            )
    if not scale_bar_valid:
        limitations.append(
            "branch-length scale metadata needs manual review before publication use"
        )
    limitations.extend(support_audit.warnings)
    audit = TreeFigureAuditReport(
        support_audit=support_audit,
        annotation_coverage=annotation_coverage,
        collapsed_clades=collapsed_summaries,
        legend_audit=legend_audit,
        table_consistency=table_consistency,
        scale_bar_valid=scale_bar_valid,
        scale_bar_note=scale_bar_note,
        reviewer_summary=reviewer_summary,
        limitations=limitations,
    )
    visible_label_rows = [
        labels.get(taxon, taxon) for taxon in visible_taxa if labels.get(taxon, taxon)
    ]
    legend_entries = build_legend_entries(
        render=render,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
    )
    legibility_audit = build_legibility_audit(
        render=render,
        visible_labels=visible_label_rows,
        has_annotation_columns=bool(
            categorical_traits
            or continuous_traits
            or metadata_strips
            or heatmap_columns
        ),
    )
    caption_draft = build_caption_draft(
        title=title,
        render=render,
        audit=audit,
        legend_entries=legend_entries,
    )

    write_taxon_rows(
        legend_path,
        columns=["surface", "label", "swatch", "detail"],
        rows=[asdict(entry) for entry in legend_entries],
    )

    manifest = {
        "title": title,
        "tree_path": str(tree_path),
        "input_checksums": {str(tree_path): _sha256(tree_path)},
        "figure_path": str(figure_path),
        "figure_checksum": _sha256(figure_path),
        "caption_path": str(caption_path),
        "legend_path": str(legend_path),
        "legend_checksum": _sha256(legend_path),
        "annotations_path": str(annotations_path),
        "annotations_checksum": _sha256(annotations_path),
        "layout": layout,
        "collapsed_clades": sorted(set(collapsed_clades or [])),
        "render": asdict(render),
        "audit": asdict(audit),
        "legibility_audit": asdict(legibility_audit),
        "legend_entries": [asdict(entry) for entry in legend_entries],
        "caption_draft": asdict(caption_draft),
    }
    caption = (
        f"# {title}\n\n"
        "## Draft Caption\n\n"
        f"{caption_draft.lead_sentence} "
        f"{caption_draft.support_sentence} "
        f"{caption_draft.scale_bar_sentence} "
        f"{caption_draft.legend_sentence} "
        f"{caption_draft.limitation_sentence}\n\n"
        "## Figure Specifications\n\n"
        f"- Layout: `{layout}`\n"
        f"- Tips rendered: `{render.visible_tip_count}` visible from `{render.tip_count}` source taxa\n"
        f"- Support labels rendered: `{render.rendered_support_count}`\n"
        f"- Metadata strips: `{render.rendered_metadata_strip_count}`\n"
        f"- Heatmap columns: `{render.rendered_heatmap_column_count}`\n"
        f"- Legend entries: `{len(legend_entries)}`\n"
        f"- Caption ready: `{caption_draft.caption_ready}`\n"
        f"- Scale bar audit: `{audit.scale_bar_note}`\n"
        f"- Legibility audit: `legible={legibility_audit.legible}`\n\n"
        "## Reviewer Summary\n\n"
        + "\n".join(f"- {line}" for line in audit.reviewer_summary)
        + "\n\n## Legend\n\n"
        + "\n".join(
            f"- {entry.surface}: {entry.label} ({entry.swatch}) — {entry.detail}"
            for entry in legend_entries
        )
        + "\n\n## Legibility\n\n"
        + "\n".join(
            [
                f"- tip label font size: `{legibility_audit.tip_label_font_size_px}px`",
                f"- vertical tip spacing: `{legibility_audit.vertical_tip_spacing_px}px`",
                f"- longest visible label length: `{legibility_audit.longest_visible_label_length}`",
                f"- estimated longest label width: `{legibility_audit.estimated_longest_label_width_px}px`",
                f"- available label lane: `{legibility_audit.available_label_lane_px}px`",
            ]
        )
        + "\n\n## Limitations\n\n"
        + "\n".join(
            f"- {line}"
            for line in (
                [*audit.limitations, *legibility_audit.warnings]
                or ["no additional limitations recorded"]
            )
        )
        + "\n"
    )
    caption_path.write_text(caption, encoding="utf-8")
    reproducibility_manifest = write_figure_reproducibility_manifest(
        reproducibility_manifest_path,
        report_kind="tree_figure_package",
        input_files=[("tree", tree_path)],
        generated_figures=[("tree_figure", figure_path)],
        generated_tables=[
            ("figure_legend", legend_path),
            ("tip_annotations", annotations_path),
        ],
        filters=[]
        if not collapsed_clades
        else [
            FigureReproducibilityFilter(
                name="collapsed_clades",
                value="|".join(sorted(set(collapsed_clades))),
                detail="these named clades were intentionally collapsed before the publication figure was rendered",
            )
        ],
        model={
            "kind": "none",
            "name": None,
            "detail": "the tree figure package renders the supplied tree directly and does not fit a new statistical model",
        },
        settings={
            "title": title,
            "layout": layout,
            "show_support_values": show_support_values,
            "collapsed_clades": sorted(set(collapsed_clades)),
            "categorical_trait_count": 0
            if categorical_traits is None
            else len(categorical_traits),
            "continuous_trait_count": 0
            if continuous_traits is None
            else len(continuous_traits),
            "metadata_strip_count": 0
            if metadata_strips is None
            else len(metadata_strips),
            "heatmap_column_count": 0
            if heatmap_columns is None
            else len(heatmap_columns),
        },
        linked_artifacts=[
            ("figure_caption", caption_path),
        ],
    )
    manifest["reproducibility_manifest_path"] = str(reproducibility_manifest_path)
    manifest["reproducibility_manifest_checksum"] = _sha256(
        reproducibility_manifest_path
    )
    manifest["reproducibility_manifest"] = reproducibility_manifest
    manifest_path.write_text(
        json.dumps(manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return TreeFigurePackageResult(
        output_dir=out_dir,
        figure_path=figure_path,
        manifest_path=manifest_path,
        reproducibility_manifest_path=reproducibility_manifest_path,
        caption_path=caption_path,
        legend_path=legend_path,
        annotations_path=annotations_path,
        render=render,
        audit=audit,
        legend_entries=legend_entries,
        legibility_audit=legibility_audit,
        caption_draft=caption_draft,
    )


__all__ = [
    "build_tree_figure_package",
]

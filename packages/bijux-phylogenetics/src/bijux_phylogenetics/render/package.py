from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.core.metadata import write_taxon_rows
from bijux_phylogenetics.core.tree import TreeNode
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.svg import (
    AnnotationStrip,
    SupportLabelRenderAudit,
    TreeRenderResult,
    audit_support_label_rendering,
    render_tree_svg,
)


@dataclass(frozen=True, slots=True)
class FigureAnnotationCoverage:
    """Coverage and drift audit for one figure annotation surface."""

    surface: str
    aligned: bool
    covered_taxa: int
    missing_taxa: list[str]
    extra_taxa: list[str]


@dataclass(frozen=True, slots=True)
class FigureCollapsedCladeSummary:
    """Reviewer-facing summary for one collapsed clade."""

    clade_name: str
    descendant_count: int
    descendant_taxa: list[str]
    metadata_summaries: list[str]


@dataclass(frozen=True, slots=True)
class FigureLegendAudit:
    """Legend completeness check for rendered figure surfaces."""

    complete: bool
    entries: list[str]
    missing_entries: list[str]
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class FigureTableConsistencyReport:
    """Consistency check between rendered annotations and tabular exports."""

    consistent: bool
    missing_from_table: list[str]
    extra_in_table: list[str]
    label_mismatches: list[str]


@dataclass(slots=True)
class TreeFigureAuditReport:
    """Combined reviewer-facing audit for a rendered figure package."""

    support_audit: SupportLabelRenderAudit
    annotation_coverage: list[FigureAnnotationCoverage]
    collapsed_clades: list[FigureCollapsedCladeSummary]
    legend_audit: FigureLegendAudit
    table_consistency: FigureTableConsistencyReport
    scale_bar_valid: bool
    scale_bar_note: str
    reviewer_summary: list[str]
    limitations: list[str]


@dataclass(slots=True)
class TreeFigurePackageResult:
    output_dir: Path
    figure_path: Path
    manifest_path: Path
    caption_path: Path
    annotations_path: Path
    render: TreeRenderResult
    audit: TreeFigureAuditReport


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _visible_tip_taxa(node: TreeNode, collapsed: set[str]) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    if node.name is not None and node.name in collapsed:
        return [node.name]
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_visible_tip_taxa(child, collapsed))
    return taxa


def _iter_collapsed_nodes(node: TreeNode, collapsed: set[str]) -> list[TreeNode]:
    if node.is_leaf():
        return []
    matches: list[TreeNode] = []
    if node.name is not None and node.name in collapsed:
        matches.append(node)
    else:
        for child in node.children:
            matches.extend(_iter_collapsed_nodes(child, collapsed))
    return matches


def _surface_coverage(
    *,
    surface: str,
    visible_taxa: list[str],
    observed_taxa: set[str],
) -> FigureAnnotationCoverage:
    visible_taxa_set = set(visible_taxa)
    missing_taxa = sorted(taxon for taxon in visible_taxa if taxon not in observed_taxa)
    extra_taxa = sorted(observed_taxa - visible_taxa_set)
    return FigureAnnotationCoverage(
        surface=surface,
        aligned=not missing_taxa and not extra_taxa,
        covered_taxa=len(visible_taxa_set & observed_taxa),
        missing_taxa=missing_taxa,
        extra_taxa=extra_taxa,
    )


def _legend_audit(
    *,
    render: TreeRenderResult,
    categorical_traits: dict[str, str],
    continuous_traits: dict[str, float],
    metadata_strips: list[AnnotationStrip],
    heatmap_columns: list[AnnotationStrip],
) -> FigureLegendAudit:
    entries: list[str] = []
    if categorical_traits:
        entries.append("categorical trait legend")
    if continuous_traits:
        entries.append("continuous trait gradient")
    if metadata_strips:
        entries.extend(f"metadata strip: {strip.name}" for strip in metadata_strips)
    if heatmap_columns:
        entries.extend(f"heatmap column: {column.name}" for column in heatmap_columns)
    if render.rendered_support_count:
        entries.append("support label audit")
    if render.has_scale_bar:
        entries.append("branch-length scale bar")

    missing_entries: list[str] = []
    warnings: list[str] = []
    if render.rendered_support_count and not render.support_labels_validated:
        missing_entries.append("validated support label note")
        warnings.extend(render.support_validation_warnings)
    complete = not missing_entries
    return FigureLegendAudit(
        complete=complete,
        entries=entries,
        missing_entries=missing_entries,
        warnings=warnings,
    )


def _collapsed_clade_summaries(
    *,
    tree,
    collapsed: set[str],
    metadata_strips: list[AnnotationStrip],
) -> list[FigureCollapsedCladeSummary]:
    summaries: list[FigureCollapsedCladeSummary] = []
    for node in _iter_collapsed_nodes(tree.root, collapsed):
        taxa = _descendant_taxa(node)
        metadata_summaries: list[str] = []
        for strip in metadata_strips:
            counts: dict[str, int] = {}
            for taxon in taxa:
                value = strip.values.get(taxon)
                if value:
                    counts[value] = counts.get(value, 0) + 1
            if counts:
                summary = ", ".join(
                    f"{value}={counts[value]}" for value in sorted(counts)
                )
                metadata_summaries.append(f"{strip.name}: {summary}")
        summaries.append(
            FigureCollapsedCladeSummary(
                clade_name=node.name or "collapsed clade",
                descendant_count=len(taxa),
                descendant_taxa=taxa,
                metadata_summaries=metadata_summaries,
            )
        )
    return summaries


def _table_consistency(
    *,
    visible_taxa: list[str],
    labels: dict[str, str],
    annotation_rows: list[dict[str, str]],
) -> FigureTableConsistencyReport:
    rows_by_taxon = {row["taxon"]: row for row in annotation_rows}
    visible_taxa_set = set(visible_taxa)
    missing_from_table = sorted(
        taxon for taxon in visible_taxa if taxon not in rows_by_taxon
    )
    extra_in_table = sorted(set(rows_by_taxon) - visible_taxa_set)
    label_mismatches: list[str] = []
    for taxon in visible_taxa:
        row = rows_by_taxon.get(taxon)
        if row is None:
            continue
        expected_label = labels.get(taxon, taxon)
        if row.get("label") != expected_label:
            label_mismatches.append(
                f"{taxon}: expected '{expected_label}' but found '{row.get('label', '')}'"
            )
    return FigureTableConsistencyReport(
        consistent=not missing_from_table
        and not extra_in_table
        and not label_mismatches,
        missing_from_table=missing_from_table,
        extra_in_table=extra_in_table,
        label_mismatches=label_mismatches,
    )


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

    tree = load_tree(tree_path)
    labels = labels or {}
    categorical_traits = categorical_traits or {}
    continuous_traits = continuous_traits or {}
    metadata_strips = metadata_strips or []
    heatmap_columns = heatmap_columns or []
    visible_taxa = _visible_tip_taxa(tree.root, collapsed_clade_names)

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
        _surface_coverage(
            surface="labels",
            visible_taxa=visible_taxa,
            observed_taxa=set(labels),
        )
    ]
    if categorical_traits:
        annotation_coverage.append(
            _surface_coverage(
                surface="categorical traits",
                visible_taxa=visible_taxa,
                observed_taxa={
                    taxon for taxon, value in categorical_traits.items() if value
                },
            )
        )
    if continuous_traits:
        annotation_coverage.append(
            _surface_coverage(
                surface="continuous traits",
                visible_taxa=visible_taxa,
                observed_taxa=set(continuous_traits),
            )
        )
    for strip in metadata_strips:
        annotation_coverage.append(
            _surface_coverage(
                surface=f"metadata strip: {strip.name}",
                visible_taxa=visible_taxa,
                observed_taxa={taxon for taxon, value in strip.values.items() if value},
            )
        )
    for column in heatmap_columns:
        annotation_coverage.append(
            _surface_coverage(
                surface=f"heatmap column: {column.name}",
                visible_taxa=visible_taxa,
                observed_taxa={
                    taxon for taxon, value in column.values.items() if value
                },
            )
        )

    collapsed_summaries = _collapsed_clade_summaries(
        tree=tree,
        collapsed=collapsed_clade_names,
        metadata_strips=metadata_strips,
    )
    legend_audit = _legend_audit(
        render=render,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
    )
    table_consistency = _table_consistency(
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

    manifest = {
        "title": title,
        "tree_path": str(tree_path),
        "input_checksums": {str(tree_path): _sha256(tree_path)},
        "figure_path": str(figure_path),
        "figure_checksum": _sha256(figure_path),
        "caption_path": str(caption_path),
        "annotations_path": str(annotations_path),
        "annotations_checksum": _sha256(annotations_path),
        "layout": layout,
        "collapsed_clades": sorted(set(collapsed_clades or [])),
        "render": asdict(render),
        "audit": asdict(audit),
    }
    manifest_path.write_text(
        json.dumps(manifest, default=str, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    caption = (
        f"# {title}\n\n"
        f"- Layout: `{layout}`\n"
        f"- Tips rendered: `{render.visible_tip_count}` visible from `{render.tip_count}` source taxa\n"
        f"- Support labels rendered: `{render.rendered_support_count}`\n"
        f"- Metadata strips: `{render.rendered_metadata_strip_count}`\n"
        f"- Heatmap columns: `{render.rendered_heatmap_column_count}`\n"
        f"- Scale bar audit: `{audit.scale_bar_note}`\n\n"
        "## Reviewer Summary\n\n"
        + "\n".join(f"- {line}" for line in audit.reviewer_summary)
        + "\n\n## Legend\n\n"
        + "\n".join(f"- {entry}" for entry in audit.legend_audit.entries)
        + "\n\n## Limitations\n\n"
        + "\n".join(
            f"- {line}"
            for line in (audit.limitations or ["no additional limitations recorded"])
        )
        + "\n"
    )
    caption_path.write_text(caption, encoding="utf-8")

    return TreeFigurePackageResult(
        output_dir=out_dir,
        figure_path=figure_path,
        manifest_path=manifest_path,
        caption_path=caption_path,
        annotations_path=annotations_path,
        render=render,
        audit=audit,
    )

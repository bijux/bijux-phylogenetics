from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
from pathlib import Path

from bijux_phylogenetics.datasets.study_inputs import write_taxon_rows
from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.render.reproducibility import (
    FigureReproducibilityFilter,
    write_figure_reproducibility_manifest,
)
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
class FigureLegendEntry:
    """One explicit legend entry emitted alongside a publication figure."""

    surface: str
    label: str
    swatch: str
    detail: str


@dataclass(frozen=True, slots=True)
class FigureLegibilityAudit:
    """Heuristic legibility review for one rendered publication tree figure."""

    legible: bool
    tip_label_font_size_px: int
    vertical_tip_spacing_px: int
    longest_visible_label_length: int
    estimated_longest_label_width_px: float
    available_label_lane_px: int
    warnings: list[str]


@dataclass(frozen=True, slots=True)
class FigureCaptionDraft:
    """Structured caption draft for one publication-oriented tree figure."""

    title: str
    lead_sentence: str
    support_sentence: str
    scale_bar_sentence: str
    legend_sentence: str
    limitation_sentence: str
    caption_ready: bool


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
    reproducibility_manifest_path: Path
    caption_path: Path
    legend_path: Path
    annotations_path: Path
    render: TreeRenderResult
    audit: TreeFigureAuditReport
    legend_entries: list[FigureLegendEntry]
    legibility_audit: FigureLegibilityAudit
    caption_draft: FigureCaptionDraft


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


def _categorical_legend_entries(
    *,
    surface: str,
    values: dict[str, str],
) -> list[FigureLegendEntry]:
    categories = sorted({value for value in values.values() if value})
    if not categories:
        return []
    colors = {
        category: color
        for category, color in zip(
            categories,
            (
                "#0f766e",
                "#1d4ed8",
                "#c2410c",
                "#7c3aed",
                "#b91c1c",
                "#047857",
                "#a16207",
                "#0f172a",
            ),
            strict=False,
        )
    }
    return [
        FigureLegendEntry(
            surface=surface,
            label=category,
            swatch=colors[category],
            detail=f"{surface} category rendered directly on the vector figure",
        )
        for category in categories
    ]


def _build_legend_entries(
    *,
    render: TreeRenderResult,
    categorical_traits: dict[str, str],
    continuous_traits: dict[str, float],
    metadata_strips: list[AnnotationStrip],
    heatmap_columns: list[AnnotationStrip],
) -> list[FigureLegendEntry]:
    entries: list[FigureLegendEntry] = []
    if render.has_scale_bar and render.scale_bar_length is not None:
        entries.append(
            FigureLegendEntry(
                surface="branch-length",
                label="scale bar",
                swatch="#0f172a",
                detail=f"scale bar represents branch length {render.scale_bar_length}",
            )
        )
    if render.rendered_support_count:
        entries.append(
            FigureLegendEntry(
                surface="support",
                label="validated support labels",
                swatch="#0f766e",
                detail=f"{render.rendered_support_count} support labels were rendered after audit",
            )
        )
    entries.extend(
        _categorical_legend_entries(
            surface="categorical trait",
            values=categorical_traits,
        )
    )
    if continuous_traits:
        minimum = min(continuous_traits.values())
        maximum = max(continuous_traits.values())
        entries.append(
            FigureLegendEntry(
                surface="continuous trait",
                label="gradient range",
                swatch=f"{minimum:.6g}..{maximum:.6g}",
                detail="continuous trait values are mapped to a low-to-high gradient",
            )
        )
    for strip in metadata_strips:
        entries.extend(
            _categorical_legend_entries(
                surface=f"metadata strip: {strip.name}",
                values=strip.values,
            )
        )
    for column in heatmap_columns:
        observed_values = [value for value in column.values.values() if value]
        if not observed_values:
            continue
        if all(_is_number(value) for value in observed_values):
            numeric_values = [float(value) for value in observed_values]
            entries.append(
                FigureLegendEntry(
                    surface=f"heatmap column: {column.name}",
                    label="numeric gradient",
                    swatch=f"{min(numeric_values):.6g}..{max(numeric_values):.6g}",
                    detail="heatmap column uses a continuous low-to-high color gradient",
                )
            )
        else:
            entries.extend(
                _categorical_legend_entries(
                    surface=f"heatmap column: {column.name}",
                    values=column.values,
                )
            )
    return entries


def _is_number(value: str) -> bool:
    try:
        float(value)
    except ValueError:
        return False
    return True


def _build_legibility_audit(
    *,
    render: TreeRenderResult,
    visible_labels: list[str],
    has_annotation_columns: bool,
) -> FigureLegibilityAudit:
    tip_label_font_size_px = 16
    vertical_tip_spacing_px = 56 if render.layout != "circular" else 20
    longest_visible_label_length = max((len(label) for label in visible_labels), default=0)
    estimated_longest_label_width_px = round(
        longest_visible_label_length * tip_label_font_size_px * 0.58,
        2,
    )
    available_label_lane_px = (
        172 if has_annotation_columns and render.layout != "circular" else 420
    )
    warnings: list[str] = []
    if estimated_longest_label_width_px > available_label_lane_px:
        warnings.append(
            "the longest visible label likely exceeds the reserved label lane and should be shortened or moved to metadata"
        )
    if render.layout != "circular" and vertical_tip_spacing_px < tip_label_font_size_px * 2:
        warnings.append(
            "vertical tip spacing fell below the publication legibility threshold"
        )
    if render.visible_tip_count > 80:
        warnings.append(
            "the visible tip count is high enough that clade collapsing or panel subdivision should be considered before journal submission"
        )
    return FigureLegibilityAudit(
        legible=not warnings,
        tip_label_font_size_px=tip_label_font_size_px,
        vertical_tip_spacing_px=vertical_tip_spacing_px,
        longest_visible_label_length=longest_visible_label_length,
        estimated_longest_label_width_px=estimated_longest_label_width_px,
        available_label_lane_px=available_label_lane_px,
        warnings=warnings,
    )


def _build_caption_draft(
    *,
    title: str,
    render: TreeRenderResult,
    audit: TreeFigureAuditReport,
    legend_entries: list[FigureLegendEntry],
) -> FigureCaptionDraft:
    lead_sentence = (
        f"{title} shows {render.visible_tip_count} rendered taxa from a source tree with {render.tip_count} total tips using a {render.layout} layout."
    )
    support_sentence = (
        f"Validated branch support labels are shown for {render.rendered_support_count} internal branches."
        if render.rendered_support_count
        else "Branch support labels were omitted because no validated support surface was available."
    )
    scale_bar_sentence = (
        f"Branch lengths are scaled directly on the figure using a scale bar of {render.scale_bar_length}."
        if render.has_scale_bar and render.scale_bar_length is not None
        else "A branch-length scale bar is not shown because the selected layout is not branch-length proportional."
    )
    legend_sentence = (
        f"The figure legend contains {len(legend_entries)} explicit entries covering rendered trait, metadata, support, and scale surfaces."
    )
    limitation_sentence = (
        "Reviewer-facing audits did not record additional publication limitations."
        if not (audit.limitations or audit.legend_audit.warnings)
        else "Publication review still requires attention to: "
        + "; ".join([*audit.limitations, *audit.legend_audit.warnings])
        + "."
    )
    return FigureCaptionDraft(
        title=title,
        lead_sentence=lead_sentence,
        support_sentence=support_sentence,
        scale_bar_sentence=scale_bar_sentence,
        legend_sentence=legend_sentence,
        limitation_sentence=limitation_sentence,
        caption_ready=audit.scale_bar_valid
        and audit.legend_audit.complete
        and audit.table_consistency.consistent,
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
    visible_label_rows = [
        labels.get(taxon, taxon) for taxon in visible_taxa if labels.get(taxon, taxon)
    ]
    legend_entries = _build_legend_entries(
        render=render,
        categorical_traits=categorical_traits,
        continuous_traits=continuous_traits,
        metadata_strips=metadata_strips,
        heatmap_columns=heatmap_columns,
    )
    legibility_audit = _build_legibility_audit(
        render=render,
        visible_labels=visible_label_rows,
        has_annotation_columns=bool(
            categorical_traits or continuous_traits or metadata_strips or heatmap_columns
        ),
    )
    caption_draft = _build_caption_draft(
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

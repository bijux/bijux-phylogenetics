from __future__ import annotations

from bijux_phylogenetics.phylo.topology.tree import TreeNode
from bijux_phylogenetics.render.tree_svg import AnnotationStrip

from .contracts import (
    FigureAnnotationCoverage,
    FigureCollapsedCladeSummary,
    FigureTableConsistencyReport,
)


def descendant_taxa(node: TreeNode) -> list[str]:
    """Collect all descendant taxa for one tree node."""
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(descendant_taxa(child))
    return sorted(taxa)


def visible_tip_taxa(node: TreeNode, collapsed: set[str]) -> list[str]:
    """Collect the visible taxa after collapsed clades are applied."""
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    if node.name is not None and node.name in collapsed:
        return [node.name]
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(visible_tip_taxa(child, collapsed))
    return taxa


def iter_collapsed_nodes(node: TreeNode, collapsed: set[str]) -> list[TreeNode]:
    """Find the named collapsed clades in the rendered tree."""
    if node.is_leaf():
        return []
    matches: list[TreeNode] = []
    if node.name is not None and node.name in collapsed:
        matches.append(node)
    else:
        for child in node.children:
            matches.extend(iter_collapsed_nodes(child, collapsed))
    return matches


def build_surface_coverage(
    *,
    surface: str,
    visible_taxa: list[str],
    observed_taxa: set[str],
) -> FigureAnnotationCoverage:
    """Summarize which visible taxa are covered by one annotation surface."""
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


def build_collapsed_clade_summaries(
    *,
    tree,
    collapsed: set[str],
    metadata_strips: list[AnnotationStrip],
) -> list[FigureCollapsedCladeSummary]:
    """Summarize the descendants hidden behind each collapsed clade."""
    summaries: list[FigureCollapsedCladeSummary] = []
    for node in iter_collapsed_nodes(tree.root, collapsed):
        taxa = descendant_taxa(node)
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


def build_table_consistency(
    *,
    visible_taxa: list[str],
    labels: dict[str, str],
    annotation_rows: list[dict[str, str]],
) -> FigureTableConsistencyReport:
    """Check the exported annotation table against the rendered labels."""
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


__all__ = [
    "build_collapsed_clade_summaries",
    "build_surface_coverage",
    "build_table_consistency",
    "descendant_taxa",
    "iter_collapsed_nodes",
    "visible_tip_taxa",
]

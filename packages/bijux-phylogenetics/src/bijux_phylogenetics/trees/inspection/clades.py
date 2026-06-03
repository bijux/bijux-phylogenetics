from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.BaseTree import Clade
from Bio.Phylo.BaseTree import Tree as BioTree

from bijux_phylogenetics.datasets.study_inputs import load_taxon_table
from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.iqtree_support import (
    parse_iqtree_branch_support_label,
    support_fraction,
)
from bijux_phylogenetics.io.trees import detect_tree_format
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


@dataclass(frozen=True, slots=True)
class CladeMetadataObservation:
    column: str
    matched_taxa: int
    missing_taxa: list[str]
    values_by_taxon: list[str]
    distinct_values: list[str]


@dataclass(frozen=True, slots=True)
class CladeTableRow:
    source_path: Path
    tree_index: int | None
    node_kind: str
    clade_id: str
    node_label: str | None
    taxon_count: int
    taxa: list[str]
    support: float | None
    support_fraction: float | None
    branch_length: float | None
    root_depth: float | None
    descendant_tip_depth_min: float | None
    descendant_tip_depth_max: float | None
    node_age: float | None
    metadata: list[CladeMetadataObservation]


@dataclass(slots=True)
class CladeTableReport:
    path: Path
    source_format: str
    tree_count: int
    metadata_path: Path | None
    taxon_column: str | None
    metadata_columns: list[str]
    rows: list[CladeTableRow]


@dataclass(slots=True)
class _MetadataContext:
    path: Path
    taxon_column: str
    columns: list[str]
    rows_by_taxon: dict[str, dict[str, str]]


def _resolve_metadata_context(
    metadata_path: Path | None,
    *,
    taxon_column: str | None,
    metadata_columns: list[str] | None,
) -> tuple[_MetadataContext | None, list[str]]:
    if metadata_path is None:
        return None, []
    table = load_taxon_table(metadata_path, taxon_column=taxon_column)
    columns = metadata_columns or [
        column for column in table.columns if column != table.taxon_column
    ]
    for column in columns:
        if column not in table.columns:
            raise ValueError(
                f"metadata column '{column}' was not found in {metadata_path}"
            )
        if column == table.taxon_column:
            raise ValueError(
                f"metadata column '{column}' is the taxon key and cannot be summarized as a trait column"
            )
    return (
        _MetadataContext(
            path=table.path,
            taxon_column=table.taxon_column,
            columns=columns,
            rows_by_taxon={row[table.taxon_column]: row for row in table.rows},
        ),
        columns,
    )


def _node_label(clade: Clade) -> str | None:
    if clade.name is not None:
        return clade.name
    if clade.confidence is not None:
        return format(float(clade.confidence), ".15g")
    return None


def _parse_support(value: object) -> float | None:
    if value is None:
        return None
    if isinstance(value, str):
        parsed = parse_iqtree_branch_support_label(value)
        if parsed is not None:
            return (
                parsed.ufboot_support
                if parsed.ufboot_support is not None
                else parsed.sh_alrt_support
            )
    try:
        return float(value)
    except ValueError:
        return None


def _metadata_observations(
    taxa: list[str],
    context: _MetadataContext | None,
) -> list[CladeMetadataObservation]:
    if context is None:
        return []
    observations: list[CladeMetadataObservation] = []
    for column in context.columns:
        values_by_taxon: list[str] = []
        missing_taxa: list[str] = []
        distinct_values: set[str] = set()
        for taxon in taxa:
            row = context.rows_by_taxon.get(taxon)
            value = "" if row is None else row.get(column, "").strip()
            if not value:
                missing_taxa.append(taxon)
                continue
            values_by_taxon.append(f"{taxon}={value}")
            distinct_values.add(value)
        observations.append(
            CladeMetadataObservation(
                column=column,
                matched_taxa=len(values_by_taxon),
                missing_taxa=missing_taxa,
                values_by_taxon=values_by_taxon,
                distinct_values=sorted(distinct_values),
            )
        )
    return observations


def _ultrametric_crown_height(
    tree: PhyloTree, *, tolerance: float = 1e-9
) -> float | None:
    lengths = tree.root_to_tip_lengths()
    if not lengths or any(length is None for length in lengths):
        return None
    numeric_lengths = [float(length) for length in lengths]
    if max(numeric_lengths) - min(numeric_lengths) > tolerance:
        return None
    return max(numeric_lengths)


def _extract_clade_rows(
    source_path: Path,
    tree: PhyloTree,
    bio_tree: BioTree,
    *,
    tree_index: int | None,
    metadata_context: _MetadataContext | None,
) -> list[CladeTableRow]:
    crown_height = _ultrametric_crown_height(tree)
    rows: list[CladeTableRow] = []

    def visit(
        node: TreeNode,
        clade: Clade,
        *,
        root_depth: float | None,
        is_root: bool,
    ) -> tuple[list[str], float | None, float | None]:
        if node.is_leaf():
            taxa = [node.name] if node.name is not None else []
            descendant_tip_depth_min = 0.0
            descendant_tip_depth_max = 0.0
        else:
            taxa = []
            child_mins: list[float] = []
            child_maxes: list[float] = []
            complete_descendant_depths = True
            for child_node, child_clade in zip(
                node.children,
                clade.clades,
                strict=True,
            ):
                child_root_depth = (
                    None
                    if root_depth is None or child_node.branch_length is None
                    else root_depth + child_node.branch_length
                )
                child_taxa, child_min, child_max = visit(
                    child_node,
                    child_clade,
                    root_depth=child_root_depth,
                    is_root=False,
                )
                taxa.extend(child_taxa)
                if (
                    child_node.branch_length is None
                    or child_min is None
                    or child_max is None
                ):
                    complete_descendant_depths = False
                    continue
                child_mins.append(child_node.branch_length + child_min)
                child_maxes.append(child_node.branch_length + child_max)
            if complete_descendant_depths and child_mins and child_maxes:
                descendant_tip_depth_min = min(child_mins)
                descendant_tip_depth_max = max(child_maxes)
            else:
                descendant_tip_depth_min = None
                descendant_tip_depth_max = None
        rounded_tip_depth_min = (
            None
            if descendant_tip_depth_min is None
            else round(descendant_tip_depth_min, 15)
        )
        rounded_tip_depth_max = (
            None
            if descendant_tip_depth_max is None
            else round(descendant_tip_depth_max, 15)
        )
        ordered_taxa = sorted(taxa)
        node_kind = "root" if is_root else ("tip" if node.is_leaf() else "internal")
        support = None
        if clade.clades:
            support = _parse_support(
                clade.confidence if clade.confidence is not None else clade.name
            )
        rows.append(
            CladeTableRow(
                source_path=source_path,
                tree_index=tree_index,
                node_kind=node_kind,
                clade_id="|".join(ordered_taxa),
                node_label=_node_label(clade),
                taxon_count=len(ordered_taxa),
                taxa=ordered_taxa,
                support=support,
                support_fraction=support_fraction(support),
                branch_length=node.branch_length,
                root_depth=None if root_depth is None else round(root_depth, 15),
                descendant_tip_depth_min=rounded_tip_depth_min,
                descendant_tip_depth_max=rounded_tip_depth_max,
                node_age=(
                    None
                    if crown_height is None or root_depth is None
                    else round(crown_height - root_depth, 15)
                ),
                metadata=_metadata_observations(ordered_taxa, metadata_context),
            )
        )
        return ordered_taxa, rounded_tip_depth_min, rounded_tip_depth_max

    visit(tree.root, bio_tree.root, root_depth=0.0, is_root=True)
    return rows


def _load_tree_with_biophylo(path: Path) -> tuple[str, BioTree, PhyloTree]:
    source_format = detect_tree_format(path)
    bio_tree = Phylo.read(path, source_format)
    return (
        source_format,
        bio_tree,
        tree_from_biophylo(bio_tree, source_format=source_format),
    )


def extract_tree_clades(
    path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    metadata_columns: list[str] | None = None,
) -> CladeTableReport:
    """Extract one row per clade from a single tree."""
    metadata_context, selected_columns = _resolve_metadata_context(
        metadata_path,
        taxon_column=taxon_column,
        metadata_columns=metadata_columns,
    )
    source_format, bio_tree, tree = _load_tree_with_biophylo(path)
    return CladeTableReport(
        path=path,
        source_format=source_format,
        tree_count=1,
        metadata_path=None if metadata_context is None else metadata_context.path,
        taxon_column=None
        if metadata_context is None
        else metadata_context.taxon_column,
        metadata_columns=selected_columns,
        rows=_extract_clade_rows(
            path,
            tree,
            bio_tree,
            tree_index=None,
            metadata_context=metadata_context,
        ),
    )


def extract_tree_set_clades(
    path: Path,
    *,
    metadata_path: Path | None = None,
    taxon_column: str | None = None,
    metadata_columns: list[str] | None = None,
) -> CladeTableReport:
    """Extract one row per clade from each tree in a tree set."""
    metadata_context, selected_columns = _resolve_metadata_context(
        metadata_path,
        taxon_column=taxon_column,
        metadata_columns=metadata_columns,
    )
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    trees = [
        tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees
    ]
    reference_taxa = sorted(trees[0].tip_names)
    for tree in trees[1:]:
        if sorted(tree.tip_names) != reference_taxa:
            raise InvalidAlignmentError(
                "clade extraction across a tree set requires identical taxon sets for every tree"
            )
    rows: list[CladeTableRow] = []
    for index, (bio_tree, tree) in enumerate(
        zip(bio_trees, trees, strict=True), start=1
    ):
        rows.extend(
            _extract_clade_rows(
                path,
                tree,
                bio_tree,
                tree_index=index,
                metadata_context=metadata_context,
            )
        )
    return CladeTableReport(
        path=path,
        source_format=source_format,
        tree_count=len(trees),
        metadata_path=None if metadata_context is None else metadata_context.path,
        taxon_column=None
        if metadata_context is None
        else metadata_context.taxon_column,
        metadata_columns=selected_columns,
        rows=rows,
    )


def write_clade_table(path: Path, report: CladeTableReport) -> Path:
    """Write extracted clade rows as a TSV table."""
    fieldnames = [
        "source_path",
        "tree_index",
        "node_kind",
        "clade_id",
        "node_label",
        "taxon_count",
        "taxa",
        "support",
        "support_fraction",
        "branch_length",
        "root_depth",
        "descendant_tip_depth_min",
        "descendant_tip_depth_max",
        "node_age",
    ]
    for column in report.metadata_columns:
        fieldnames.extend(
            [
                f"{column}_values",
                f"{column}_distinct_values",
                f"{column}_missing_taxa",
            ]
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in report.rows:
            metadata_by_column = {item.column: item for item in row.metadata}
            rendered_row = {
                "source_path": str(row.source_path),
                "tree_index": "" if row.tree_index is None else row.tree_index,
                "node_kind": row.node_kind,
                "clade_id": row.clade_id,
                "node_label": "" if row.node_label is None else row.node_label,
                "taxon_count": row.taxon_count,
                "taxa": "|".join(row.taxa),
                "support": "" if row.support is None else row.support,
                "support_fraction": (
                    "" if row.support_fraction is None else row.support_fraction
                ),
                "branch_length": "" if row.branch_length is None else row.branch_length,
                "root_depth": "" if row.root_depth is None else row.root_depth,
                "descendant_tip_depth_min": (
                    ""
                    if row.descendant_tip_depth_min is None
                    else row.descendant_tip_depth_min
                ),
                "descendant_tip_depth_max": (
                    ""
                    if row.descendant_tip_depth_max is None
                    else row.descendant_tip_depth_max
                ),
                "node_age": "" if row.node_age is None else row.node_age,
            }
            for column in report.metadata_columns:
                metadata = metadata_by_column[column]
                rendered_row[f"{column}_values"] = "|".join(metadata.values_by_taxon)
                rendered_row[f"{column}_distinct_values"] = "|".join(
                    metadata.distinct_values
                )
                rendered_row[f"{column}_missing_taxa"] = "|".join(metadata.missing_taxa)
            writer.writerow(rendered_row)
    return path

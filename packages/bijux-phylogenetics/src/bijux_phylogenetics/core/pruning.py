from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.alignment import AlignmentRecord
from bijux_phylogenetics.core.metadata import load_taxon_table
from bijux_phylogenetics.core.tree import PhyloTree, TreeNode
from bijux_phylogenetics.errors import MetadataJoinError
from bijux_phylogenetics.io.fasta import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class TreePruningReport:
    """Explicit record of tree pruning against a taxon table."""

    tree_path: Path
    keep_from_path: Path
    taxon_column: str
    original_tip_count: int
    kept_taxa: list[str]
    removed_taxa: list[str]


@dataclass(slots=True)
class AlignmentPruningReport:
    """Explicit record of alignment pruning against a tree tip set."""

    alignment_path: Path
    tree_path: Path
    original_sequence_count: int
    kept_ids: list[str]
    removed_ids: list[str]


def _merge_branch_lengths(left: float | None, right: float | None) -> float | None:
    if left is None:
        return right
    if right is None:
        return left
    return left + right


def _collapse_unary(node: TreeNode) -> TreeNode:
    if len(node.children) != 1:
        return node
    only_child = _collapse_unary(node.children[0])
    return TreeNode(
        name=only_child.name,
        branch_length=_merge_branch_lengths(node.branch_length, only_child.branch_length),
        children=only_child.children,
    )


def _prune_node(node: TreeNode, keep_taxa: set[str]) -> TreeNode | None:
    if node.is_leaf():
        if node.name is None or node.name not in keep_taxa:
            return None
        return TreeNode(name=node.name, branch_length=node.branch_length, children=[])

    children = [child for child in (_prune_node(child, keep_taxa) for child in node.children) if child is not None]
    if not children:
        return None
    pruned_node = TreeNode(name=node.name, branch_length=node.branch_length, children=children)
    if len(children) == 1:
        return _collapse_unary(pruned_node)
    return pruned_node


def _prune_tree_against_taxa(tree: PhyloTree, keep_taxa: set[str]) -> tuple[PhyloTree, list[str], list[str]]:
    retained_tips = sorted(name for name in tree.tip_names if name in keep_taxa)
    removed_tips = sorted(name for name in tree.tip_names if name not in keep_taxa)

    pruned_root = _prune_node(tree.root, keep_taxa)
    if pruned_root is None or not retained_tips:
        raise MetadataJoinError("no overlapping taxa remain after pruning request")

    pruned_root.branch_length = None
    pruned_tree = PhyloTree(root=pruned_root, source_format=tree.source_format)
    return pruned_tree, retained_tips, removed_tips


def prune_tree_to_taxa(
    tree_path: Path,
    keep_from_path: Path,
    *,
    taxon_column: str | None = None,
) -> tuple[PhyloTree, TreePruningReport]:
    """Prune a tree to the taxa present in a metadata or traits table."""
    tree = load_tree(tree_path)
    keep_table = load_taxon_table(keep_from_path, taxon_column=taxon_column)
    pruned_tree, retained_tips, removed_tips = _prune_tree_against_taxa(tree, set(keep_table.taxa))
    return pruned_tree, TreePruningReport(
        tree_path=tree_path,
        keep_from_path=keep_from_path,
        taxon_column=keep_table.taxon_column,
        original_tip_count=tree.tip_count,
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
    )


def write_pruned_taxa(path: Path, removed_taxa: list[str]) -> Path:
    """Write the removed taxa for an explicit pruning run."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["taxon"]
    lines.extend(removed_taxa)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def prune_alignment_to_tree(
    alignment_path: Path,
    tree_path: Path,
) -> tuple[list[AlignmentRecord], AlignmentPruningReport]:
    """Prune an alignment to the taxa present in a tree."""
    records = load_fasta_alignment(alignment_path)
    tree = load_tree(tree_path)
    tree_taxa = set(tree.tip_names)
    kept_records = [record for record in records if record.identifier in tree_taxa]
    if not kept_records:
        raise MetadataJoinError("no overlapping taxa remain after alignment pruning request")

    kept_ids = [record.identifier for record in kept_records]
    removed_ids = [record.identifier for record in records if record.identifier not in tree_taxa]
    return kept_records, AlignmentPruningReport(
        alignment_path=alignment_path,
        tree_path=tree_path,
        original_sequence_count=len(records),
        kept_ids=kept_ids,
        removed_ids=removed_ids,
    )


def prune_tree_to_alignment(
    tree_path: Path,
    alignment_path: Path,
) -> tuple[PhyloTree, TreePruningReport]:
    """Prune a tree to the identifiers present in an alignment."""
    tree = load_tree(tree_path)
    records = load_fasta_alignment(alignment_path)
    pruned_tree, retained_tips, removed_tips = _prune_tree_against_taxa(
        tree,
        {record.identifier for record in records},
    )
    return pruned_tree, TreePruningReport(
        tree_path=tree_path,
        keep_from_path=alignment_path,
        taxon_column="identifier",
        original_tip_count=tree.tip_count,
        kept_taxa=retained_tips,
        removed_taxa=removed_tips,
    )

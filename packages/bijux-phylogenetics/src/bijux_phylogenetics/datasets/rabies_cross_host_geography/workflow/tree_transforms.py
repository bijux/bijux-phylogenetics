from __future__ import annotations

from pathlib import Path

from Bio import Phylo

from bijux_phylogenetics.io.biopython import tree_from_biophylo, tree_to_biophylo
from bijux_phylogenetics.io.newick import write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.simulation import write_tree_set
from bijux_phylogenetics.trees import CladeTableReport, CladeTableRow

from ..models import RabiesComparativeBranchRepair


def _write_comparative_tree(
    rooted_tree_path: Path,
    *,
    out_path: Path,
    branch_length_floor: float,
) -> tuple[Path, list[RabiesComparativeBranchRepair]]:
    tree = load_tree(rooted_tree_path)
    repairs = _apply_branch_length_floor(tree.root, floor=branch_length_floor)
    write_newick(out_path, tree)
    return out_path, repairs


def _write_comparative_tree_set(
    rooted_tree_set_path: Path,
    *,
    out_path: Path,
    reference_tree_path: Path,
    branch_length_floor: float,
) -> Path:
    reference_length_lookup = _build_branch_length_lookup(
        load_tree(reference_tree_path)
    )
    rooted_trees = _load_tree_set_trees(rooted_tree_set_path)
    adjusted_trees = []
    for tree in rooted_trees:
        _overlay_branch_lengths_from_reference(
            tree.root,
            reference_length_lookup=reference_length_lookup,
            floor=branch_length_floor,
        )
        _apply_branch_length_floor(tree.root, floor=branch_length_floor)
        adjusted_trees.append(tree)
    return write_tree_set(out_path, adjusted_trees)


def _write_rooted_tree_set_on_outgroup(
    tree_set_path: Path,
    *,
    out_path: Path,
    outgroup_taxa: list[str],
) -> Path:
    rooted_trees = [
        _root_tree_on_outgroup_from_tree(tree, outgroup_taxa=outgroup_taxa)
        for tree in _load_tree_set_trees(tree_set_path)
    ]
    return write_tree_set(out_path, rooted_trees)


def _load_tree_set_trees(path: Path) -> list[PhyloTree]:
    source_format = "newick"
    return [
        tree_from_biophylo(tree, source_format=source_format)
        for tree in Phylo.parse(path, source_format)
    ]


def _root_tree_on_outgroup_from_tree(
    tree: PhyloTree, *, outgroup_taxa: list[str]
) -> PhyloTree:
    biophylo_tree = tree_to_biophylo(tree)
    matched = [
        next(biophylo_tree.find_clades(name=taxon), None) for taxon in outgroup_taxa
    ]
    if not any(clade is not None for clade in matched):
        raise ValueError(
            "none of the requested outgroup taxa were found while rooting a tree set"
        )
    biophylo_tree.root_with_outgroup(*[clade for clade in matched if clade is not None])
    return tree_from_biophylo(biophylo_tree, source_format="newick")


def _build_branch_length_lookup(tree: PhyloTree) -> dict[str, float]:
    lookup: dict[str, float] = {}
    for node in tree.iter_nodes():
        if node is tree.root or node.branch_length is None:
            continue
        lookup[_node_signature(node)] = node.branch_length
    return lookup


def _overlay_branch_lengths_from_reference(
    node: TreeNode,
    *,
    reference_length_lookup: dict[str, float],
    floor: float,
) -> None:
    for child in node.children:
        signature = _node_signature(child)
        child.branch_length = reference_length_lookup.get(signature, floor)
        _overlay_branch_lengths_from_reference(
            child,
            reference_length_lookup=reference_length_lookup,
            floor=floor,
        )


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed>"
    return "|".join(sorted(_descendant_taxa(node)))


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return taxa


def _apply_branch_length_floor(
    root: TreeNode,
    *,
    floor: float,
) -> list[RabiesComparativeBranchRepair]:
    repairs: list[RabiesComparativeBranchRepair] = []

    def visit(node: TreeNode, *, is_root: bool) -> None:
        if not is_root and node.branch_length is not None and node.branch_length <= 0.0:
            repairs.append(
                RabiesComparativeBranchRepair(
                    node_label=node.name or "<internal>",
                    original_branch_length=float(node.branch_length),
                    repaired_branch_length=floor,
                    reason=(
                        "comparative branch-length methods require strictly positive "
                        "nonroot branch lengths"
                    ),
                )
            )
            node.branch_length = floor
        for child in node.children:
            visit(child, is_root=False)

    visit(root, is_root=True)
    return repairs


def _stabilize_clade_report(
    report: CladeTableReport,
    *,
    stable_source_path: Path,
) -> CladeTableReport:
    return CladeTableReport(
        path=stable_source_path,
        source_format=report.source_format,
        tree_count=report.tree_count,
        metadata_path=report.metadata_path,
        taxon_column=report.taxon_column,
        metadata_columns=list(report.metadata_columns),
        rows=[
            CladeTableRow(
                source_path=stable_source_path,
                tree_index=row.tree_index,
                node_kind=row.node_kind,
                clade_id=row.clade_id,
                node_label=row.node_label,
                taxon_count=row.taxon_count,
                taxa=list(row.taxa),
                support=row.support,
                support_fraction=row.support_fraction,
                branch_length=row.branch_length,
                root_depth=row.root_depth,
                descendant_tip_depth_min=row.descendant_tip_depth_min,
                descendant_tip_depth_max=row.descendant_tip_depth_max,
                node_age=row.node_age,
                metadata=list(row.metadata),
            )
            for row in report.rows
        ],
    )

from __future__ import annotations

from collections import Counter
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.node_identity import (
    ape_node_id_for_node,
    build_ape_internal_node_map,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .models import (
    CladeExtractionReport,
    SubtreeExtractionReport,
    TreeMonophylyReport,
    TreeMrcaReport,
)
from .transformation import _clone_node, _descendant_taxa, _summarize_transformation


def _find_named_nodes(node: TreeNode, *, clade_name: str) -> list[TreeNode]:
    matches: list[TreeNode] = []
    if node.name == clade_name:
        matches.append(node)
    for child in node.children:
        matches.extend(_find_named_nodes(child, clade_name=clade_name))
    return matches


def _interpreted_rooted_state(tree: PhyloTree) -> bool:
    if tree.rooted is True:
        return True
    return len(tree.root.children) == 2


def _build_subtree(
    node: TreeNode, *, source_format: str, rooted: bool | None
) -> PhyloTree:
    subtree_root = _clone_node(node)
    subtree_root.branch_length = None
    return PhyloTree(root=subtree_root, source_format=source_format, rooted=rooted)


def _root_to_tip_paths(tree: PhyloTree) -> dict[str, list[TreeNode]]:
    paths: dict[str, list[TreeNode]] = {}

    def visit(node: TreeNode, path: list[TreeNode]) -> None:
        next_path = [*path, node]
        if node.is_leaf():
            if node.name is not None:
                paths[node.name] = next_path
            return
        for child in node.children:
            visit(child, next_path)

    visit(tree.root, [])
    return paths


def _common_prefix_length(left: list[TreeNode], right: list[TreeNode]) -> int:
    length = 0
    for left_node, right_node in zip(left, right, strict=False):
        if left_node is not right_node:
            break
        length += 1
    return length


def _mrca_node_from_taxa(tree: PhyloTree, taxa: list[str]) -> TreeNode:
    paths = _root_to_tip_paths(tree)
    reference_path = paths[taxa[0]]
    prefix_length = len(reference_path)
    for taxon in taxa[1:]:
        prefix_length = min(
            prefix_length,
            _common_prefix_length(reference_path, paths[taxon]),
        )
    return reference_path[prefix_length - 1]


def _monophyly_report_from_node(
    tree: PhyloTree,
    *,
    tree_path: Path,
    requested_taxa: list[str],
    unique_requested_taxa: list[str],
    duplicate_requested_taxa: list[str],
    missing_requested_taxa: list[str],
    present_requested_taxa: list[str],
    reroot: bool,
    monophyletic: bool,
    complementary_clade_used: bool,
    matched_node: TreeNode | None,
) -> TreeMonophylyReport:
    matched_taxa = [] if matched_node is None else _descendant_taxa(matched_node)
    matched_extra_taxa = sorted(set(matched_taxa) - set(present_requested_taxa))
    matched_node_id = None
    matched_node_name = None
    is_root = None
    if matched_node is not None:
        matched_node_id = ape_node_id_for_node(tree, matched_node)
        matched_node_name = matched_node.name
        is_root = matched_node is tree.root
    return TreeMonophylyReport(
        tree_path=tree_path,
        requested_taxa=requested_taxa,
        unique_requested_taxa=unique_requested_taxa,
        duplicate_requested_taxa=duplicate_requested_taxa,
        missing_requested_taxa=missing_requested_taxa,
        present_requested_taxa=present_requested_taxa,
        reroot=reroot,
        rooted=_interpreted_rooted_state(tree),
        monophyletic=monophyletic,
        complementary_clade_used=complementary_clade_used,
        matched_node_id=matched_node_id,
        matched_node_name=matched_node_name,
        matched_taxa=matched_taxa,
        matched_extra_taxa=matched_extra_taxa,
        matched_tip_count=len(matched_taxa),
        is_root=is_root,
    )


def _extract_subtree_report(
    tree: PhyloTree,
    subtree: PhyloTree,
    *,
    tree_path: Path,
    selector_kind: str,
    requested_node_id: int | None,
    matched_node_id: int,
    requested_taxa: list[str],
    matched_node_name: str | None,
    expected_taxa: list[str],
) -> SubtreeExtractionReport:
    observed_taxa = sorted(subtree.tip_names)
    summary = _summarize_transformation(
        tree, subtree, transformation="extract-tree-clade"
    )
    return SubtreeExtractionReport(
        tree_path=tree_path,
        selector_kind=selector_kind,
        requested_node_id=requested_node_id,
        matched_node_id=matched_node_id,
        requested_taxa=requested_taxa,
        matched_node_name=matched_node_name,
        tip_count=subtree.tip_count,
        taxa=observed_taxa,
        retained_all_requested_descendants=observed_taxa == expected_taxa,
        missing_requested_descendants=sorted(set(expected_taxa) - set(observed_taxa)),
        unexpected_retained_taxa=sorted(set(observed_taxa) - set(expected_taxa)),
        summary=summary,
    )


def extract_named_clade(
    tree_path: Path,
    *,
    clade_name: str,
) -> tuple[PhyloTree, CladeExtractionReport]:
    """Extract a named clade as a standalone subtree."""
    tree = load_tree(tree_path)
    matches = _find_named_nodes(tree.root, clade_name=clade_name)
    if not matches:
        raise ValueError(f"clade '{clade_name}' was not found in {tree_path}")
    if len(matches) > 1:
        raise ValueError(f"clade '{clade_name}' is ambiguous in {tree_path}")

    subtree = _build_subtree(
        matches[0],
        source_format=tree.source_format,
        rooted=tree.rooted,
    )
    expected_taxa = sorted(_descendant_taxa(matches[0]))
    observed_taxa = sorted(subtree.tip_names)
    summary = _summarize_transformation(
        tree, subtree, transformation="extract-named-clade"
    )
    return subtree, CladeExtractionReport(
        tree_path=tree_path,
        clade_name=clade_name,
        tip_count=subtree.tip_count,
        taxa=observed_taxa,
        retained_all_requested_descendants=observed_taxa == expected_taxa,
        missing_requested_descendants=sorted(set(expected_taxa) - set(observed_taxa)),
        unexpected_retained_taxa=sorted(set(observed_taxa) - set(expected_taxa)),
        summary=summary,
    )


def extract_tree_clade_by_node_id(
    tree_path: Path,
    *,
    node_id: int,
) -> tuple[PhyloTree, SubtreeExtractionReport]:
    """Extract a subtree using ape-style internal node numbering."""
    tree = load_tree(tree_path)
    if node_id <= tree.tip_count:
        raise ValueError("node number must be greater than the number of tips")
    node_map = build_ape_internal_node_map(tree)
    if node_id not in node_map:
        raise IndexError(f"node id {node_id} is out of bounds for {tree_path}")

    source_node = node_map[node_id]
    subtree = _build_subtree(
        source_node,
        source_format=tree.source_format,
        rooted=True,
    )
    expected_taxa = sorted(_descendant_taxa(source_node))
    return subtree, _extract_subtree_report(
        tree,
        subtree,
        tree_path=tree_path,
        selector_kind="node-id",
        requested_node_id=node_id,
        matched_node_id=node_id,
        requested_taxa=[],
        matched_node_name=source_node.name,
        expected_taxa=expected_taxa,
    )


def extract_tree_clade_by_descendant_taxa(
    tree_path: Path,
    *,
    descendant_taxa: list[str],
) -> tuple[PhyloTree, SubtreeExtractionReport]:
    """Extract a subtree whose descendants match the requested taxa exactly."""
    tree = load_tree(tree_path)
    requested_taxa = sorted(set(descendant_taxa))
    if len(requested_taxa) < 2:
        raise ValueError("descendant taxa must contain at least two distinct taxa")
    tree_taxa = set(tree.tip_names)
    missing_taxa = sorted(set(requested_taxa) - tree_taxa)
    if missing_taxa:
        raise ValueError(
            "requested descendant taxa are not present in the tree: "
            + ", ".join(missing_taxa)
        )

    matches: list[tuple[int, TreeNode]] = []
    for matched_node_id, node in build_ape_internal_node_map(tree).items():
        if _descendant_taxa(node) == requested_taxa:
            matches.append((matched_node_id, node))
    if not matches:
        raise ValueError("requested descendant taxa do not define an internal clade")
    if len(matches) > 1:
        raise ValueError("requested descendant taxa are ambiguous in this tree")

    matched_node_id, source_node = matches[0]
    subtree = _build_subtree(
        source_node,
        source_format=tree.source_format,
        rooted=True,
    )
    return subtree, _extract_subtree_report(
        tree,
        subtree,
        tree_path=tree_path,
        selector_kind="descendant-taxa",
        requested_node_id=None,
        matched_node_id=matched_node_id,
        requested_taxa=requested_taxa,
        matched_node_name=source_node.name,
        expected_taxa=requested_taxa,
    )


def find_tree_mrca(
    tree_path: Path,
    *,
    taxa: list[str],
) -> TreeMrcaReport:
    """Resolve the MRCA of two or more taxa using ape-style internal node ids."""
    tree = load_tree(tree_path)
    requested_taxa = sorted(taxa)
    unique_requested_taxa = sorted(set(requested_taxa))
    duplicate_requested_taxa = sorted(
        taxon for taxon, count in Counter(requested_taxa).items() if count > 1
    )
    missing_requested_taxa = sorted(set(unique_requested_taxa) - set(tree.tip_names))
    if missing_requested_taxa:
        raise ValueError(
            "requested taxa are not present in the tree: "
            + ", ".join(missing_requested_taxa)
        )
    if len(unique_requested_taxa) < 2:
        raise ValueError("mrca requires at least two distinct taxa")

    matched_node = _mrca_node_from_taxa(tree, unique_requested_taxa)
    matched_taxa = _descendant_taxa(matched_node)
    matched_extra_taxa = sorted(set(matched_taxa) - set(unique_requested_taxa))
    matched_node_id = next(
        node_id
        for node_id, node in build_ape_internal_node_map(tree).items()
        if node is matched_node
    )
    return TreeMrcaReport(
        tree_path=tree_path,
        requested_taxa=requested_taxa,
        unique_requested_taxa=unique_requested_taxa,
        duplicate_requested_taxa=duplicate_requested_taxa,
        matched_node_id=matched_node_id,
        matched_node_name=matched_node.name,
        matched_taxa=matched_taxa,
        matched_extra_taxa=matched_extra_taxa,
        matched_tip_count=len(matched_taxa),
        is_root=matched_node is tree.root,
        rooted=tree.rooted,
    )


def assess_tree_monophyly(
    tree_path: Path,
    *,
    taxa: list[str],
    reroot: bool = False,
) -> TreeMonophylyReport:
    """Assess whether requested taxa form a monophyletic group."""
    tree = load_tree(tree_path)
    requested_taxa = sorted(taxa)
    unique_requested_taxa = sorted(set(requested_taxa))
    duplicate_requested_taxa = sorted(
        taxon for taxon, count in Counter(requested_taxa).items() if count > 1
    )
    if not unique_requested_taxa:
        raise ValueError("monophyly assessment requires at least one requested taxon")

    tree_taxa = set(tree.tip_names)
    missing_requested_taxa = sorted(set(unique_requested_taxa) - tree_taxa)
    present_requested_taxa = sorted(set(unique_requested_taxa) & tree_taxa)
    if not present_requested_taxa:
        if reroot:
            raise ValueError("specified outgroup not in labels of the tree")
        return _monophyly_report_from_node(
            tree,
            tree_path=tree_path,
            requested_taxa=requested_taxa,
            unique_requested_taxa=unique_requested_taxa,
            duplicate_requested_taxa=duplicate_requested_taxa,
            missing_requested_taxa=missing_requested_taxa,
            present_requested_taxa=present_requested_taxa,
            reroot=reroot,
            monophyletic=False,
            complementary_clade_used=False,
            matched_node=None,
        )

    if len(present_requested_taxa) == 1:
        matched_node = next(
            node
            for node in tree.iter_leaves()
            if node.name == present_requested_taxa[0]
        )
        return _monophyly_report_from_node(
            tree,
            tree_path=tree_path,
            requested_taxa=requested_taxa,
            unique_requested_taxa=unique_requested_taxa,
            duplicate_requested_taxa=duplicate_requested_taxa,
            missing_requested_taxa=missing_requested_taxa,
            present_requested_taxa=present_requested_taxa,
            reroot=reroot,
            monophyletic=True,
            complementary_clade_used=False,
            matched_node=matched_node,
        )

    matched_node = _mrca_node_from_taxa(tree, present_requested_taxa)
    matched_taxa = _descendant_taxa(matched_node)
    if matched_taxa == present_requested_taxa:
        return _monophyly_report_from_node(
            tree,
            tree_path=tree_path,
            requested_taxa=requested_taxa,
            unique_requested_taxa=unique_requested_taxa,
            duplicate_requested_taxa=duplicate_requested_taxa,
            missing_requested_taxa=missing_requested_taxa,
            present_requested_taxa=present_requested_taxa,
            reroot=reroot,
            monophyletic=True,
            complementary_clade_used=False,
            matched_node=matched_node,
        )

    complementary_taxa = sorted(tree_taxa - set(present_requested_taxa))
    complementary_clade_used = reroot and len(complementary_taxa) == 1
    return _monophyly_report_from_node(
        tree,
        tree_path=tree_path,
        requested_taxa=requested_taxa,
        unique_requested_taxa=unique_requested_taxa,
        duplicate_requested_taxa=duplicate_requested_taxa,
        missing_requested_taxa=missing_requested_taxa,
        present_requested_taxa=present_requested_taxa,
        reroot=reroot,
        monophyletic=complementary_clade_used,
        complementary_clade_used=complementary_clade_used,
        matched_node=matched_node,
    )

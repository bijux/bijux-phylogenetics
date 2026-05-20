from __future__ import annotations

from collections.abc import Iterator

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def iter_internal_nodes_preorder(node: TreeNode) -> Iterator[TreeNode]:
    """Yield internal nodes in one deterministic preorder traversal."""
    if not node.is_leaf():
        yield node
    for child in node.children:
        yield from iter_internal_nodes_preorder(child)


def build_ape_internal_node_map(tree: PhyloTree) -> dict[int, TreeNode]:
    """Return one ape-style internal-node map with the root first."""
    start_node_id = tree.tip_count + 1
    return {
        start_node_id + offset: node
        for offset, node in enumerate(iter_internal_nodes_preorder(tree.root))
    }


def build_ape_tip_node_map(tree: PhyloTree) -> dict[int, TreeNode]:
    """Return one ape-style tip-node map in deterministic leaf order."""
    return dict(enumerate(tree.iter_leaves(), start=1))


def ape_node_id_for_node(tree: PhyloTree, node: TreeNode) -> int:
    """Resolve one ape-style node id for a node already present in the tree."""
    for node_id, candidate in build_ape_tip_node_map(tree).items():
        if candidate is node:
            return node_id
    for node_id, candidate in build_ape_internal_node_map(tree).items():
        if candidate is node:
            return node_id
    raise ValueError("node is not part of the supplied tree")

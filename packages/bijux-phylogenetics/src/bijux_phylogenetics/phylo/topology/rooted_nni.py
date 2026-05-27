from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id

from .tree import PhyloTree, TreeNode, descendant_taxa


@dataclass(frozen=True, slots=True)
class RootedNniMoveCandidate:
    """One deterministic rooted NNI move over a binary tree."""

    parent_node_id: str
    child_node_id: str
    sibling_node_id: str
    exchanged_child_node_id: str
    pivot_branch_id: str
    sibling_clade_id: str
    exchanged_clade_id: str


def rooted_nni_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort candidate NNI branches deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def rooted_nni_clade_id(node: TreeNode) -> str:
    """Render one branch endpoint as a stable descendant-clade identifier."""
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def require_rooted_nni_node_id(node: TreeNode) -> str:
    """Require one refreshed stable node identifier on a rooted NNI tree."""
    if node.node_id is None:
        raise AssertionError("rooted NNI search requires refreshed node identities")
    return node.node_id


def iter_rooted_nni_move_candidates(tree: PhyloTree):
    """Yield deterministic rooted NNI candidates for one rooted binary tree."""
    for parent in tree.iter_internal_nodes(order="preorder"):
        if len(parent.children) != 2:
            continue
        sorted_parent_children = sorted(parent.children, key=rooted_nni_node_sort_key)
        for child in sorted_parent_children:
            if child.is_leaf() or len(child.children) != 2:
                continue
            sibling = next(
                candidate
                for candidate in sorted_parent_children
                if candidate is not child
            )
            sorted_child_children = sorted(
                child.children,
                key=rooted_nni_node_sort_key,
            )
            for exchanged_child in sorted_child_children:
                yield RootedNniMoveCandidate(
                    parent_node_id=require_rooted_nni_node_id(parent),
                    child_node_id=require_rooted_nni_node_id(child),
                    sibling_node_id=require_rooted_nni_node_id(sibling),
                    exchanged_child_node_id=require_rooted_nni_node_id(
                        exchanged_child
                    ),
                    pivot_branch_id=rooted_nni_clade_id(child),
                    sibling_clade_id=rooted_nni_clade_id(sibling),
                    exchanged_clade_id=rooted_nni_clade_id(exchanged_child),
                )


def apply_rooted_nni_move(
    tree: PhyloTree,
    candidate: RootedNniMoveCandidate,
) -> PhyloTree:
    """Return one copied rooted tree with the selected NNI move applied."""
    swapped_tree = tree.copy().refresh()
    parent = swapped_tree.node_by_id(candidate.parent_node_id)
    child = swapped_tree.node_by_id(candidate.child_node_id)
    sibling = swapped_tree.node_by_id(candidate.sibling_node_id)
    exchanged_child = swapped_tree.node_by_id(candidate.exchanged_child_node_id)
    remaining_child = next(
        branch for branch in child.children if branch is not exchanged_child
    )
    child.replace_children([remaining_child, sibling])
    parent.replace_children([child, exchanged_child])
    return swapped_tree.refresh()

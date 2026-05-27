from __future__ import annotations

from dataclasses import dataclass

from bijux_phylogenetics.phylo.topology.clades import canonical_clade_id

from .tree import PhyloTree, TreeNode, descendant_taxa

_ROOT_REGRAFT_BRANCH_ID = "root"


@dataclass(frozen=True, slots=True)
class RootedSprMoveCandidate:
    """One deterministic rooted SPR move over a binary rooted tree."""

    pruned_node_id: str
    pruned_clade_id: str
    pruned_descendant_taxa: tuple[str, ...]
    regraft_target_branch_id: str
    regraft_target_descendant_taxa: tuple[str, ...] | None


def rooted_spr_node_sort_key(node: TreeNode) -> tuple[int, tuple[str, ...]]:
    """Sort rooted SPR candidates deterministically by descendant taxa."""
    descendants = tuple(descendant_taxa(node))
    return (len(descendants), descendants)


def rooted_spr_clade_id(node: TreeNode) -> str:
    """Render one rooted SPR branch endpoint as a stable descendant-clade identifier."""
    return canonical_clade_id(frozenset(descendant_taxa(node)))


def require_rooted_spr_node_id(node: TreeNode) -> str:
    """Require one refreshed stable node identifier on a rooted SPR tree."""
    if node.node_id is None:
        raise AssertionError("rooted SPR search requires refreshed node identities")
    return node.node_id


def iter_rooted_spr_move_candidates(tree: PhyloTree):
    """Yield deterministic rooted SPR candidates for one rooted tree."""
    seen_candidates: set[tuple[str, str]] = set()
    sorted_prune_nodes = sorted(
        (node for node in tree.iter_nodes(order="preorder") if node is not tree.root),
        key=rooted_spr_node_sort_key,
    )
    for prune_node in sorted_prune_nodes:
        pruned_clade_id = rooted_spr_clade_id(prune_node)
        remainder_tree, _pruned_subtree = prune_rooted_spr_subtree(
            tree,
            require_rooted_spr_node_id(prune_node),
        )
        regraft_targets = [
            (_ROOT_REGRAFT_BRANCH_ID, None),
            *[
                (
                    rooted_spr_clade_id(target_node),
                    tuple(sorted(target_node.descendant_taxa)),
                )
                for target_node in sorted(
                    remainder_tree.iter_nodes(order="preorder"),
                    key=rooted_spr_node_sort_key,
                )
            ],
        ]
        for regraft_target_branch_id, regraft_target_descendant_taxa in regraft_targets:
            signature = (pruned_clade_id, regraft_target_branch_id)
            if signature in seen_candidates:
                continue
            seen_candidates.add(signature)
            yield RootedSprMoveCandidate(
                pruned_node_id=require_rooted_spr_node_id(prune_node),
                pruned_clade_id=pruned_clade_id,
                pruned_descendant_taxa=tuple(sorted(prune_node.descendant_taxa)),
                regraft_target_branch_id=regraft_target_branch_id,
                regraft_target_descendant_taxa=regraft_target_descendant_taxa,
            )


def apply_rooted_spr_move(
    tree: PhyloTree,
    candidate: RootedSprMoveCandidate,
) -> PhyloTree:
    """Return one copied rooted tree with the selected SPR move applied."""
    tree_has_explicit_branch_lengths = any(
        child.branch_length is not None for _parent, child in tree.iter_edges()
    )
    remainder_tree, pruned_subtree = prune_rooted_spr_subtree(
        tree,
        candidate.pruned_node_id,
    )
    if candidate.regraft_target_branch_id == _ROOT_REGRAFT_BRANCH_ID:
        if tree_has_explicit_branch_lengths:
            if remainder_tree.root.branch_length is None:
                remainder_tree.root.branch_length = _seed_missing_rooted_spr_branch_length(
                    pruned_subtree
                )
            if pruned_subtree.branch_length is None:
                pruned_subtree.branch_length = _seed_missing_rooted_spr_branch_length(
                    remainder_tree.root
                )
        return PhyloTree(
            root=TreeNode(children=[remainder_tree.root, pruned_subtree]),
            source_format=remainder_tree.source_format,
            rooted=True,
        ).refresh()
    regraft_target = _find_node_by_clade_id(
        remainder_tree,
        candidate.regraft_target_branch_id,
    )
    target_parent = regraft_target.parent
    original_target_branch_length = regraft_target.branch_length
    inserted = TreeNode(children=[regraft_target, pruned_subtree])
    if tree_has_explicit_branch_lengths:
        if original_target_branch_length is not None:
            inserted.branch_length = original_target_branch_length / 2.0
            regraft_target.branch_length = original_target_branch_length / 2.0
        else:
            seeded_length = _seed_missing_rooted_spr_branch_length(pruned_subtree)
            inserted.branch_length = seeded_length
            regraft_target.branch_length = seeded_length
    if target_parent is None:
        return PhyloTree(
            root=inserted,
            source_format=remainder_tree.source_format,
            rooted=True,
        ).refresh()
    target_parent.replace_children(
        [
            child if child is not regraft_target else inserted
            for child in target_parent.children
        ]
    )
    return remainder_tree.refresh()


def prune_rooted_spr_subtree(
    tree: PhyloTree,
    prune_node_id: str,
) -> tuple[PhyloTree, TreeNode]:
    """Detach one non-root subtree and return the reduced tree plus detached clade."""
    working_tree = tree.copy().refresh()
    prune_node = working_tree.node_by_id(prune_node_id)
    parent = prune_node.parent
    if parent is None:
        raise AssertionError("rooted SPR prune node must not be the root")
    sibling = next(child for child in parent.children if child is not prune_node)
    grandparent = parent.parent
    detached_subtree = prune_node.copy()
    if grandparent is None:
        sibling.parent = None
        remainder_tree = PhyloTree(
            root=sibling,
            source_format=working_tree.source_format,
            rooted=True,
        )
        return remainder_tree.refresh(), detached_subtree
    grandparent.replace_children(
        [child if child is not parent else sibling for child in grandparent.children]
    )
    return working_tree.refresh(), detached_subtree


def _find_node_by_clade_id(tree: PhyloTree, clade_id: str) -> TreeNode:
    for node in tree.iter_nodes(order="preorder"):
        if rooted_spr_clade_id(node) == clade_id:
            return node
    raise KeyError(f"tree does not contain clade_id '{clade_id}'")


def _seed_missing_rooted_spr_branch_length(node: TreeNode) -> float:
    branch_length = node.branch_length
    if branch_length is not None:
        return branch_length
    return 0.0

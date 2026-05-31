from __future__ import annotations

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode

from .fixed_topology_policy import validate_fixed_topology_distance_input


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


def score_balanced_minimum_evolution(
    tree: PhyloTree,
    identifiers: list[str],
    distance_lookup: dict[tuple[str, str], float],
) -> float:
    """Score one fixed topology by the Pauplin balanced minimum-evolution criterion."""
    validate_fixed_topology_distance_input(tree, identifiers, distance_lookup)
    paths = _root_to_tip_paths(tree)
    total = 0.0
    for left_index, left_identifier in enumerate(identifiers):
        left_path = paths[left_identifier]
        for right_identifier in identifiers[left_index + 1 :]:
            right_path = paths[right_identifier]
            prefix_length = _common_prefix_length(left_path, right_path)
            non_root_internal_count = sum(
                1
                for node in left_path[prefix_length - 1 : -1]
                if node is not tree.root and not node.is_leaf()
            ) + sum(
                1
                for node in right_path[prefix_length:-1]
                if node is not tree.root and not node.is_leaf()
            )
            total += distance_lookup[(left_identifier, right_identifier)] * (
                0.5**non_root_internal_count
            )
    return round(total, 12)

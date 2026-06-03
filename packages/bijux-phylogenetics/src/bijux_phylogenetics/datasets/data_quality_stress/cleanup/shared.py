from __future__ import annotations

from bijux_phylogenetics.phylo.topology.tree import TreeNode


def apply_branch_length_floor(root: TreeNode, *, floor: float) -> list[str]:
    repaired_nodes: list[str] = []

    def visit(node: TreeNode, *, is_root: bool) -> None:
        if not is_root and node.branch_length is not None and node.branch_length <= 0.0:
            node.branch_length = floor
            repaired_nodes.append(node.name or "<internal>")
        for child in node.children:
            visit(child, is_root=False)

    visit(root, is_root=True)
    return repaired_nodes

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.node_identity import (
    build_ape_internal_node_map,
    build_ape_tip_node_map,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError


@dataclass(slots=True)
class TreeNodeDepthRow:
    """One ape-style node-depth row for a tree node."""

    node_id: int
    node_kind: str
    node_label: str | None
    descendant_taxa: list[str]
    branch_length_depth: float
    branch_length: float | None


@dataclass(slots=True)
class TreeNodeDepthReport:
    """Deterministic ape-style branch-length depths for all nodes in one tree."""

    tree_path: Path
    tip_labels: list[str]
    rooted: bool | None
    tree_is_ultrametric: bool
    node_count: int
    tip_count: int
    internal_node_count: int
    branch_length_count: int
    expected_branch_length_count: int
    complete_branch_lengths: bool
    zero_branch_length_count: int
    minimum_tip_depth: float
    maximum_tip_depth: float
    minimum_internal_depth: float
    maximum_internal_depth: float
    rows: list[TreeNodeDepthRow]


def compute_tree_node_depths(path: Path) -> TreeNodeDepthReport:
    """Compute ape-style branch-length node depths for one tree."""
    tree = load_tree(path)
    return _summarize_tree_node_depths(tree, tree_path=path)


def write_tree_node_depth_table(path: Path, report: TreeNodeDepthReport) -> Path:
    """Write one deterministic long-form ape-style node-depth ledger."""
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        (
            "node_id\tnode_kind\tnode_label\tdescendant_taxa\tbranch_length_depth\t"
            "branch_length\trooted\ttree_is_ultrametric\ttip_count\tinternal_node_count\t"
            "node_count\tzero_branch_length_count\tminimum_tip_depth\tmaximum_tip_depth\t"
            "minimum_internal_depth\tmaximum_internal_depth"
        )
    ]
    for row in report.rows:
        lines.append(
            "\t".join(
                [
                    str(row.node_id),
                    row.node_kind,
                    row.node_label or "",
                    "|".join(row.descendant_taxa),
                    format(row.branch_length_depth, ".15g"),
                    ""
                    if row.branch_length is None
                    else format(row.branch_length, ".15g"),
                    str(report.rooted),
                    str(report.tree_is_ultrametric),
                    str(report.tip_count),
                    str(report.internal_node_count),
                    str(report.node_count),
                    str(report.zero_branch_length_count),
                    format(report.minimum_tip_depth, ".15g"),
                    format(report.maximum_tip_depth, ".15g"),
                    format(report.minimum_internal_depth, ".15g"),
                    format(report.maximum_internal_depth, ".15g"),
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _summarize_tree_node_depths(
    tree: PhyloTree, *, tree_path: Path
) -> TreeNodeDepthReport:
    branch_length_count, expected_branch_length_count = _branch_length_counts(tree)
    if branch_length_count != expected_branch_length_count:
        raise InvalidBranchLengthError(
            "tree requires complete branch lengths for node-depth calculations",
            code="tree_node_depth_missing_branch_lengths",
            details={
                "branch_length_count": branch_length_count,
                "expected_branch_length_count": expected_branch_length_count,
            },
        )
    depth_lookup = _node_depth_lookup(tree)
    rows = _build_node_depth_rows(tree, depth_lookup)
    tip_depths = [row.branch_length_depth for row in rows if row.node_kind == "tip"]
    internal_depths = [
        row.branch_length_depth for row in rows if row.node_kind != "tip"
    ]
    zero_branch_length_count = sum(
        1
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length == 0.0
    )
    minimum_tip_depth = min(tip_depths)
    maximum_tip_depth = max(tip_depths)
    minimum_internal_depth = min(internal_depths)
    maximum_internal_depth = max(internal_depths)
    return TreeNodeDepthReport(
        tree_path=tree_path,
        tip_labels=list(tree.tip_names),
        rooted=_interpreted_rooted_state(tree),
        tree_is_ultrametric=abs(maximum_tip_depth - minimum_tip_depth) <= 1e-12,
        node_count=len(rows),
        tip_count=len(tip_depths),
        internal_node_count=len(internal_depths),
        branch_length_count=branch_length_count,
        expected_branch_length_count=expected_branch_length_count,
        complete_branch_lengths=True,
        zero_branch_length_count=zero_branch_length_count,
        minimum_tip_depth=minimum_tip_depth,
        maximum_tip_depth=maximum_tip_depth,
        minimum_internal_depth=minimum_internal_depth,
        maximum_internal_depth=maximum_internal_depth,
        rows=rows,
    )


def _build_node_depth_rows(
    tree: PhyloTree, depth_lookup: dict[str, float]
) -> list[TreeNodeDepthRow]:
    rows: list[TreeNodeDepthRow] = []
    root_node_id = tree.tip_count + 1
    for node_id, node in build_ape_tip_node_map(tree).items():
        rows.append(
            TreeNodeDepthRow(
                node_id=node_id,
                node_kind="tip",
                node_label=node.name,
                descendant_taxa=[node.name] if node.name is not None else [],
                branch_length_depth=depth_lookup[node.node_id or ""],
                branch_length=node.branch_length,
            )
        )
    for node_id, node in build_ape_internal_node_map(tree).items():
        rows.append(
            TreeNodeDepthRow(
                node_id=node_id,
                node_kind="root" if node_id == root_node_id else "internal",
                node_label=node.name,
                descendant_taxa=_descendant_taxa(node),
                branch_length_depth=depth_lookup[node.node_id or ""],
                branch_length=node.branch_length,
            )
        )
    return rows


def _node_depth_lookup(tree: PhyloTree) -> dict[str, float]:
    root_id = tree.root.node_id or ""
    depths: dict[str, float] = {root_id: 0.0}

    def visit(node: TreeNode) -> None:
        node_id = node.node_id or ""
        base_depth = depths[node_id]
        for child in node.children:
            depths[child.node_id or ""] = base_depth + float(child.branch_length or 0.0)
            visit(child)

    visit(tree.root)
    return depths


def _branch_length_counts(tree: PhyloTree) -> tuple[int, int]:
    branch_length_count = 0
    expected_branch_length_count = 0
    for node in tree.iter_nodes():
        if node is tree.root:
            continue
        expected_branch_length_count += 1
        if node.branch_length is not None:
            branch_length_count += 1
    return branch_length_count, expected_branch_length_count


def _descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _interpreted_rooted_state(tree: PhyloTree) -> bool:
    if tree.rooted is True:
        return True
    return len(tree.root.children) == 2

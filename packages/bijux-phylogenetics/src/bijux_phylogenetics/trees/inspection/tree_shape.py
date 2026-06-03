from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import math
from pathlib import Path
from statistics import median

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BioTree

from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.trees import detect_tree_format, load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


@dataclass(frozen=True, slots=True)
class TreeShapeRow:
    source_path: Path
    tree_index: int | None
    rooted: bool | None
    tip_count: int
    internal_node_count: int
    is_binary: bool
    cherry_count: int
    sackin_imbalance_index: int
    colless_imbalance_index: float | None
    normalized_colless_imbalance: float | None
    tree_height_edges: int
    tree_height_branch_length: float | None
    mean_tip_depth_edges: float
    mean_root_to_tip_branch_length: float | None
    imbalance_summary: str
    ladderized: bool
    star_like: bool
    comb_like: bool
    unusually_imbalanced: bool | None


@dataclass(frozen=True, slots=True)
class TreeShapeAggregate:
    balanced_tree_count: int
    skewed_tree_count: int
    ladderized_tree_count: int
    star_like_tree_count: int
    comb_like_tree_count: int
    binary_tree_count: int
    mean_cherry_count: float
    mean_sackin_imbalance_index: float
    median_sackin_imbalance_index: float
    maximum_sackin_imbalance_index: int
    mean_tree_height_edges: float
    maximum_tree_height_edges: int
    colless_defined_tree_count: int
    mean_colless_imbalance_index: float | None
    mean_normalized_colless_imbalance: float | None
    branch_length_height_defined_tree_count: int
    mean_tree_height_branch_length: float | None


@dataclass(slots=True)
class TreeShapeReport:
    path: Path
    source_format: str
    tree_count: int
    rows: list[TreeShapeRow]
    aggregate: TreeShapeAggregate


def _round_float(value: float) -> float:
    return round(value, 15)


def _mean(values: list[float]) -> float:
    return _round_float(sum(values) / len(values))


def _load_tree_set(path: Path) -> tuple[str, list[BioTree], list[PhyloTree]]:
    if not path.exists():
        raise FileNotFoundError(f"tree-set file not found: {path}")
    source_format = detect_tree_format(path)
    bio_trees = list(Phylo.parse(path, source_format))
    if not bio_trees:
        raise InvalidAlignmentError(f"tree set contains no trees: {path}")
    trees = [
        tree_from_biophylo(tree, source_format=source_format) for tree in bio_trees
    ]
    return source_format, bio_trees, trees


def _leaf_depths(tree: PhyloTree) -> list[int]:
    depths: list[int] = []

    def visit(node: TreeNode, depth: int) -> None:
        if node.is_leaf():
            depths.append(depth)
            return
        for child in node.children:
            visit(child, depth + 1)

    visit(tree.root, 0)
    return depths


def _is_binary(tree: PhyloTree) -> bool:
    return all(node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes())


def _cherry_count(tree: PhyloTree) -> int:
    return sum(
        1
        for node in tree.iter_nodes()
        if len(node.children) == 2 and all(child.is_leaf() for child in node.children)
    )


def _tree_height_edges(tree: PhyloTree) -> int:
    return max(_leaf_depths(tree), default=0)


def _tree_height_branch_length(tree: PhyloTree) -> float | None:
    lengths = tree.root_to_tip_lengths()
    if not lengths or any(length is None for length in lengths):
        return None
    defined_lengths = [float(length) for length in lengths if length is not None]
    return _round_float(max(defined_lengths))


def _mean_root_to_tip_branch_length(tree: PhyloTree) -> float | None:
    lengths = tree.root_to_tip_lengths()
    if not lengths or any(length is None for length in lengths):
        return None
    return _mean([float(length) for length in lengths if length is not None])


def _colless_imbalance_index(tree: PhyloTree) -> float | None:
    if len(tree.root.children) != 2:
        return None

    def visit(node: TreeNode) -> tuple[int, float | None]:
        if node.is_leaf():
            return 1, 0.0
        if len(node.children) != 2:
            return sum(visit(child)[0] for child in node.children), None
        left_count, left_score = visit(node.children[0])
        right_count, right_score = visit(node.children[1])
        if left_score is None or right_score is None:
            return left_count + right_count, None
        return left_count + right_count, left_score + right_score + abs(
            left_count - right_count
        )

    _, score = visit(tree.root)
    return None if score is None else float(score)


def _normalized_colless_imbalance(tree: PhyloTree) -> float | None:
    raw_score = _colless_imbalance_index(tree)
    if raw_score is None:
        return None
    if tree.tip_count <= 2:
        return 0.0
    maximum = ((tree.tip_count - 1) * (tree.tip_count - 2)) / 2
    return _round_float(raw_score / maximum)


def _unusually_imbalanced(tree: PhyloTree, *, threshold: float = 0.75) -> bool | None:
    normalized = _normalized_colless_imbalance(tree)
    if normalized is None:
        return None
    return normalized >= threshold


def _imbalance_summary(tree: PhyloTree) -> str:
    def visit(node: TreeNode) -> tuple[int, int]:
        if node.is_leaf():
            return 1, 0
        child_counts: list[int] = []
        total_score = 0
        for child in node.children:
            tip_count, child_score = visit(child)
            child_counts.append(tip_count)
            total_score += child_score
        if len(child_counts) >= 2:
            total_score += max(child_counts) - min(child_counts)
        return sum(child_counts), total_score

    _, score = visit(tree.root)
    normalized_score = score / max(tree.tip_count - 1, 1)
    if normalized_score == 0:
        return "balanced"
    if normalized_score >= 1:
        return "ladderized"
    return "skewed"


def _star_like(tree: PhyloTree) -> bool:
    threshold = max(4, math.ceil(tree.tip_count * 0.75))
    for node in tree.iter_nodes():
        leaf_children = sum(1 for child in node.children if child.is_leaf())
        if leaf_children == len(node.children) and leaf_children >= threshold:
            return True
    return False


def _comb_like(tree: PhyloTree) -> bool:
    if len(tree.root.children) != 2:
        return False
    if any(
        not node.is_leaf() and len(node.children) != 2 for node in tree.iter_nodes()
    ):
        return False
    return _cherry_count(tree) == 1 and _tree_height_edges(tree) == tree.tip_count - 1


def summarize_tree_shape_from_tree(
    tree: PhyloTree,
    *,
    source_path: Path,
    tree_index: int | None = None,
) -> TreeShapeRow:
    depths = _leaf_depths(tree)
    colless_imbalance_index = _colless_imbalance_index(tree)
    normalized_colless_imbalance = _normalized_colless_imbalance(tree)
    imbalance_summary = _imbalance_summary(tree)
    star_like = _star_like(tree)
    comb_like = _comb_like(tree)
    return TreeShapeRow(
        source_path=source_path,
        tree_index=tree_index,
        rooted=tree.rooted,
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        is_binary=_is_binary(tree),
        cherry_count=_cherry_count(tree),
        sackin_imbalance_index=sum(depths),
        colless_imbalance_index=colless_imbalance_index,
        normalized_colless_imbalance=normalized_colless_imbalance,
        tree_height_edges=max(depths, default=0),
        tree_height_branch_length=_tree_height_branch_length(tree),
        mean_tip_depth_edges=_mean([float(depth) for depth in depths]),
        mean_root_to_tip_branch_length=_mean_root_to_tip_branch_length(tree),
        imbalance_summary=imbalance_summary,
        ladderized=imbalance_summary == "ladderized",
        star_like=star_like,
        comb_like=comb_like,
        unusually_imbalanced=_unusually_imbalanced(tree),
    )


def _aggregate_tree_shapes(rows: list[TreeShapeRow]) -> TreeShapeAggregate:
    sackin_values = [row.sackin_imbalance_index for row in rows]
    cherry_values = [float(row.cherry_count) for row in rows]
    height_edges = [float(row.tree_height_edges) for row in rows]
    colless_values = [
        row.colless_imbalance_index
        for row in rows
        if row.colless_imbalance_index is not None
    ]
    normalized_colless_values = [
        row.normalized_colless_imbalance
        for row in rows
        if row.normalized_colless_imbalance is not None
    ]
    branch_length_heights = [
        row.tree_height_branch_length
        for row in rows
        if row.tree_height_branch_length is not None
    ]
    return TreeShapeAggregate(
        balanced_tree_count=sum(
            1 for row in rows if row.imbalance_summary == "balanced"
        ),
        skewed_tree_count=sum(1 for row in rows if row.imbalance_summary == "skewed"),
        ladderized_tree_count=sum(1 for row in rows if row.ladderized),
        star_like_tree_count=sum(1 for row in rows if row.star_like),
        comb_like_tree_count=sum(1 for row in rows if row.comb_like),
        binary_tree_count=sum(1 for row in rows if row.is_binary),
        mean_cherry_count=_mean(cherry_values),
        mean_sackin_imbalance_index=_mean([float(value) for value in sackin_values]),
        median_sackin_imbalance_index=float(median(sackin_values)),
        maximum_sackin_imbalance_index=max(sackin_values),
        mean_tree_height_edges=_mean(height_edges),
        maximum_tree_height_edges=max(row.tree_height_edges for row in rows),
        colless_defined_tree_count=len(colless_values),
        mean_colless_imbalance_index=_mean(colless_values) if colless_values else None,
        mean_normalized_colless_imbalance=(
            _mean(normalized_colless_values) if normalized_colless_values else None
        ),
        branch_length_height_defined_tree_count=len(branch_length_heights),
        mean_tree_height_branch_length=(
            _mean(branch_length_heights) if branch_length_heights else None
        ),
    )


def summarize_tree_shape(
    path: Path, *, source_format: str | None = None
) -> TreeShapeReport:
    tree = load_tree(path, source_format=source_format)
    row = summarize_tree_shape_from_tree(tree, source_path=path)
    return TreeShapeReport(
        path=path,
        source_format=tree.source_format,
        tree_count=1,
        rows=[row],
        aggregate=_aggregate_tree_shapes([row]),
    )


def summarize_tree_set_shapes(path: Path) -> TreeShapeReport:
    source_format, _, trees = _load_tree_set(path)
    rows = [
        summarize_tree_shape_from_tree(tree, source_path=path, tree_index=index)
        for index, tree in enumerate(trees, start=1)
    ]
    return TreeShapeReport(
        path=path,
        source_format=source_format,
        tree_count=len(rows),
        rows=rows,
        aggregate=_aggregate_tree_shapes(rows),
    )


def write_tree_shape_table(path: Path, report: TreeShapeReport) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = list(asdict(report.rows[0]).keys())
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in report.rows:
            payload = asdict(row)
            payload["source_path"] = str(row.source_path)
            writer.writerow(payload)
    return path

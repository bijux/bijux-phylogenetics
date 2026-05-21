from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from pathlib import Path
from statistics import median

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree as BioTree

from bijux_phylogenetics.io.biopython import tree_from_biophylo
from bijux_phylogenetics.io.trees import detect_tree_format, load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


@dataclass(frozen=True, slots=True)
class BranchLengthRow:
    source_path: Path
    tree_index: int | None
    node: str
    branch_type: str
    tip_taxon: str | None
    descendant_tip_count: int
    branch_length: float | None
    root_depth: float | None
    tree_positive_branch_median: float | None
    zero_length: bool
    negative_length: bool
    missing_length: bool
    long_outlier: bool
    short_outlier: bool
    outlier_class: str


@dataclass(frozen=True, slots=True)
class BranchLengthAggregate:
    branch_count: int
    defined_branch_count: int
    missing_branch_count: int
    zero_length_branch_count: int
    negative_branch_count: int
    positive_branch_count: int
    long_outlier_count: int
    short_outlier_count: int
    minimum_branch_length: float | None
    maximum_branch_length: float | None
    mean_branch_length: float | None
    median_branch_length: float | None
    positive_branch_median: float | None


@dataclass(slots=True)
class BranchLengthDistributionReport:
    path: Path
    source_format: str
    tree_count: int
    rows: list[BranchLengthRow]
    aggregate: BranchLengthAggregate


def _round_float(value: float) -> float:
    return round(value, 15)


def _mean(values: list[float]) -> float:
    return _round_float(sum(values) / len(values))


def _node_signature(node: TreeNode) -> str:
    if node.is_leaf():
        return node.name or "<unnamed-tip>"
    tips = sorted(child.name for child in node.iter_leaves() if child.name)
    return "|".join(tips) if tips else "<internal>"


def _descendant_tip_count(node: TreeNode) -> int:
    return sum(1 for _ in node.iter_leaves())


def _branch_root_depths(tree: PhyloTree) -> dict[str, float | None]:
    depths: dict[str, float | None] = {tree.root.node_id or "": 0.0}

    def visit(node: TreeNode, depth: float | None) -> None:
        depths[node.node_id or ""] = None if depth is None else _round_float(depth)
        for child in node.children:
            if depth is None or child.branch_length is None:
                child_depth = None
            else:
                child_depth = depth + child.branch_length
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return depths


def _positive_branch_median(tree: PhyloTree) -> float | None:
    positive = sorted(
        float(node.branch_length)
        for node in tree.iter_nodes()
        if node is not tree.root
        and node.branch_length is not None
        and node.branch_length > 0
    )
    if not positive:
        return None
    return _round_float(float(median(positive)))


def _classify_outlier(
    branch_length: float | None,
    *,
    positive_branch_median: float | None,
    long_factor: float = 3.0,
    short_factor: float = 0.1,
) -> tuple[bool, bool, str]:
    if branch_length is None:
        return False, False, "missing"
    if branch_length < 0:
        return False, False, "negative"
    if branch_length == 0:
        return False, False, "zero"
    if positive_branch_median is None:
        return False, False, "typical"
    long_threshold = positive_branch_median * long_factor
    short_threshold = positive_branch_median * short_factor
    if branch_length > long_threshold:
        return True, False, "long"
    if 0 < branch_length < short_threshold:
        return False, True, "short"
    return False, False, "typical"


def _rows_from_tree(
    tree: PhyloTree,
    *,
    source_path: Path,
    tree_index: int | None = None,
) -> list[BranchLengthRow]:
    root_depths = _branch_root_depths(tree)
    positive_branch_median = _positive_branch_median(tree)
    rows: list[BranchLengthRow] = []
    for node in tree.iter_nodes():
        if node is tree.root:
            continue
        branch_length = (
            None
            if node.branch_length is None
            else _round_float(float(node.branch_length))
        )
        long_outlier, short_outlier, outlier_class = _classify_outlier(
            branch_length,
            positive_branch_median=positive_branch_median,
        )
        rows.append(
            BranchLengthRow(
                source_path=source_path,
                tree_index=tree_index,
                node=_node_signature(node),
                branch_type="terminal" if node.is_leaf() else "internal",
                tip_taxon=node.name if node.is_leaf() else None,
                descendant_tip_count=_descendant_tip_count(node),
                branch_length=branch_length,
                root_depth=root_depths.get(node.node_id or ""),
                tree_positive_branch_median=positive_branch_median,
                zero_length=branch_length == 0 if branch_length is not None else False,
                negative_length=branch_length < 0
                if branch_length is not None
                else False,
                missing_length=branch_length is None,
                long_outlier=long_outlier,
                short_outlier=short_outlier,
                outlier_class=outlier_class,
            )
        )
    return rows


def _aggregate_rows(rows: list[BranchLengthRow]) -> BranchLengthAggregate:
    defined = [row.branch_length for row in rows if row.branch_length is not None]
    positive = [
        row.branch_length
        for row in rows
        if row.branch_length is not None and row.branch_length > 0
    ]
    return BranchLengthAggregate(
        branch_count=len(rows),
        defined_branch_count=len(defined),
        missing_branch_count=sum(1 for row in rows if row.missing_length),
        zero_length_branch_count=sum(1 for row in rows if row.zero_length),
        negative_branch_count=sum(1 for row in rows if row.negative_length),
        positive_branch_count=len(positive),
        long_outlier_count=sum(1 for row in rows if row.long_outlier),
        short_outlier_count=sum(1 for row in rows if row.short_outlier),
        minimum_branch_length=min(defined) if defined else None,
        maximum_branch_length=max(defined) if defined else None,
        mean_branch_length=_mean(defined) if defined else None,
        median_branch_length=_round_float(float(median(defined))) if defined else None,
        positive_branch_median=_round_float(float(median(positive)))
        if positive
        else None,
    )


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


def analyze_branch_length_distribution(
    path: Path,
    *,
    source_format: str | None = None,
) -> BranchLengthDistributionReport:
    tree = load_tree(path, source_format=source_format)
    rows = _rows_from_tree(tree, source_path=path)
    return BranchLengthDistributionReport(
        path=path,
        source_format=tree.source_format,
        tree_count=1,
        rows=rows,
        aggregate=_aggregate_rows(rows),
    )


def analyze_tree_set_branch_lengths(path: Path) -> BranchLengthDistributionReport:
    source_format, _, trees = _load_tree_set(path)
    rows: list[BranchLengthRow] = []
    for index, tree in enumerate(trees, start=1):
        rows.extend(_rows_from_tree(tree, source_path=path, tree_index=index))
    return BranchLengthDistributionReport(
        path=path,
        source_format=source_format,
        tree_count=len(trees),
        rows=rows,
        aggregate=_aggregate_rows(rows),
    )


def write_branch_length_table(
    path: Path, report: BranchLengthDistributionReport
) -> Path:
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

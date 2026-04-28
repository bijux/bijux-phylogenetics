from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from pathlib import Path

from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.errors import (
    DuplicateTaxonError,
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    UnnamedTipError,
    UnrootedTreeError,
)
from bijux_phylogenetics.io.trees import load_tree


@dataclass(slots=True)
class TreeValidationReport:
    path: Path
    source_format: str
    tip_count: int
    internal_node_count: int
    rooted: bool
    has_complete_branch_lengths: bool
    branch_length_status: str
    total_branch_length: float
    zero_length_branches: int
    negative_branch_lengths: int
    polytomy_count: int
    polytomy_nodes: list[str]
    missing_taxa: int
    duplicate_taxa: list[str]
    ultrametric: bool | None
    warnings: list[str]


@dataclass(slots=True)
class TreeInspectionReport:
    path: Path
    source_format: str
    tip_count: int
    node_count: int
    internal_node_count: int
    edge_count: int
    clade_count: int
    rooted: bool
    is_binary: bool
    polytomy_count: int
    polytomy_nodes: list[str]
    has_branch_lengths: bool
    branch_length_status: str
    total_branch_length: float
    zero_length_branch_count: int
    min_root_to_tip: float | None
    max_root_to_tip: float | None
    max_depth: int
    mean_depth: float
    imbalance_summary: str
    cherry_count: int
    taxa: list[str]


@dataclass(slots=True)
class TreeDiagnosticReport:
    path: Path
    inspection: TreeInspectionReport
    validation: TreeValidationReport


def _load_tree(path: Path, *, source_format: str | None = None) -> PhyloTree:
    return load_tree(path, source_format=source_format)


def _count_polytomies(tree: PhyloTree) -> int:
    return sum(1 for node in tree.iter_nodes() if not node.is_leaf() and len(node.children) > 2)


def _polytomy_nodes(tree: PhyloTree) -> list[str]:
    nodes: list[str] = []

    def visit(node) -> list[str]:
        if node.is_leaf():
            return [node.name] if node.name is not None else []
        taxa: list[str] = []
        for child in node.children:
            taxa.extend(visit(child))
        if len(node.children) > 2:
            nodes.append("|".join(sorted(taxa)))
        return taxa

    visit(tree.root)
    return nodes


def _edge_count(tree: PhyloTree) -> int:
    return sum(1 for node in tree.iter_nodes() if node is not tree.root)


def _node_count(tree: PhyloTree) -> int:
    return sum(1 for _ in tree.iter_nodes())


def _branch_length_health(tree: PhyloTree) -> tuple[bool, int, int]:
    branch_lengths = [node.branch_length for node in tree.iter_nodes() if node is not tree.root]
    has_complete = all(length is not None for length in branch_lengths) if branch_lengths else False
    zero_count = sum(1 for length in branch_lengths if length == 0)
    negative_count = sum(1 for length in branch_lengths if length is not None and length < 0)
    return has_complete, zero_count, negative_count


def _branch_length_status(tree: PhyloTree) -> str:
    branch_lengths = [node.branch_length for node in tree.iter_nodes() if node is not tree.root]
    if not branch_lengths or all(length is None for length in branch_lengths):
        return "absent"
    if all(length is not None for length in branch_lengths):
        return "complete"
    return "partial"


def _duplicate_taxa(tree: PhyloTree) -> tuple[int, list[str]]:
    names = [name for name in tree.tip_names if name]
    counts = Counter(names)
    duplicates = sorted(name for name, count in counts.items() if count > 1)
    missing = sum(1 for node in tree.iter_leaves() if not node.name)
    return missing, duplicates


def _ultrametric(tree: PhyloTree, *, tolerance: float = 1e-9) -> bool | None:
    lengths = tree.root_to_tip_lengths()
    if not lengths or any(length is None for length in lengths):
        return None
    numeric_lengths = [float(length) for length in lengths if length is not None]
    return max(numeric_lengths) - min(numeric_lengths) <= tolerance


def _max_depth(tree: PhyloTree) -> int:
    def visit(node, depth: int) -> int:
        if node.is_leaf():
            return depth
        return max(visit(child, depth + 1) for child in node.children)

    return visit(tree.root, 0)


def _leaf_depths(tree: PhyloTree) -> list[int]:
    depths: list[int] = []

    def visit(node, depth: int) -> None:
        if node.is_leaf():
            depths.append(depth)
            return
        for child in node.children:
            visit(child, depth + 1)

    visit(tree.root, 0)
    return depths


def _imbalance_summary(tree: PhyloTree) -> str:
    def visit(node) -> tuple[int, int]:
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


def _cherry_count(tree: PhyloTree) -> int:
    return sum(
        1
        for node in tree.iter_nodes()
        if len(node.children) == 2 and all(child.is_leaf() for child in node.children)
    )


def validate_tree_path(
    path: Path,
    *,
    source_format: str | None = None,
    allow_duplicates: bool = False,
    strict: bool = False,
    allow_negative_branch_lengths: bool = False,
    require_rooted: bool = False,
    require_ultrametric: bool = False,
) -> TreeValidationReport:
    """Validate a tree file and produce a diagnostic report."""
    tree = _load_tree(path, source_format=source_format)
    rooted = len(tree.root.children) == 2
    has_complete, zero_count, negative_count = _branch_length_health(tree)
    branch_length_status = _branch_length_status(tree)
    missing_taxa, duplicate_taxa = _duplicate_taxa(tree)
    if duplicate_taxa and not allow_duplicates:
        raise DuplicateTaxonError(f"duplicate tip labels found: {', '.join(duplicate_taxa)}")
    if missing_taxa and strict:
        raise UnnamedTipError(f"tree contains {missing_taxa} unnamed tip labels")
    if negative_count and not allow_negative_branch_lengths:
        raise InvalidBranchLengthError(f"tree contains {negative_count} negative branch lengths")
    ultrametric = _ultrametric(tree)
    if require_rooted and not rooted:
        raise UnrootedTreeError(f"tree is not rooted: {path}")
    if require_ultrametric and ultrametric is not True:
        raise NonUltrametricTreeError(f"tree is not ultrametric within default validation tolerance: {path}")
    polytomy_nodes = _polytomy_nodes(tree)
    warnings: list[str] = []
    if duplicate_taxa:
        warnings.append("tree contains duplicate tip labels")
    if missing_taxa:
        warnings.append("tree contains unnamed tips")
    if negative_count:
        warnings.append("tree contains negative branch lengths")
    if zero_count:
        warnings.append("tree contains zero-length branches")
    if polytomy_nodes:
        warnings.append("tree contains one or more polytomies")
    return TreeValidationReport(
        path=path,
        source_format=tree.source_format,
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        rooted=rooted,
        has_complete_branch_lengths=has_complete,
        branch_length_status=branch_length_status,
        total_branch_length=tree.total_branch_length(),
        zero_length_branches=zero_count,
        negative_branch_lengths=negative_count,
        polytomy_count=len(polytomy_nodes),
        polytomy_nodes=polytomy_nodes,
        missing_taxa=missing_taxa,
        duplicate_taxa=duplicate_taxa,
        ultrametric=ultrametric,
        warnings=warnings,
    )


def inspect_tree_path(path: Path, *, source_format: str | None = None) -> TreeInspectionReport:
    """Inspect a tree file and return lightweight summary metrics."""
    tree = _load_tree(path, source_format=source_format)
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    branch_lengths = [node.branch_length for node in tree.iter_nodes() if node is not tree.root]
    polytomy_nodes = _polytomy_nodes(tree)
    branch_length_status = _branch_length_status(tree)
    depths = _leaf_depths(tree)
    return TreeInspectionReport(
        path=path,
        source_format=tree.source_format,
        tip_count=tree.tip_count,
        node_count=_node_count(tree),
        internal_node_count=tree.internal_node_count,
        edge_count=_edge_count(tree),
        clade_count=tree.internal_node_count,
        rooted=len(tree.root.children) == 2,
        is_binary=all(node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes()),
        polytomy_count=len(polytomy_nodes),
        polytomy_nodes=polytomy_nodes,
        has_branch_lengths=all(length is not None for length in branch_lengths) if branch_lengths else False,
        branch_length_status=branch_length_status,
        total_branch_length=tree.total_branch_length(),
        zero_length_branch_count=sum(1 for length in branch_lengths if length == 0),
        min_root_to_tip=min(lengths) if lengths else None,
        max_root_to_tip=max(lengths) if lengths else None,
        max_depth=_max_depth(tree),
        mean_depth=sum(depths) / len(depths),
        imbalance_summary=_imbalance_summary(tree),
        cherry_count=_cherry_count(tree),
        taxa=sorted(tree.tip_names),
    )


def diagnose_tree_path(path: Path, *, source_format: str | None = None) -> TreeDiagnosticReport:
    """Return a combined inspection and validation report for one tree."""
    return TreeDiagnosticReport(
        path=path,
        inspection=inspect_tree_path(path, source_format=source_format),
        validation=validate_tree_path(path, source_format=source_format),
    )

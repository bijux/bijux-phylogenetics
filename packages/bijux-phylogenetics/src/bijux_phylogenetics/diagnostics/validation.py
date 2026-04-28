from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from pathlib import Path
from statistics import mean, median

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
    missing_internal_branch_nodes: list[str]
    missing_terminal_branch_taxa: list[str]
    singleton_internal_nodes: list[str]
    missing_taxa: int
    duplicate_taxa: list[str]
    ultrametric: bool | None
    warnings: list[str]


@dataclass(slots=True)
class BranchLengthSummary:
    count: int
    minimum: float
    maximum: float
    mean: float
    median: float
    first_quartile: float
    third_quartile: float


@dataclass(slots=True)
class TreeQualityWarning:
    code: str
    message: str
    penalty: float
    affected_taxa: list[str]
    affected_nodes: list[str]


@dataclass(slots=True)
class InternalNodeChildCount:
    node: str
    child_count: int


@dataclass(slots=True)
class BranchLengthOutlier:
    node: str
    branch_length: float
    branch_type: str


@dataclass(slots=True)
class InternalLabelInterpretation:
    node: str
    label: str
    interpretation: str
    numeric_value: float | None


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
    internal_child_counts: list[InternalNodeChildCount]
    singleton_internal_nodes: list[str]
    polytomy_count: int
    polytomy_nodes: list[str]
    has_branch_lengths: bool
    branch_length_status: str
    missing_internal_branch_nodes: list[str]
    missing_terminal_branch_taxa: list[str]
    is_ultrametric: bool | None
    total_branch_length: float
    branch_length_summary: BranchLengthSummary | None
    tree_diameter: float | None
    zero_length_branch_count: int
    min_root_to_tip: float | None
    max_root_to_tip: float | None
    max_depth: int
    mean_depth: float
    colless_imbalance_index: float | None
    normalized_colless_imbalance: float | None
    sackin_imbalance_index: int
    unusually_imbalanced: bool | None
    long_branch_taxa: list[str]
    long_branch_outliers: list[BranchLengthOutlier]
    short_branch_outliers: list[BranchLengthOutlier]
    suspicious_support_value_ranges: list[str]
    mixed_support_scales: bool
    likely_support_labels: list[InternalLabelInterpretation]
    likely_named_internal_labels: list[InternalLabelInterpretation]
    star_like: bool
    comb_like: bool
    tree_quality_score: float
    tree_quality_warnings: list[TreeQualityWarning]
    imbalance_summary: str
    cherry_count: int
    taxa: list[str]
    warnings: list[str]


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


def _descendant_taxa(node) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(_descendant_taxa(child))
    return sorted(taxa)


def _node_signature(node) -> str:
    taxa = _descendant_taxa(node)
    if taxa:
        return "|".join(taxa)
    return node.name or "<unnamed>"


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


def _internal_node_child_counts(tree: PhyloTree) -> list[InternalNodeChildCount]:
    return [
        InternalNodeChildCount(node=_node_signature(node), child_count=len(node.children))
        for node in tree.iter_nodes()
        if not node.is_leaf()
    ]


def _singleton_internal_nodes(tree: PhyloTree) -> list[str]:
    return sorted(
        _node_signature(node)
        for node in tree.iter_nodes()
        if not node.is_leaf() and len(node.children) == 1
    )


def _missing_internal_branch_nodes(tree: PhyloTree) -> list[str]:
    return sorted(
        _node_signature(node)
        for node in tree.iter_nodes()
        if node is not tree.root and not node.is_leaf() and node.branch_length is None
    )


def _missing_terminal_branch_taxa(tree: PhyloTree) -> list[str]:
    return sorted(
        node.name
        for node in tree.iter_leaves()
        if node.name is not None and node.branch_length is None
    )


def _zero_length_branch_nodes(tree: PhyloTree) -> list[str]:
    return sorted(
        _node_signature(node)
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length == 0
    )


def _negative_branch_nodes(tree: PhyloTree) -> list[str]:
    return sorted(
        _node_signature(node)
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None and node.branch_length < 0
    )


def _branch_length_status(tree: PhyloTree) -> str:
    branch_lengths = [node.branch_length for node in tree.iter_nodes() if node is not tree.root]
    if not branch_lengths or all(length is None for length in branch_lengths):
        return "absent"
    if all(length is not None for length in branch_lengths):
        return "complete"
    return "partial"


def _median(values: list[float]) -> float:
    return float(median(values))


def _quartiles(values: list[float]) -> tuple[float, float]:
    ordered = sorted(values)
    midpoint = len(ordered) // 2
    if len(ordered) % 2 == 0:
        lower = ordered[:midpoint]
        upper = ordered[midpoint:]
    else:
        lower = ordered[:midpoint]
        upper = ordered[midpoint + 1 :]
    if not lower:
        lower = ordered
    if not upper:
        upper = ordered
    return _median(lower), _median(upper)


def _branch_length_summary(tree: PhyloTree) -> BranchLengthSummary | None:
    values = [float(length) for length in tree.branch_lengths() if length is not None]
    if not values:
        return None
    first_quartile, third_quartile = _quartiles(values)

    def stable(value: float) -> float:
        return round(value, 15)

    return BranchLengthSummary(
        count=len(values),
        minimum=stable(min(values)),
        maximum=stable(max(values)),
        mean=stable(float(mean(values))),
        median=stable(_median(values)),
        first_quartile=stable(first_quartile),
        third_quartile=stable(third_quartile),
    )


def _tree_diameter(tree: PhyloTree) -> float | None:
    branch_lengths = tree.branch_lengths()
    if any(length is None for length in branch_lengths):
        return None
    leaves = [node for node in tree.iter_leaves() if node.name is not None]
    if len(leaves) < 2:
        return 0.0

    adjacency: dict[int, list[tuple[int, float]]] = {}

    def connect(parent, child) -> None:
        parent_id = id(parent)
        child_id = id(child)
        adjacency.setdefault(parent_id, []).append((child_id, float(child.branch_length or 0.0)))
        adjacency.setdefault(child_id, []).append((parent_id, float(child.branch_length or 0.0)))
        for grandchild in child.children:
            connect(child, grandchild)

    for child in tree.root.children:
        connect(tree.root, child)

    leaf_ids = [id(node) for node in leaves]

    def distances_from(start_id: int) -> dict[int, float]:
        distances: dict[int, float] = {start_id: 0.0}
        stack: list[int] = [start_id]
        while stack:
            current = stack.pop()
            for neighbor, length in adjacency.get(current, []):
                if neighbor in distances:
                    continue
                distances[neighbor] = distances[current] + length
                stack.append(neighbor)
        return distances

    diameter = 0.0
    for index, leaf_id in enumerate(leaf_ids):
        distances = distances_from(leaf_id)
        for other_id in leaf_ids[index + 1 :]:
            diameter = max(diameter, distances[other_id])
    return round(diameter, 15)


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


def _colless_imbalance_index(tree: PhyloTree) -> float | None:
    if len(tree.root.children) != 2:
        return None

    def visit(node) -> tuple[int, float | None]:
        if node.is_leaf():
            return 1, 0.0
        if len(node.children) != 2:
            return sum(visit(child)[0] for child in node.children), None
        left_count, left_score = visit(node.children[0])
        right_count, right_score = visit(node.children[1])
        if left_score is None or right_score is None:
            return left_count + right_count, None
        return left_count + right_count, left_score + right_score + abs(left_count - right_count)

    _, score = visit(tree.root)
    return None if score is None else float(score)


def _normalized_colless_imbalance(tree: PhyloTree) -> float | None:
    raw_score = _colless_imbalance_index(tree)
    if raw_score is None:
        return None
    if tree.tip_count <= 2:
        return 0.0
    maximum = ((tree.tip_count - 1) * (tree.tip_count - 2)) / 2
    return round(raw_score / maximum, 15)


def _unusually_imbalanced(tree: PhyloTree, *, threshold: float = 0.75) -> bool | None:
    normalized = _normalized_colless_imbalance(tree)
    if normalized is None:
        return None
    return normalized >= threshold


def _comb_like(tree: PhyloTree) -> bool:
    if len(tree.root.children) != 2:
        return False
    if any(not node.is_leaf() and len(node.children) != 2 for node in tree.iter_nodes()):
        return False
    return _cherry_count(tree) == 1 and _max_depth(tree) == tree.tip_count - 1


def _long_branch_taxa(tree: PhyloTree, *, factor: float = 3.0) -> list[str]:
    terminal_lengths = [(name, length) for name, length in tree.terminal_branch_lengths() if length is not None and length > 0]
    if len(terminal_lengths) < 2:
        return []
    baseline = _median(sorted(length for _, length in terminal_lengths))
    threshold = baseline * factor
    return sorted(name for name, length in terminal_lengths if length > threshold)


def _branch_outliers(
    tree: PhyloTree,
    *,
    long_factor: float = 3.0,
    short_factor: float = 0.1,
) -> tuple[list[BranchLengthOutlier], list[BranchLengthOutlier]]:
    positive_branches = [
        node
        for node in tree.iter_nodes()
        if node is not tree.root and node.branch_length is not None and node.branch_length > 0
    ]
    if len(positive_branches) < 2:
        return [], []

    baseline = _median(sorted(float(node.branch_length) for node in positive_branches if node.branch_length is not None))
    long_threshold = baseline * long_factor
    short_threshold = baseline * short_factor

    def item(node) -> BranchLengthOutlier:
        return BranchLengthOutlier(
            node=_node_signature(node),
            branch_length=round(float(node.branch_length), 15),
            branch_type="terminal" if node.is_leaf() else "internal",
        )

    long_outliers = sorted(
        [item(node) for node in positive_branches if float(node.branch_length) > long_threshold],
        key=lambda row: (-row.branch_length, row.node),
    )
    short_outliers = sorted(
        [item(node) for node in positive_branches if 0 < float(node.branch_length) < short_threshold],
        key=lambda row: (row.branch_length, row.node),
    )
    return long_outliers, short_outliers


def _star_like(tree: PhyloTree) -> bool:
    threshold = max(4, math.ceil(tree.tip_count * 0.75))
    for node in tree.iter_nodes():
        leaf_children = sum(1 for child in node.children if child.is_leaf())
        if leaf_children == len(node.children) and leaf_children >= threshold:
            return True
    return False


def _parse_internal_label_numeric(label: str) -> float | None:
    try:
        return float(label)
    except ValueError:
        return None


def _internal_label_diagnostics(
    tree: PhyloTree,
) -> tuple[list[InternalLabelInterpretation], list[InternalLabelInterpretation], list[str], bool]:
    likely_support_labels: list[InternalLabelInterpretation] = []
    likely_named_internal_labels: list[InternalLabelInterpretation] = []

    for node in tree.iter_nodes():
        if node.is_leaf() or node.name is None or not node.name.strip():
            continue
        numeric_value = _parse_internal_label_numeric(node.name.strip())
        interpretation = InternalLabelInterpretation(
            node=_node_signature(node),
            label=node.name,
            interpretation="support" if numeric_value is not None else "name",
            numeric_value=None if numeric_value is None else round(numeric_value, 15),
        )
        if numeric_value is not None:
            likely_support_labels.append(interpretation)
        else:
            likely_named_internal_labels.append(interpretation)

    suspicious_ranges: list[str] = []
    support_values = [row.numeric_value for row in likely_support_labels if row.numeric_value is not None]
    for row in likely_support_labels:
        if row.numeric_value is None:
            continue
        if row.numeric_value < 0:
            suspicious_ranges.append(f"support value {row.numeric_value:g} at {row.node} is negative")
        elif row.numeric_value > 100:
            suspicious_ranges.append(f"support value {row.numeric_value:g} at {row.node} exceeds 100")

    has_fraction_scale = any(0 <= value <= 1 for value in support_values)
    has_percent_scale = any(1 < value <= 100 for value in support_values)
    mixed_scales = has_fraction_scale and has_percent_scale
    return (
        sorted(likely_support_labels, key=lambda row: (row.node, row.label)),
        sorted(likely_named_internal_labels, key=lambda row: (row.node, row.label)),
        suspicious_ranges,
        mixed_scales,
    )


def _tree_quality_warnings(
    tree: PhyloTree,
    *,
    branch_length_status: str,
    zero_length_branch_count: int,
    negative_branch_count: int,
    polytomy_nodes: list[str],
    unusually_imbalanced: bool | None,
    long_branch_taxa: list[str],
    short_branch_outliers: list[BranchLengthOutlier],
    suspicious_support_value_ranges: list[str],
    mixed_support_scales: bool,
    star_like: bool,
    comb_like: bool,
) -> list[TreeQualityWarning]:
    warnings: list[TreeQualityWarning] = []
    if branch_length_status == "partial":
        warnings.append(
            TreeQualityWarning(
                code="partial_branch_lengths",
                message="tree has partial branch lengths, so weighted diagnostics are incomplete",
                penalty=20.0,
                affected_taxa=[],
                affected_nodes=[],
            )
        )
    if branch_length_status == "absent":
        warnings.append(
            TreeQualityWarning(
                code="missing_branch_lengths",
                message="tree has no branch lengths, so weighted diagnostics are unavailable",
                penalty=35.0,
                affected_taxa=[],
                affected_nodes=[],
            )
        )
    if negative_branch_count:
        warnings.append(
            TreeQualityWarning(
                code="negative_branch_lengths",
                message="tree contains negative branch lengths, which violate common phylogenetic assumptions",
                penalty=35.0,
                affected_taxa=[],
                affected_nodes=_negative_branch_nodes(tree),
            )
        )
    if zero_length_branch_count:
        warnings.append(
            TreeQualityWarning(
                code="zero_length_branches",
                message="tree contains zero-length branches that may represent unresolved or collapsed splits",
                penalty=10.0,
                affected_taxa=[],
                affected_nodes=_zero_length_branch_nodes(tree),
            )
        )
    if polytomy_nodes:
        warnings.append(
            TreeQualityWarning(
                code="polytomies",
                message="tree contains multifurcations that reduce binary-resolution diagnostics",
                penalty=10.0,
                affected_taxa=[],
                affected_nodes=polytomy_nodes,
            )
        )
    if unusually_imbalanced:
        warnings.append(
            TreeQualityWarning(
                code="unusually_imbalanced",
                message="tree shape is unusually imbalanced relative to a fully balanced binary tree",
                penalty=15.0,
                affected_taxa=[],
                affected_nodes=[_node_signature(tree.root)],
            )
        )
    if long_branch_taxa:
        warnings.append(
            TreeQualityWarning(
                code="long_branches",
                message="one or more taxa sit on unusually long terminal branches",
                penalty=15.0,
                affected_taxa=long_branch_taxa,
                affected_nodes=[],
            )
        )
    if short_branch_outliers:
        warnings.append(
            TreeQualityWarning(
                code="short_branches",
                message="tree contains unusually short nonzero branches that may represent numerically fragile splits",
                penalty=5.0,
                affected_taxa=[row.node for row in short_branch_outliers if row.branch_type == "terminal"],
                affected_nodes=[row.node for row in short_branch_outliers if row.branch_type == "internal"],
            )
        )
    if suspicious_support_value_ranges:
        warnings.append(
            TreeQualityWarning(
                code="suspicious_support_ranges",
                message="tree contains support-like internal labels outside standard probability or percentage ranges",
                penalty=10.0,
                affected_taxa=[],
                affected_nodes=[],
            )
        )
    if mixed_support_scales:
        warnings.append(
            TreeQualityWarning(
                code="mixed_support_scales",
                message="tree mixes support-like internal labels on probability and percentage scales",
                penalty=5.0,
                affected_taxa=[],
                affected_nodes=[],
            )
        )
    if star_like:
        warnings.append(
            TreeQualityWarning(
                code="star_like",
                message="tree has a star-like shape with one node connected directly to most tips",
                penalty=10.0,
                affected_taxa=[],
                affected_nodes=[_node_signature(tree.root)],
            )
        )
    if comb_like:
        warnings.append(
            TreeQualityWarning(
                code="comb_like",
                message="tree has an extreme comb-like topology with a single ladderized chain",
                penalty=10.0,
                affected_taxa=[],
                affected_nodes=[_node_signature(tree.root)],
            )
        )
    return warnings


def _tree_quality_score(warnings: list[TreeQualityWarning]) -> float:
    return max(0.0, round(100.0 - sum(warning.penalty for warning in warnings), 15))


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
    missing_internal_branch_nodes = _missing_internal_branch_nodes(tree)
    missing_terminal_branch_taxa = _missing_terminal_branch_taxa(tree)
    singleton_internal_nodes = _singleton_internal_nodes(tree)
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
    if missing_internal_branch_nodes:
        warnings.append("tree contains internal branches without lengths")
    if missing_terminal_branch_taxa:
        warnings.append("tree contains terminal branches without lengths")
    if singleton_internal_nodes:
        warnings.append("tree contains singleton internal nodes")
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
        missing_internal_branch_nodes=missing_internal_branch_nodes,
        missing_terminal_branch_taxa=missing_terminal_branch_taxa,
        singleton_internal_nodes=singleton_internal_nodes,
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
    internal_child_counts = _internal_node_child_counts(tree)
    singleton_internal_nodes = _singleton_internal_nodes(tree)
    missing_internal_branch_nodes = _missing_internal_branch_nodes(tree)
    missing_terminal_branch_taxa = _missing_terminal_branch_taxa(tree)
    depths = _leaf_depths(tree)
    zero_length_branch_count = sum(1 for length in branch_lengths if length == 0)
    ultrametric = _ultrametric(tree)
    branch_length_summary = _branch_length_summary(tree)
    colless_imbalance_index = _colless_imbalance_index(tree)
    normalized_colless_imbalance = _normalized_colless_imbalance(tree)
    unusually_imbalanced = _unusually_imbalanced(tree)
    long_branch_taxa = _long_branch_taxa(tree)
    long_branch_outliers, short_branch_outliers = _branch_outliers(tree)
    likely_support_labels, likely_named_internal_labels, suspicious_support_value_ranges, mixed_support_scales = (
        _internal_label_diagnostics(tree)
    )
    star_like = _star_like(tree)
    comb_like = _comb_like(tree)
    _, _, negative_branch_count = _branch_length_health(tree)
    tree_quality_warnings = _tree_quality_warnings(
        tree,
        branch_length_status=branch_length_status,
        zero_length_branch_count=zero_length_branch_count,
        negative_branch_count=negative_branch_count,
        polytomy_nodes=polytomy_nodes,
        unusually_imbalanced=unusually_imbalanced,
        long_branch_taxa=long_branch_taxa,
        short_branch_outliers=short_branch_outliers,
        suspicious_support_value_ranges=suspicious_support_value_ranges,
        mixed_support_scales=mixed_support_scales,
        star_like=star_like,
        comb_like=comb_like,
    )
    warnings: list[str] = []
    if zero_length_branch_count:
        warnings.append("tree contains zero-length branches")
    if missing_internal_branch_nodes:
        warnings.append("tree contains internal branches without lengths")
    if missing_terminal_branch_taxa:
        warnings.append("tree contains terminal branches without lengths")
    if singleton_internal_nodes:
        warnings.append("tree contains singleton internal nodes")
    if polytomy_nodes:
        warnings.append("tree contains one or more polytomies")
    if branch_length_status == "partial":
        warnings.append("tree contains partial branch lengths")
    if branch_length_status == "absent":
        warnings.append("tree contains no branch lengths")
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
        internal_child_counts=internal_child_counts,
        singleton_internal_nodes=singleton_internal_nodes,
        polytomy_count=len(polytomy_nodes),
        polytomy_nodes=polytomy_nodes,
        has_branch_lengths=any(length is not None for length in branch_lengths),
        branch_length_status=branch_length_status,
        missing_internal_branch_nodes=missing_internal_branch_nodes,
        missing_terminal_branch_taxa=missing_terminal_branch_taxa,
        is_ultrametric=ultrametric,
        total_branch_length=tree.total_branch_length(),
        branch_length_summary=branch_length_summary,
        tree_diameter=_tree_diameter(tree),
        zero_length_branch_count=zero_length_branch_count,
        min_root_to_tip=min(lengths) if lengths else None,
        max_root_to_tip=max(lengths) if lengths else None,
        max_depth=_max_depth(tree),
        mean_depth=sum(depths) / len(depths),
        colless_imbalance_index=colless_imbalance_index,
        normalized_colless_imbalance=normalized_colless_imbalance,
        sackin_imbalance_index=sum(depths),
        unusually_imbalanced=unusually_imbalanced,
        long_branch_taxa=long_branch_taxa,
        long_branch_outliers=long_branch_outliers,
        short_branch_outliers=short_branch_outliers,
        suspicious_support_value_ranges=suspicious_support_value_ranges,
        mixed_support_scales=mixed_support_scales,
        likely_support_labels=likely_support_labels,
        likely_named_internal_labels=likely_named_internal_labels,
        star_like=star_like,
        comb_like=comb_like,
        tree_quality_score=_tree_quality_score(tree_quality_warnings),
        tree_quality_warnings=tree_quality_warnings,
        imbalance_summary=_imbalance_summary(tree),
        cherry_count=_cherry_count(tree),
        taxa=sorted(tree.tip_names),
        warnings=warnings,
    )


def diagnose_tree_path(path: Path, *, source_format: str | None = None) -> TreeDiagnosticReport:
    """Return a combined inspection and validation report for one tree."""
    return TreeDiagnosticReport(
        path=path,
        inspection=inspect_tree_path(path, source_format=source_format),
        validation=validate_tree_path(path, source_format=source_format),
    )

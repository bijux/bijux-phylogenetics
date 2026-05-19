from __future__ import annotations

from collections import Counter
import math
from pathlib import Path
from statistics import mean, median

from bijux_phylogenetics.core.taxonomy import (
    inspect_tree_taxa_safety,
    inspect_tree_taxon_identity,
)
from bijux_phylogenetics.core.tree import PhyloTree
from bijux_phylogenetics.core.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.runtime.errors import (
    DuplicateTaxonError,
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    UnnamedTipError,
    UnrootedTreeError,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.trees import summarize_tree_shape_from_tree
from .models import (
    BranchLengthContextAssessment,
    BranchLengthOutlier,
    BranchLengthRepairSuggestion,
    BranchLengthSummary,
    InternalLabelInterpretation,
    InternalNodeChildCount,
    InternalNodeLabelConflict,
    RootStateConfidenceReport,
    StableNodeIdentity,
    TreeDiagnosticReport,
    TreeFinding,
    TreeForensicReport,
    TreeInspectionReport,
    TreeIntegrityIssue,
    TreeQualityWarning,
    TreeValidationReport,
    UnsafeExternalLabel,
)

LONG_BRANCH_OUTLIER_FACTOR = 3.0
SHORT_BRANCH_OUTLIER_FACTOR = 0.1
TREE_IMBALANCE_WARNING_THRESHOLD = 0.75
STAR_LIKE_FRACTION_THRESHOLD = 0.75
ROOT_BRANCH_LENGTH_IMBALANCE_THRESHOLD = 10.0
ROOT_CHILD_BALANCE_RATIO_THRESHOLD = 0.1


def _load_tree(path: Path, *, source_format: str | None = None) -> PhyloTree:
    return load_tree(path, source_format=source_format)


def _count_polytomies(tree: PhyloTree) -> int:
    return sum(
        1 for node in tree.iter_nodes() if not node.is_leaf() and len(node.children) > 2
    )


def _integrity_issues(tree: PhyloTree) -> list[TreeIntegrityIssue]:
    issues: list[TreeIntegrityIssue] = []
    seen: set[int] = set()
    stack: set[int] = set()
    parent_counts: Counter[int] = Counter()
    node_ids: dict[int, str] = {}

    def visit(node) -> None:
        identifier = id(node)
        node_ids.setdefault(identifier, _node_signature(node))
        if identifier in stack:
            issues.append(
                TreeIntegrityIssue(
                    code="cycle",
                    message="tree contains a cycle in parent-child traversal",
                    severity="fatal",
                    affected_nodes=[node_ids[identifier]],
                )
            )
            return
        if identifier in seen:
            parent_counts[identifier] += 1
            return
        seen.add(identifier)
        stack.add(identifier)
        for child in node.children:
            if child is node:
                issues.append(
                    TreeIntegrityIssue(
                        code="self_child",
                        message="tree contains a node that references itself as a child",
                        severity="fatal",
                        affected_nodes=[_node_signature(node)],
                    )
                )
                continue
            parent_counts[id(child)] += 1
            visit(child)
        stack.remove(identifier)

    visit(tree.root)
    for identifier, count in sorted(
        parent_counts.items(), key=lambda item: (node_ids.get(item[0], ""), item[1])
    ):
        if count > 1:
            issues.append(
                TreeIntegrityIssue(
                    code="duplicate_parentage",
                    message="tree contains a node referenced by more than one parent",
                    severity="fatal",
                    affected_nodes=[node_ids.get(identifier, "<unknown>")],
                )
            )
    if tree.tip_count == 0:
        issues.append(
            TreeIntegrityIssue(
                code="no_tips",
                message="tree does not contain any terminal taxa",
                severity="fatal",
                affected_nodes=[_node_signature(tree.root)],
            )
        )
    if not tree.root.children and tree.root.name is None:
        issues.append(
            TreeIntegrityIssue(
                code="empty_root",
                message="tree root is empty and has no descendants",
                severity="fatal",
                affected_nodes=["<unnamed>"],
            )
        )
    if len(tree.root.children) == 1:
        issues.append(
            TreeIntegrityIssue(
                code="degenerate_root",
                message="tree root has only one child and represents a degenerate rooted structure",
                severity="blocker",
                affected_nodes=[_node_signature(tree.root)],
            )
        )
    return issues


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


def _stable_node_identities(tree: PhyloTree) -> list[StableNodeIdentity]:
    identities: list[StableNodeIdentity] = []
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = _descendant_taxa(node)
        identities.append(
            StableNodeIdentity(
                node_id=_node_signature(node),
                descendant_taxa=descendant_taxa,
            )
        )
    return sorted(identities, key=lambda row: (len(row.descendant_taxa), row.node_id))


def _root_state_confidence(tree: PhyloTree) -> RootStateConfidenceReport:
    rationale: list[str] = []
    suspicious_reasons: list[str] = []
    explicit_rooted = tree.rooted is True
    apparent_rooted = len(tree.root.children) == 2
    if explicit_rooted:
        classification = "explicitly_rooted"
        rationale.append("source parser preserved an explicit rooted flag")
    elif apparent_rooted:
        classification = "apparently_rooted"
        rationale.append(
            "root has exactly two child clades but no explicit rooted flag was preserved"
        )
    elif len(tree.root.children) > 2:
        classification = "unrooted"
        rationale.append(
            "root has more than two child clades and behaves like an unrooted representation"
        )
    else:
        classification = "ambiguous"
        rationale.append(
            "root structure is degenerate or not informative enough to classify confidently"
        )

    if apparent_rooted and tree.tip_count >= 4:
        child_tip_counts = sorted(
            len(_descendant_taxa(child)) for child in tree.root.children
        )
        if (
            child_tip_counts
            and child_tip_counts[0] == 1
            and child_tip_counts[-1] >= max(3, tree.tip_count - 1)
        ):
            suspicious_reasons.append(
                "root isolates a single tip against almost the entire tree"
            )
        if len(child_tip_counts) == 2 and child_tip_counts[1] > 0:
            imbalance_ratio = child_tip_counts[0] / child_tip_counts[1]
            if imbalance_ratio <= ROOT_CHILD_BALANCE_RATIO_THRESHOLD:
                suspicious_reasons.append(
                    "root creates an extreme basal imbalance between its two child clades"
                )
        child_lengths = [
            float(child.branch_length)
            for child in tree.root.children
            if child.branch_length is not None and child.branch_length > 0
        ]
        if len(child_lengths) == 2:
            shorter = min(child_lengths)
            longer = max(child_lengths)
            if (
                shorter > 0
                and longer / shorter >= ROOT_BRANCH_LENGTH_IMBALANCE_THRESHOLD
            ):
                suspicious_reasons.append(
                    "one basal branch is more than ten times longer than its sister branch"
                )
    return RootStateConfidenceReport(
        classification=classification,
        rationale=rationale,
        suspicious_placement=bool(suspicious_reasons),
        suspicious_reasons=suspicious_reasons,
    )


def _edge_count(tree: PhyloTree) -> int:
    return sum(1 for node in tree.iter_nodes() if node is not tree.root)


def _node_count(tree: PhyloTree) -> int:
    return sum(1 for _ in tree.iter_nodes())


def _branch_length_health(tree: PhyloTree) -> tuple[bool, int, int]:
    branch_lengths = [
        node.branch_length for node in tree.iter_nodes() if node is not tree.root
    ]
    has_complete = (
        all(length is not None for length in branch_lengths)
        if branch_lengths
        else False
    )
    zero_count = sum(1 for length in branch_lengths if length == 0)
    negative_count = sum(
        1 for length in branch_lengths if length is not None and length < 0
    )
    return has_complete, zero_count, negative_count


def _internal_node_child_counts(tree: PhyloTree) -> list[InternalNodeChildCount]:
    return [
        InternalNodeChildCount(
            node=_node_signature(node), child_count=len(node.children)
        )
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
        if node is not tree.root
        and node.branch_length is not None
        and node.branch_length < 0
    )


def _branch_length_status(tree: PhyloTree) -> str:
    branch_lengths = [
        node.branch_length for node in tree.iter_nodes() if node is not tree.root
    ]
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
        adjacency.setdefault(parent_id, []).append(
            (child_id, float(child.branch_length or 0.0))
        )
        adjacency.setdefault(child_id, []).append(
            (parent_id, float(child.branch_length or 0.0))
        )
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


def _ultrametric(
    tree: PhyloTree,
    *,
    tolerance: float = APE_ULTRAMETRIC_TOLERANCE,
) -> bool | None:
    lengths = tree.root_to_tip_lengths()
    if not lengths or any(length is None for length in lengths):
        return None
    return summarize_ultrametric_tip_depths(
        {
            tip_name or f"unnamed-tip-{index}": float(length)
            for index, (tip_name, length) in enumerate(
                tree.root_to_tip_pairs(),
                start=1,
            )
            if length is not None
        },
        tolerance=tolerance,
    ).ultrametric


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
    return round(raw_score / maximum, 15)


def _unusually_imbalanced(
    tree: PhyloTree, *, threshold: float = TREE_IMBALANCE_WARNING_THRESHOLD
) -> bool | None:
    normalized = _normalized_colless_imbalance(tree)
    if normalized is None:
        return None
    return normalized >= threshold


def _comb_like(tree: PhyloTree) -> bool:
    if len(tree.root.children) != 2:
        return False
    if any(
        not node.is_leaf() and len(node.children) != 2 for node in tree.iter_nodes()
    ):
        return False
    return _cherry_count(tree) == 1 and _max_depth(tree) == tree.tip_count - 1


def _long_branch_taxa(
    tree: PhyloTree, *, factor: float = LONG_BRANCH_OUTLIER_FACTOR
) -> list[str]:
    terminal_lengths = [
        (name, length)
        for name, length in tree.terminal_branch_lengths()
        if length is not None and length > 0
    ]
    if len(terminal_lengths) < 2:
        return []
    baseline = _median(sorted(length for _, length in terminal_lengths))
    threshold = baseline * factor
    return sorted(name for name, length in terminal_lengths if length > threshold)


def _branch_outliers(
    tree: PhyloTree,
    *,
    long_factor: float = LONG_BRANCH_OUTLIER_FACTOR,
    short_factor: float = SHORT_BRANCH_OUTLIER_FACTOR,
) -> tuple[list[BranchLengthOutlier], list[BranchLengthOutlier]]:
    positive_branches = [
        node
        for node in tree.iter_nodes()
        if node is not tree.root
        and node.branch_length is not None
        and node.branch_length > 0
    ]
    if len(positive_branches) < 2:
        return [], []

    baseline = _median(
        sorted(
            float(node.branch_length)
            for node in positive_branches
            if node.branch_length is not None
        )
    )
    long_threshold = baseline * long_factor
    short_threshold = baseline * short_factor

    def item(node) -> BranchLengthOutlier:
        return BranchLengthOutlier(
            node=_node_signature(node),
            branch_length=round(float(node.branch_length), 15),
            branch_type="terminal" if node.is_leaf() else "internal",
        )

    long_outliers = sorted(
        [
            item(node)
            for node in positive_branches
            if float(node.branch_length) > long_threshold
        ],
        key=lambda row: (-row.branch_length, row.node),
    )
    short_outliers = sorted(
        [
            item(node)
            for node in positive_branches
            if 0 < float(node.branch_length) < short_threshold
        ],
        key=lambda row: (row.branch_length, row.node),
    )
    return long_outliers, short_outliers


def _star_like(tree: PhyloTree) -> bool:
    threshold = max(4, math.ceil(tree.tip_count * STAR_LIKE_FRACTION_THRESHOLD))
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
) -> tuple[
    list[InternalLabelInterpretation],
    list[InternalLabelInterpretation],
    list[str],
    bool,
]:
    likely_support_labels: list[InternalLabelInterpretation] = []
    likely_named_internal_labels: list[InternalLabelInterpretation] = []

    for node in tree.iter_nodes():
        if node.is_leaf() or node.name is None or not node.name.strip():
            continue
        stripped = node.name.strip()
        numeric_value = _parse_internal_label_numeric(stripped)
        if numeric_value is None:
            interpretation = (
                "ambiguous_alphanumeric_label"
                if any(character.isdigit() for character in stripped)
                else "named_internal_label"
            )
        elif numeric_value < 0 or numeric_value > 100:
            interpretation = "out_of_range_support"
        elif 0.0 <= numeric_value <= 1.0:
            interpretation = "fractional_support"
        else:
            interpretation = "percentage_support"
        interpretation = InternalLabelInterpretation(
            node=_node_signature(node),
            node_id=_node_signature(node),
            label=node.name,
            interpretation=interpretation,
            numeric_value=None if numeric_value is None else round(numeric_value, 15),
        )
        if numeric_value is not None:
            likely_support_labels.append(interpretation)
        else:
            likely_named_internal_labels.append(interpretation)

    suspicious_ranges: list[str] = []
    support_values = [
        row.numeric_value
        for row in likely_support_labels
        if row.numeric_value is not None
    ]
    for row in likely_support_labels:
        if row.numeric_value is None:
            continue
        if row.numeric_value < 0:
            suspicious_ranges.append(
                f"support value {row.numeric_value:g} at {row.node} is negative"
            )
        elif row.numeric_value > 100:
            suspicious_ranges.append(
                f"support value {row.numeric_value:g} at {row.node} exceeds 100"
            )

    has_fraction_scale = any(0 <= value <= 1 for value in support_values)
    has_percent_scale = any(1 < value <= 100 for value in support_values)
    mixed_scales = has_fraction_scale and has_percent_scale
    return (
        sorted(likely_support_labels, key=lambda row: (row.node, row.label)),
        sorted(likely_named_internal_labels, key=lambda row: (row.node, row.label)),
        suspicious_ranges,
        mixed_scales,
    )


def _internal_label_conflicts(
    likely_support_labels: list[InternalLabelInterpretation],
    likely_named_internal_labels: list[InternalLabelInterpretation],
    suspicious_support_value_ranges: list[str],
    mixed_support_scales: bool,
) -> list[InternalNodeLabelConflict]:
    conflicts: list[InternalNodeLabelConflict] = []
    for row in likely_support_labels:
        if row.numeric_value is None:
            continue
        if row.numeric_value < 0 or row.numeric_value > 100:
            conflicts.append(
                InternalNodeLabelConflict(
                    node_id=row.node_id,
                    label=row.label,
                    conflict_type="support_out_of_range",
                    detail="numeric internal label falls outside accepted support scales",
                )
            )
        elif 0.0 <= row.numeric_value <= 1.0:
            conflicts.append(
                InternalNodeLabelConflict(
                    node_id=row.node_id,
                    label=row.label,
                    conflict_type="ambiguous_fraction_support",
                    detail="numeric internal label could be posterior support or a named code on a 0-1 scale",
                )
            )
    for row in likely_named_internal_labels:
        stripped = row.label.strip()
        if any(character.isdigit() for character in stripped):
            conflicts.append(
                InternalNodeLabelConflict(
                    node_id=row.node_id,
                    label=row.label,
                    conflict_type="alphanumeric_internal_label",
                    detail="non-numeric internal label contains digits and may conflate clade naming with support annotation",
                )
            )
    if mixed_support_scales:
        conflicts.append(
            InternalNodeLabelConflict(
                node_id="<tree>",
                label="mixed-scales",
                conflict_type="mixed_support_scales",
                detail="tree mixes fraction-like and percentage-like support labels",
            )
        )
    if suspicious_support_value_ranges:
        conflicts.append(
            InternalNodeLabelConflict(
                node_id="<tree>",
                label="support-range",
                conflict_type="support_range_summary",
                detail="one or more internal support-like labels fall outside standard ranges",
            )
        )
    return sorted(
        conflicts, key=lambda row: (row.node_id, row.conflict_type, row.label)
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
                affected_taxa=[
                    row.node
                    for row in short_branch_outliers
                    if row.branch_type == "terminal"
                ],
                affected_nodes=[
                    row.node
                    for row in short_branch_outliers
                    if row.branch_type == "internal"
                ],
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


def _unsafe_external_labels(tree: PhyloTree) -> list[UnsafeExternalLabel]:
    report = inspect_tree_taxa_safety(tree, policy="spaces-to-underscores")
    labels: list[UnsafeExternalLabel] = []
    default_engines = ["iqtree", "raxml", "mrbayes", "beast", "r", "shell"]
    for entry in report.unsafe_taxa:
        engines = list(default_engines)
        labels.append(
            UnsafeExternalLabel(
                raw_label=entry.raw_label,
                normalized_label=entry.normalized_label,
                engines=engines,
                reasons=entry.reasons,
            )
        )
    for collision in report.collisions:
        for raw_label in collision.raw_labels:
            if raw_label in {item.raw_label for item in labels}:
                continue
            labels.append(
                UnsafeExternalLabel(
                    raw_label=raw_label,
                    normalized_label=collision.normalized_label,
                    engines=list(default_engines),
                    reasons=["collides with another label after normalization"],
                )
            )
    return sorted(labels, key=lambda row: row.raw_label)


def _branch_length_contexts(
    inspection: TreeInspectionReport,
) -> list[BranchLengthContextAssessment]:
    contexts: list[BranchLengthContextAssessment] = []
    zero_length = inspection.zero_length_branch_count > 0
    negative = any(
        warning.code == "negative_branch_lengths"
        for warning in inspection.tree_quality_warnings
    )
    long_outliers = bool(inspection.long_branch_outliers)

    contexts.append(
        BranchLengthContextAssessment(
            context="topology_only",
            allowed=True,
            blocked_by=[],
            warnings=[
                "topology-only workflows ignore branch lengths, but unresolved splits still affect interpretation"
                if inspection.polytomy_count
                else ""
            ],
        )
    )

    substitution_blockers: list[str] = []
    substitution_warnings: list[str] = []
    if inspection.branch_length_status != "complete":
        substitution_blockers.append("complete branch lengths are required")
    if negative:
        substitution_blockers.append(
            "negative branch lengths violate substitution-distance assumptions"
        )
    if zero_length:
        substitution_warnings.append(
            "zero-length branches may represent collapsed or unresolved substitutions"
        )
    if long_outliers:
        substitution_warnings.append(
            "extreme branch-length outliers may indicate saturation or model mismatch"
        )
    contexts.append(
        BranchLengthContextAssessment(
            context="substitution_tree",
            allowed=not substitution_blockers,
            blocked_by=substitution_blockers,
            warnings=substitution_warnings,
        )
    )

    time_blockers = list(substitution_blockers)
    time_warnings = list(substitution_warnings)
    if inspection.is_ultrametric is not True:
        time_blockers.append("time trees require ultrametric root-to-tip distances")
    contexts.append(
        BranchLengthContextAssessment(
            context="time_tree",
            allowed=not time_blockers,
            blocked_by=time_blockers,
            warnings=time_warnings,
        )
    )

    comparative_blockers = list(substitution_blockers)
    comparative_warnings = list(substitution_warnings)
    if not inspection.rooted:
        comparative_blockers.append("comparative methods require a rooted tree")
    if inspection.polytomy_count:
        comparative_warnings.append(
            "polytomies reduce comparative-model identifiability"
        )
    contexts.append(
        BranchLengthContextAssessment(
            context="comparative_methods",
            allowed=not comparative_blockers,
            blocked_by=comparative_blockers,
            warnings=comparative_warnings,
        )
    )
    for context in contexts:
        context.warnings = [warning for warning in context.warnings if warning]
    return contexts


def _branch_length_repair_suggestions(
    inspection: TreeInspectionReport,
) -> list[BranchLengthRepairSuggestion]:
    suggestions: list[BranchLengthRepairSuggestion] = []
    if inspection.branch_length_status == "absent":
        suggestions.append(
            BranchLengthRepairSuggestion(
                issue_code="missing_branch_lengths",
                summary="branch lengths are absent throughout the tree",
                blocked_analyses=[
                    "substitution_tree",
                    "time_tree",
                    "comparative_methods",
                ],
                suggested_action="re-export the inference tree with branch lengths or estimate them before weighted analyses",
            )
        )
    if inspection.branch_length_status == "partial":
        suggestions.append(
            BranchLengthRepairSuggestion(
                issue_code="partial_branch_lengths",
                summary="some branches are missing lengths",
                blocked_analyses=[
                    "substitution_tree",
                    "time_tree",
                    "comparative_methods",
                ],
                suggested_action="fill or re-estimate missing branch lengths consistently; partial lengths are unsafe for weighted methods",
            )
        )
    if inspection.zero_length_branch_count:
        suggestions.append(
            BranchLengthRepairSuggestion(
                issue_code="zero_length_branches",
                summary="tree contains zero-length branches",
                blocked_analyses=["time_tree"],
                suggested_action="confirm whether zero lengths indicate unresolved splits, then collapse or re-estimate those branches",
            )
        )
    negative_nodes = [
        warning
        for warning in inspection.tree_quality_warnings
        if warning.code == "negative_branch_lengths"
    ]
    if negative_nodes:
        suggestions.append(
            BranchLengthRepairSuggestion(
                issue_code="negative_branch_lengths",
                summary="tree contains negative branch lengths",
                blocked_analyses=[
                    "substitution_tree",
                    "time_tree",
                    "comparative_methods",
                    "publication",
                ],
                suggested_action="recompute the tree or sanitize the source export; negative branch lengths should not be silently truncated",
            )
        )
    if inspection.long_branch_outliers:
        suggestions.append(
            BranchLengthRepairSuggestion(
                issue_code="extreme_branch_lengths",
                summary="tree contains branch-length outliers",
                blocked_analyses=[],
                suggested_action="review alignment quality, saturation, and model fit before interpreting long-branch placement",
            )
        )
    return suggestions


def _findings_from_reports(
    integrity_issues: list[TreeIntegrityIssue],
    inspection: TreeInspectionReport,
    duplicate_taxa: list[str],
    missing_taxa: int,
    root_state: RootStateConfidenceReport,
    internal_label_conflicts: list[InternalNodeLabelConflict],
    unsafe_labels: list[UnsafeExternalLabel],
) -> list[TreeFinding]:
    findings: list[TreeFinding] = [
        TreeFinding(
            code=issue.code,
            message=issue.message,
            severity=issue.severity,
            affected_taxa=[],
            affected_nodes=issue.affected_nodes,
        )
        for issue in integrity_issues
    ]
    if duplicate_taxa:
        findings.append(
            TreeFinding(
                code="duplicate_taxa",
                message="tree contains duplicate tip labels and is biologically unsafe until taxa are disambiguated",
                severity="blocker",
                affected_taxa=duplicate_taxa,
                affected_nodes=[],
            )
        )
    if missing_taxa:
        findings.append(
            TreeFinding(
                code="unnamed_tips",
                message="tree contains unnamed tips and cannot be safely reconciled to external data",
                severity="blocker",
                affected_taxa=[],
                affected_nodes=[],
            )
        )
    for warning in inspection.tree_quality_warnings:
        severity = "warning"
        if warning.code in {
            "negative_branch_lengths",
            "missing_branch_lengths",
            "partial_branch_lengths",
        }:
            severity = "blocker"
        findings.append(
            TreeFinding(
                code=warning.code,
                message=warning.message,
                severity=severity,
                affected_taxa=warning.affected_taxa,
                affected_nodes=warning.affected_nodes,
            )
        )
    if root_state.suspicious_placement:
        findings.append(
            TreeFinding(
                code="suspicious_root_placement",
                message="root placement appears biologically suspicious",
                severity="warning",
                affected_taxa=[],
                affected_nodes=["<root>"],
            )
        )
    for conflict in internal_label_conflicts:
        findings.append(
            TreeFinding(
                code=conflict.conflict_type,
                message=conflict.detail,
                severity="warning",
                affected_taxa=[],
                affected_nodes=[]
                if conflict.node_id == "<tree>"
                else [conflict.node_id],
            )
        )
    if unsafe_labels:
        findings.append(
            TreeFinding(
                code="unsafe_external_labels",
                message="one or more taxon labels are unsafe across common phylogenetics engines or shell workflows",
                severity="warning",
                affected_taxa=[row.raw_label for row in unsafe_labels],
                affected_nodes=[],
            )
        )
    if (
        inspection.taxon_identity_audit.whitespace_variants
        or inspection.taxon_identity_audit.underscore_space_collisions
        or inspection.taxon_identity_audit.case_collisions
        or inspection.taxon_identity_audit.suspicious_near_duplicates
    ):
        findings.append(
            TreeFinding(
                code="taxon_identity_conflicts",
                message="one or more taxon labels are suspiciously similar and may not represent distinct biological identities",
                severity="warning",
                affected_taxa=sorted(
                    {
                        pair.left_label
                        for pair in (
                            inspection.taxon_identity_audit.whitespace_variants
                            + inspection.taxon_identity_audit.underscore_space_collisions
                            + inspection.taxon_identity_audit.case_collisions
                            + inspection.taxon_identity_audit.suspicious_near_duplicates
                        )
                    }
                    | {
                        pair.right_label
                        for pair in (
                            inspection.taxon_identity_audit.whitespace_variants
                            + inspection.taxon_identity_audit.underscore_space_collisions
                            + inspection.taxon_identity_audit.case_collisions
                            + inspection.taxon_identity_audit.suspicious_near_duplicates
                        )
                    }
                ),
                affected_nodes=[],
            )
        )
    severity_order = {"fatal": 0, "blocker": 1, "warning": 2, "info": 3}
    return sorted(
        findings,
        key=lambda row: (severity_order.get(row.severity, 99), row.code, row.message),
    )


def _validity_decision(findings: list[TreeFinding]) -> tuple[bool, bool, str]:
    has_fatal = any(finding.severity == "fatal" for finding in findings)
    has_blocker = any(finding.severity == "blocker" for finding in findings)
    has_warning = any(finding.severity == "warning" for finding in findings)
    syntax_valid = not has_fatal
    biologically_safe = syntax_valid and not has_blocker
    if not syntax_valid or has_blocker:
        return syntax_valid, biologically_safe, "invalid"
    if has_warning:
        return syntax_valid, biologically_safe, "valid_with_warnings"
    return syntax_valid, biologically_safe, "valid"


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
    inspection = inspect_tree_path(path, source_format=source_format)
    rooted = inspection.rooted
    has_complete, zero_count, negative_count = _branch_length_health(tree)
    branch_length_status = _branch_length_status(tree)
    missing_internal_branch_nodes = _missing_internal_branch_nodes(tree)
    missing_terminal_branch_taxa = _missing_terminal_branch_taxa(tree)
    singleton_internal_nodes = _singleton_internal_nodes(tree)
    missing_taxa, duplicate_taxa = _duplicate_taxa(tree)
    integrity_issues = _integrity_issues(tree)
    if duplicate_taxa and not allow_duplicates:
        raise DuplicateTaxonError(
            f"duplicate tip labels found: {', '.join(duplicate_taxa)}"
        )
    if missing_taxa and strict:
        raise UnnamedTipError(f"tree contains {missing_taxa} unnamed tip labels")
    if negative_count and not allow_negative_branch_lengths:
        raise InvalidBranchLengthError(
            f"tree contains {negative_count} negative branch lengths"
        )
    ultrametric = _ultrametric(tree)
    if require_rooted and not rooted:
        raise UnrootedTreeError(f"tree is not rooted: {path}")
    if require_ultrametric and ultrametric is not True:
        raise NonUltrametricTreeError(
            f"tree is not ultrametric within default validation tolerance: {path}"
        )
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
    if inspection.root_state_confidence.suspicious_placement:
        warnings.append("tree root placement appears biologically suspicious")
    if inspection.unsafe_external_labels:
        warnings.append(
            "tree contains taxon labels unsafe across common external engines"
        )
    if (
        inspection.taxon_identity_audit.whitespace_variants
        or inspection.taxon_identity_audit.underscore_space_collisions
        or inspection.taxon_identity_audit.case_collisions
        or inspection.taxon_identity_audit.suspicious_near_duplicates
    ):
        warnings.append("tree contains potentially ambiguous taxon identity variants")
    findings = _findings_from_reports(
        integrity_issues,
        inspection,
        duplicate_taxa,
        missing_taxa,
        inspection.root_state_confidence,
        inspection.internal_label_conflicts,
        inspection.unsafe_external_labels,
    )
    syntax_valid, biologically_safe, validity_decision = _validity_decision(findings)
    return TreeValidationReport(
        path=path,
        source_format=tree.source_format,
        tip_count=tree.tip_count,
        internal_node_count=tree.internal_node_count,
        rooted=rooted,
        syntax_valid=syntax_valid,
        biologically_safe=biologically_safe,
        validity_decision=validity_decision,
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
        integrity_issues=integrity_issues,
        warning_details=findings,
        root_state_confidence=inspection.root_state_confidence,
        branch_length_contexts=_branch_length_contexts(inspection),
        branch_length_repair_suggestions=_branch_length_repair_suggestions(inspection),
        internal_label_conflicts=inspection.internal_label_conflicts,
        stable_node_identities=inspection.stable_node_identities,
        unsafe_external_labels=inspection.unsafe_external_labels,
        taxon_identity_audit=inspection.taxon_identity_audit,
        warnings=warnings,
    )


def inspect_tree_path(
    path: Path, *, source_format: str | None = None
) -> TreeInspectionReport:
    """Inspect a tree file and return lightweight summary metrics."""
    tree = _load_tree(path, source_format=source_format)
    shape = summarize_tree_shape_from_tree(tree, source_path=path)
    lengths = [length for length in tree.root_to_tip_lengths() if length is not None]
    branch_lengths = [
        node.branch_length for node in tree.iter_nodes() if node is not tree.root
    ]
    polytomy_nodes = _polytomy_nodes(tree)
    branch_length_status = _branch_length_status(tree)
    internal_child_counts = _internal_node_child_counts(tree)
    singleton_internal_nodes = _singleton_internal_nodes(tree)
    missing_internal_branch_nodes = _missing_internal_branch_nodes(tree)
    missing_terminal_branch_taxa = _missing_terminal_branch_taxa(tree)
    zero_length_branch_count = sum(1 for length in branch_lengths if length == 0)
    ultrametric = _ultrametric(tree)
    branch_length_summary = _branch_length_summary(tree)
    long_branch_taxa = _long_branch_taxa(tree)
    long_branch_outliers, short_branch_outliers = _branch_outliers(tree)
    (
        likely_support_labels,
        likely_named_internal_labels,
        suspicious_support_value_ranges,
        mixed_support_scales,
    ) = _internal_label_diagnostics(tree)
    internal_label_conflicts = _internal_label_conflicts(
        likely_support_labels,
        likely_named_internal_labels,
        suspicious_support_value_ranges,
        mixed_support_scales,
    )
    root_state_confidence = _root_state_confidence(tree)
    stable_node_identities = _stable_node_identities(tree)
    unsafe_external_labels = _unsafe_external_labels(tree)
    taxon_identity_audit = inspect_tree_taxon_identity(tree)
    _, _, negative_branch_count = _branch_length_health(tree)
    tree_quality_warnings = _tree_quality_warnings(
        tree,
        branch_length_status=branch_length_status,
        zero_length_branch_count=zero_length_branch_count,
        negative_branch_count=negative_branch_count,
        polytomy_nodes=polytomy_nodes,
        unusually_imbalanced=shape.unusually_imbalanced,
        long_branch_taxa=long_branch_taxa,
        short_branch_outliers=short_branch_outliers,
        suspicious_support_value_ranges=suspicious_support_value_ranges,
        mixed_support_scales=mixed_support_scales,
        star_like=shape.star_like,
        comb_like=shape.comb_like,
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
    if (
        taxon_identity_audit.whitespace_variants
        or taxon_identity_audit.underscore_space_collisions
    ):
        warnings.append(
            "tree contains taxon labels with whitespace or underscore identity collisions"
        )
    if (
        taxon_identity_audit.case_collisions
        or taxon_identity_audit.suspicious_near_duplicates
    ):
        warnings.append(
            "tree contains potentially ambiguous near-duplicate taxon labels"
        )
    return TreeInspectionReport(
        path=path,
        source_format=tree.source_format,
        tip_count=tree.tip_count,
        node_count=_node_count(tree),
        internal_node_count=tree.internal_node_count,
        edge_count=_edge_count(tree),
        clade_count=tree.internal_node_count,
        rooted=root_state_confidence.classification
        in {"explicitly_rooted", "apparently_rooted"},
        root_state_confidence=root_state_confidence,
        is_binary=all(
            node.is_leaf() or len(node.children) == 2 for node in tree.iter_nodes()
        ),
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
        max_depth=shape.tree_height_edges,
        mean_depth=shape.mean_tip_depth_edges,
        colless_imbalance_index=shape.colless_imbalance_index,
        normalized_colless_imbalance=shape.normalized_colless_imbalance,
        sackin_imbalance_index=shape.sackin_imbalance_index,
        unusually_imbalanced=shape.unusually_imbalanced,
        long_branch_taxa=long_branch_taxa,
        long_branch_outliers=long_branch_outliers,
        short_branch_outliers=short_branch_outliers,
        suspicious_support_value_ranges=suspicious_support_value_ranges,
        mixed_support_scales=mixed_support_scales,
        likely_support_labels=likely_support_labels,
        likely_named_internal_labels=likely_named_internal_labels,
        internal_label_conflicts=internal_label_conflicts,
        stable_node_identities=stable_node_identities,
        unsafe_external_labels=unsafe_external_labels,
        taxon_identity_audit=taxon_identity_audit,
        star_like=shape.star_like,
        comb_like=shape.comb_like,
        tree_quality_score=_tree_quality_score(tree_quality_warnings),
        tree_quality_warnings=tree_quality_warnings,
        imbalance_summary=shape.imbalance_summary,
        cherry_count=shape.cherry_count,
        taxa=sorted(tree.tip_names),
        warnings=warnings,
    )


def forensic_tree_path(
    path: Path, *, source_format: str | None = None
) -> TreeForensicReport:
    """Build a reviewer-facing forensic summary of whether a tree is safe for downstream use."""
    inspection = inspect_tree_path(path, source_format=source_format)
    validation = validate_tree_path(path, source_format=source_format)
    context_lookup = {
        context.context: context for context in validation.branch_length_contexts
    }
    safe_for_topology_comparison = (
        validation.syntax_valid
        and not validation.duplicate_taxa
        and validation.missing_taxa == 0
    )
    safe_for_time_tree_analysis = (
        context_lookup["time_tree"].allowed and validation.biologically_safe
    )
    safe_for_comparative_methods = (
        context_lookup["comparative_methods"].allowed and validation.biologically_safe
    )
    safe_for_visualization = validation.syntax_valid
    safe_for_publication = (
        validation.biologically_safe and not inspection.internal_label_conflicts
    )
    warnings = list(dict.fromkeys([*validation.warnings, *inspection.warnings]))
    return TreeForensicReport(
        path=path,
        source_format=validation.source_format,
        syntax_valid=validation.syntax_valid,
        biologically_safe=validation.biologically_safe,
        validity_decision=validation.validity_decision,
        integrity_issues=validation.integrity_issues,
        findings=validation.warning_details,
        root_state_confidence=validation.root_state_confidence,
        branch_length_contexts=validation.branch_length_contexts,
        branch_length_repair_suggestions=validation.branch_length_repair_suggestions,
        internal_label_conflicts=validation.internal_label_conflicts,
        stable_node_identities=validation.stable_node_identities,
        unsafe_external_labels=validation.unsafe_external_labels,
        taxon_identity_audit=validation.taxon_identity_audit,
        safe_for_topology_comparison=safe_for_topology_comparison,
        safe_for_time_tree_analysis=safe_for_time_tree_analysis,
        safe_for_comparative_methods=safe_for_comparative_methods,
        safe_for_visualization=safe_for_visualization,
        safe_for_publication=safe_for_publication,
        warnings=warnings,
    )


def diagnose_tree_path(
    path: Path, *, source_format: str | None = None
) -> TreeDiagnosticReport:
    """Return a combined inspection and validation report for one tree."""
    return TreeDiagnosticReport(
        path=path,
        inspection=inspect_tree_path(path, source_format=source_format),
        validation=validate_tree_path(path, source_format=source_format),
        forensic=forensic_tree_path(path, source_format=source_format),
    )

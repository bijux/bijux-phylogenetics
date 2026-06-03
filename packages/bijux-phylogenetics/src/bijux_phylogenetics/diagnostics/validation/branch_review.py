from __future__ import annotations

from statistics import mean, median

from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import (
    BranchLengthContextAssessment,
    BranchLengthOutlier,
    BranchLengthRepairSuggestion,
    BranchLengthSummary,
    InternalNodeChildCount,
    TreeInspectionReport,
    TreeQualityWarning,
)
from .structure import _node_signature

LONG_BRANCH_OUTLIER_FACTOR = 3.0
SHORT_BRANCH_OUTLIER_FACTOR = 0.1


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

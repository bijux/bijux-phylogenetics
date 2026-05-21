from __future__ import annotations

from collections import Counter
from pathlib import Path

from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.branch_lengths.ultrametric import (
    APE_ULTRAMETRIC_TOLERANCE,
    summarize_ultrametric_tip_depths,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree

from .models import RootStateConfidenceReport, StableNodeIdentity, TreeIntegrityIssue

ROOT_BRANCH_LENGTH_IMBALANCE_THRESHOLD = 10.0
ROOT_CHILD_BALANCE_RATIO_THRESHOLD = 0.1


def _load_tree(path: Path, *, source_format: str | None = None) -> PhyloTree:
    return load_tree(path, source_format=source_format)


def _count_polytomies(tree: PhyloTree) -> int:
    return sum(
        1 for node in tree.iter_nodes() if not node.is_leaf() and len(node.children) > 2
    )


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

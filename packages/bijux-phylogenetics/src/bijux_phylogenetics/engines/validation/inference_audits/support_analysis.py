from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.engines.artifacts.support import (
    BootstrapSupportNode,
    BootstrapSupportSummaryReport,
    FastTreeSupportNode,
    FastTreeSupportSummaryReport,
    ShAlrtSupportNode,
    ShAlrtSupportSummaryReport,
    WeakBackboneReport,
)
from bijux_phylogenetics.io.fasttree_support import (
    FastTreeBranchSupportLabel,
    parse_fasttree_branch_support_label,
)
from bijux_phylogenetics.io.iqtree_support import (
    IqtreeBranchSupportLabel,
    parse_iqtree_branch_support_label,
    support_fraction,
)
from bijux_phylogenetics.io.trees import load_tree


def summarize_bootstrap_support_distribution(
    tree_path: Path,
    *,
    weak_support_threshold: float = 70.0,
) -> BootstrapSupportSummaryReport:
    """Summarize internal-node support values and their distribution across one tree."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

    tree = load_tree(tree_path)
    nodes: list[BootstrapSupportNode] = []
    warnings: list[str] = []
    total_tip_count = tree.tip_count
    internal_node_count = sum(
        1 for node in tree.iter_nodes() if node is not tree.root and not node.is_leaf()
    )
    for node in tree.iter_nodes():
        if node is tree.root or node.is_leaf():
            continue
        descendant_taxa = node_descendant_taxa(node)
        support_label = _parse_iqtree_support_label(node.name)
        support = None if support_label is None else support_label.ufboot_support
        if support is None and _is_collapsed_consensus_branch_without_label(node):
            support = 0.0
        if support is None:
            continue
        normalized_support_fraction = support_fraction(support)
        if normalized_support_fraction is None:
            continue
        nodes.append(
            BootstrapSupportNode(
                node="|".join(descendant_taxa)
                if descendant_taxa
                else (node.name or "<unnamed>"),
                descendant_taxa=descendant_taxa,
                support=support,
                support_fraction=normalized_support_fraction,
                is_backbone=len(descendant_taxa) >= max(2, total_tip_count // 2),
            )
        )
    histogram = {
        "lt50": sum(1 for node in nodes if node.support < 50.0),
        "50to69": sum(1 for node in nodes if 50.0 <= node.support < 70.0),
        "70to89": sum(1 for node in nodes if 70.0 <= node.support < 90.0),
        "ge90": sum(1 for node in nodes if node.support >= 90.0),
    }
    supports = sorted(node.support for node in nodes)
    if len(nodes) < internal_node_count:
        warnings.append(
            "one or more internal nodes did not expose numeric support labels"
        )
    if any(node.support < weak_support_threshold for node in nodes):
        warnings.append("one or more internal clades remain weakly supported")
    return BootstrapSupportSummaryReport(
        tree_path=tree_path,
        internal_node_count=internal_node_count,
        supported_node_count=len(nodes),
        minimum_support=None if not supports else supports[0],
        maximum_support=None if not supports else supports[-1],
        median_support=_median_support(supports),
        weakly_supported_clade_count=sum(
            1 for node in nodes if node.support < weak_support_threshold
        ),
        support_histogram=histogram,
        nodes=nodes,
        warnings=warnings,
    )


def _is_collapsed_consensus_branch_without_label(node: object) -> bool:
    branch_length = float(getattr(node, "branch_length", 0.0) or 0.0)
    if branch_length > 1e-5:
        return False
    children = tuple(getattr(node, "children", ()))
    if len(children) < 2:
        return False
    return all(
        float(getattr(child, "branch_length", 0.0) or 0.0) <= 1e-5 for child in children
    )


def summarize_fasttree_support_distribution(
    tree_path: Path,
    *,
    weak_support_threshold: float = 0.7,
) -> FastTreeSupportSummaryReport:
    """Summarize FastTree SH-like local support labels across one tree."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

    tree = load_tree(tree_path)
    nodes: list[FastTreeSupportNode] = []
    warnings: list[str] = []
    total_tip_count = tree.tip_count
    internal_node_count = sum(1 for node in tree.iter_nodes() if not node.is_leaf())
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = node_descendant_taxa(node)
        support_label = _parse_fasttree_support_label(node.name)
        if support_label is None:
            continue
        local_support = support_label.local_support
        nodes.append(
            FastTreeSupportNode(
                node="|".join(descendant_taxa)
                if descendant_taxa
                else (node.name or "<unnamed>"),
                descendant_taxa=descendant_taxa,
                local_support=local_support,
                support_fraction=local_support,
                is_backbone=len(descendant_taxa) >= max(2, total_tip_count // 2),
            )
        )
    histogram = {
        "lt0p5": sum(1 for node in nodes if node.local_support < 0.5),
        "0p5to0p69": sum(1 for node in nodes if 0.5 <= node.local_support < 0.7),
        "0p7to0p89": sum(1 for node in nodes if 0.7 <= node.local_support < 0.9),
        "ge0p9": sum(1 for node in nodes if node.local_support >= 0.9),
    }
    supports = sorted(node.local_support for node in nodes)
    if len(nodes) < internal_node_count:
        warnings.append(
            "one or more internal nodes did not expose parsable FastTree local support labels"
        )
    if any(node.local_support < weak_support_threshold for node in nodes):
        warnings.append(
            "one or more internal clades remain weakly supported by FastTree local support"
        )
    warnings.append(
        "FastTree is an approximately maximum-likelihood workflow and its SH-like local support values should be reviewed as approximate evidence"
    )
    return FastTreeSupportSummaryReport(
        tree_path=tree_path,
        internal_node_count=internal_node_count,
        annotated_node_count=len(nodes),
        minimum_local_support=None if not supports else supports[0],
        maximum_local_support=None if not supports else supports[-1],
        median_local_support=_median_support(supports),
        weakly_supported_clade_count=sum(
            1 for node in nodes if node.local_support < weak_support_threshold
        ),
        support_histogram=histogram,
        approximate_method=True,
        support_label_kind="sh-like-local-support",
        support_scale="proportion-0-to-1",
        nodes=nodes,
        warnings=warnings,
    )


def detect_weakly_supported_backbone(
    tree_path: Path,
    *,
    threshold: float = 70.0,
) -> WeakBackboneReport:
    """Flag broad internal clades whose support falls below a declared backbone threshold."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

    summary = summarize_bootstrap_support_distribution(
        tree_path, weak_support_threshold=threshold
    )
    tree = load_tree(tree_path)
    total_tip_count = tree.tip_count
    backbone_nodes = list(summary.nodes)
    observed_backbone_labels = {
        node.node for node in backbone_nodes if node.is_backbone
    }
    root = tree.root
    if not root.is_leaf():
        descendant_taxa = node_descendant_taxa(root)
        support_label = _parse_iqtree_support_label(root.name)
        support = None if support_label is None else support_label.ufboot_support
        normalized_support_fraction = support_fraction(support)
        if normalized_support_fraction is not None:
            root_label = (
                "|".join(descendant_taxa)
                if descendant_taxa
                else (root.name or "<unnamed>")
            )
            if root_label not in observed_backbone_labels:
                backbone_nodes.append(
                    BootstrapSupportNode(
                        node=root_label,
                        descendant_taxa=descendant_taxa,
                        support=support,
                        support_fraction=normalized_support_fraction,
                        is_backbone=len(descendant_taxa)
                        >= max(2, total_tip_count // 2),
                    )
                )
    weak_nodes = [
        node for node in backbone_nodes if node.is_backbone and node.support < threshold
    ]
    warnings = list(summary.warnings)
    if weak_nodes:
        warnings.append(
            "major internal branches remain weakly supported along the backbone"
        )
    return WeakBackboneReport(
        tree_path=tree_path,
        threshold=threshold,
        evaluated_backbone_node_count=sum(
            1 for node in backbone_nodes if node.is_backbone
        ),
        weak_backbone_node_count=len(weak_nodes),
        weak_nodes=weak_nodes,
        warnings=warnings,
    )


def summarize_sh_alrt_support_distribution(
    tree_path: Path,
    *,
    sh_alrt_strong_threshold: float = 80.0,
    ufboot_strong_threshold: float = 95.0,
) -> ShAlrtSupportSummaryReport:
    """Summarize compound SH-aLRT/UFBoot branch labels across one tree."""
    from bijux_phylogenetics.ancestral.common import node_descendant_taxa

    tree = load_tree(tree_path)
    nodes: list[ShAlrtSupportNode] = []
    warnings: list[str] = []
    total_tip_count = tree.tip_count
    for node in tree.iter_nodes():
        if node.is_leaf():
            continue
        descendant_taxa = node_descendant_taxa(node)
        support_label = _parse_iqtree_support_label(node.name)
        if support_label is None:
            continue
        sh_alrt_value = support_label.sh_alrt_support
        ufboot_value = support_label.ufboot_support
        sh_alrt_strong = (
            sh_alrt_value is not None and sh_alrt_value >= sh_alrt_strong_threshold
        )
        ufboot_strong = (
            ufboot_value is not None and ufboot_value >= ufboot_strong_threshold
        )
        support_agreement = _support_agreement(
            sh_alrt_value=sh_alrt_value,
            ufboot_value=ufboot_value,
            sh_alrt_strong=sh_alrt_strong,
            ufboot_strong=ufboot_strong,
        )
        nodes.append(
            ShAlrtSupportNode(
                node="|".join(descendant_taxa)
                if descendant_taxa
                else (node.name or "<unnamed>"),
                descendant_taxa=descendant_taxa,
                sh_alrt_support=sh_alrt_value,
                sh_alrt_support_fraction=support_fraction(sh_alrt_value),
                ufboot_support=ufboot_value,
                ufboot_support_fraction=support_fraction(ufboot_value),
                is_backbone=len(descendant_taxa) >= max(2, total_tip_count // 2),
                sh_alrt_strong=sh_alrt_strong,
                ufboot_strong=ufboot_strong,
                conflicting_support_signal=support_agreement
                in {"sh_alrt_only", "ufboot_only"},
                support_agreement=support_agreement,
            )
        )
    if len(nodes) < sum(1 for node in tree.iter_nodes() if not node.is_leaf()):
        warnings.append(
            "one or more internal nodes did not expose parsable sh-alrt or ufboot labels"
        )
    if any(node.conflicting_support_signal for node in nodes):
        warnings.append(
            "one or more internal clades show conflicting sh-alrt and ufboot support signals"
        )
    if any(
        not node.sh_alrt_strong for node in nodes if node.sh_alrt_support is not None
    ):
        warnings.append(
            "one or more internal clades remain weakly supported by sh-alrt"
        )
    if any(not node.ufboot_strong for node in nodes if node.ufboot_support is not None):
        warnings.append("one or more internal clades remain weakly supported by ufboot")
    sh_alrt_values = sorted(
        node.sh_alrt_support for node in nodes if node.sh_alrt_support is not None
    )
    ufboot_values = sorted(
        node.ufboot_support for node in nodes if node.ufboot_support is not None
    )
    return ShAlrtSupportSummaryReport(
        tree_path=tree_path,
        internal_node_count=sum(1 for node in tree.iter_nodes() if not node.is_leaf()),
        annotated_node_count=len(nodes),
        fully_scored_node_count=sum(
            1
            for node in nodes
            if node.sh_alrt_support is not None and node.ufboot_support is not None
        ),
        minimum_sh_alrt_support=None if not sh_alrt_values else sh_alrt_values[0],
        maximum_sh_alrt_support=None if not sh_alrt_values else sh_alrt_values[-1],
        minimum_ufboot_support=None if not ufboot_values else ufboot_values[0],
        maximum_ufboot_support=None if not ufboot_values else ufboot_values[-1],
        weak_sh_alrt_clade_count=sum(
            1
            for node in nodes
            if node.sh_alrt_support is not None and not node.sh_alrt_strong
        ),
        weak_ufboot_clade_count=sum(
            1
            for node in nodes
            if node.ufboot_support is not None and not node.ufboot_strong
        ),
        conflicting_support_signal_count=sum(
            1 for node in nodes if node.conflicting_support_signal
        ),
        nodes=nodes,
        warnings=warnings,
    )


def _parse_iqtree_support_label(
    value: str | None,
) -> IqtreeBranchSupportLabel | None:
    return parse_iqtree_branch_support_label(value)


def _parse_fasttree_support_label(
    value: str | None,
) -> FastTreeBranchSupportLabel | None:
    return parse_fasttree_branch_support_label(value)


def _support_agreement(
    *,
    sh_alrt_value: float | None,
    ufboot_value: float | None,
    sh_alrt_strong: bool,
    ufboot_strong: bool,
) -> str:
    if sh_alrt_value is None or ufboot_value is None:
        return "incomplete"
    if sh_alrt_strong and ufboot_strong:
        return "both_strong"
    if not sh_alrt_strong and not ufboot_strong:
        return "both_weak"
    if sh_alrt_strong:
        return "sh_alrt_only"
    return "ufboot_only"


def _median_support(values: list[float]) -> float | None:
    if not values:
        return None
    midpoint = len(values) // 2
    if len(values) % 2 == 1:
        return values[midpoint]
    return (values[midpoint - 1] + values[midpoint]) / 2.0

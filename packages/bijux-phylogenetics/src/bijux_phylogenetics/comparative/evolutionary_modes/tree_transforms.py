from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.comparative.common import node_signature
from bijux_phylogenetics.comparative.evolutionary_modes.models import (
    ComparativeTreeRescalingReport,
    EvolutionaryModeBranchLengthRow,
)
from bijux_phylogenetics.comparative.evolutionary_modes.numeric import stable_float
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import ComparativeMethodError


def rescale_tree_ornstein_uhlenbeck(
    tree_path: Path,
    *,
    alpha: float,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style OU branch rescaling to a rooted tree."""
    return build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="ornstein-uhlenbeck",
        parameter_name="alpha",
        parameter_value=alpha,
        sigsq=sigsq,
    )


def rescale_tree_early_burst(
    tree_path: Path,
    *,
    rate_change: float,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style early-burst branch rescaling to a rooted tree."""
    return build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="early-burst",
        parameter_name="rate_change",
        parameter_value=rate_change,
        sigsq=sigsq,
    )


def rescale_tree_pagel_lambda(
    tree_path: Path,
    *,
    lambda_value: float,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style Pagel-lambda rescaling to a rooted tree."""
    return build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="pagel-lambda",
        parameter_name="lambda",
        parameter_value=lambda_value,
        sigsq=1.0,
    )


def rescale_tree_pagel_kappa(
    tree_path: Path,
    *,
    kappa: float,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style Pagel-kappa branch-length rescaling to a rooted tree."""
    return build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="pagel-kappa",
        parameter_name="kappa",
        parameter_value=kappa,
        sigsq=1.0,
    )


def rescale_tree_pagel_delta(
    tree_path: Path,
    *,
    delta: float,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style Pagel-delta depth rescaling to a rooted tree."""
    return build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="pagel-delta",
        parameter_name="delta",
        parameter_value=delta,
        sigsq=1.0,
    )


def rescale_tree_white_noise(
    tree_path: Path,
    *,
    sigsq: float = 1.0,
) -> ComparativeTreeRescalingReport:
    """Apply the geiger-style white no-phylogeny tree rescaling."""
    return build_tree_rescaling_report(
        load_tree(tree_path),
        tree_path,
        mode="white-noise",
        parameter_name="sigsq",
        parameter_value=sigsq,
        sigsq=sigsq,
    )


def transform_tree_for_evolutionary_mode(
    tree: PhyloTree,
    *,
    mode: str,
    parameter_value: float,
    sigsq: float = 1.0,
) -> PhyloTree:
    """Transform an in-memory tree under a governed continuous-mode branch rule."""
    return transform_tree(
        tree,
        mode=mode,
        parameter_value=parameter_value,
        sigsq=sigsq,
    )


def build_tree_rescaling_report(
    tree: PhyloTree,
    tree_path: Path,
    *,
    mode: str,
    parameter_name: str,
    parameter_value: float,
    sigsq: float,
) -> ComparativeTreeRescalingReport:
    transformed_tree = transform_tree(
        tree,
        mode=mode,
        parameter_value=parameter_value,
        sigsq=sigsq,
    )
    branch_rows = branch_length_rows(
        original_tree=tree,
        transformed_tree=transformed_tree,
    )
    return ComparativeTreeRescalingReport(
        tree_path=tree_path,
        mode=mode,
        parameter_name=parameter_name,
        parameter_value=stable_float(parameter_value),
        tip_count=tree.tip_count,
        original_total_branch_length=stable_float(tree.total_branch_length()),
        transformed_total_branch_length=stable_float(
            transformed_tree.total_branch_length()
        ),
        transformed_tree_newick=dumps_newick(transformed_tree),
        branch_rows=branch_rows,
    )


def transform_tree(
    tree: PhyloTree,
    *,
    mode: str,
    parameter_value: float,
    sigsq: float = 1.0,
) -> PhyloTree:
    if mode not in {
        "ornstein-uhlenbeck",
        "early-burst",
        "pagel-lambda",
        "pagel-kappa",
        "pagel-delta",
        "white-noise",
    }:
        raise ComparativeMethodError(
            "tree transformation mode must be 'ornstein-uhlenbeck', 'early-burst', 'pagel-lambda', 'pagel-kappa', 'pagel-delta', or 'white-noise'"
        )
    reject_negative_transform_branch_lengths(tree, mode=mode)
    if mode == "ornstein-uhlenbeck" and parameter_value < 0.0:
        raise ComparativeMethodError("OU alpha must be non-negative")
    if mode == "pagel-lambda" and not 0.0 <= parameter_value <= 1.0:
        raise ComparativeMethodError("Pagel lambda must lie within [0, 1]")
    if mode == "pagel-kappa" and parameter_value < 0.0:
        raise ComparativeMethodError("Pagel kappa must be non-negative")
    if mode == "pagel-delta" and parameter_value < 0.0:
        raise ComparativeMethodError("Pagel delta must be non-negative")
    if mode == "white-noise" and sigsq < 0.0:
        raise ComparativeMethodError("White-noise sigsq must be non-negative")
    cloned_root = clone_node(tree.root)
    if mode == "pagel-lambda":

        def visit_pagel_lambda(node: TreeNode, depth: float) -> None:
            for child in node.children:
                original_length = float(child.branch_length or 0.0)
                if child.is_leaf():
                    child.branch_length = original_length + (
                        (1.0 - parameter_value) * depth
                    )
                else:
                    child.branch_length = original_length * parameter_value
                visit_pagel_lambda(child, depth + original_length)

        visit_pagel_lambda(cloned_root, 0.0)
        return PhyloTree(
            root=cloned_root,
            source_format=tree.source_format,
            rooted=tree.rooted,
        )
    if mode == "pagel-kappa":

        def visit_pagel_kappa(node: TreeNode) -> None:
            for child in node.children:
                original_length = float(child.branch_length or 0.0)
                child.branch_length = kappa_branch_length(
                    original_length,
                    kappa=parameter_value,
                )
                visit_pagel_kappa(child)

        visit_pagel_kappa(cloned_root)
        return PhyloTree(
            root=cloned_root,
            source_format=tree.source_format,
            rooted=tree.rooted,
        )
    if mode == "pagel-delta":
        total_depth = max_tip_depth(tree.root, depth=0.0)

        def visit_pagel_delta(
            node: TreeNode,
            depth: float,
            transformed_depth: float,
        ) -> None:
            for child in node.children:
                original_length = float(child.branch_length or 0.0)
                child_depth = depth + original_length
                transformed_child_depth = delta_transformed_depth(
                    child_depth,
                    total_depth=total_depth,
                    delta=parameter_value,
                )
                child.branch_length = max(
                    0.0,
                    transformed_child_depth - transformed_depth,
                )
                visit_pagel_delta(child, child_depth, transformed_child_depth)

        visit_pagel_delta(cloned_root, 0.0, 0.0)
        return PhyloTree(
            root=cloned_root,
            source_format=tree.source_format,
            rooted=tree.rooted,
        )
    if mode == "white-noise":

        def visit_white(node: TreeNode) -> None:
            for child in node.children:
                child.branch_length = sigsq if child.is_leaf() else 0.0
                visit_white(child)

        visit_white(cloned_root)
        return PhyloTree(
            root=cloned_root,
            source_format=tree.source_format,
            rooted=tree.rooted,
        )

    total_depth = max_tip_depth(tree.root, depth=0.0)

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            original_length = float(child.branch_length or 0.0)
            child_depth = depth + original_length
            if mode == "ornstein-uhlenbeck":
                child.branch_length = ou_branch_length(
                    parent_depth=depth,
                    child_depth=child_depth,
                    total_depth=total_depth,
                    alpha=parameter_value,
                    sigsq=sigsq,
                )
            else:
                child.branch_length = early_burst_branch_length(
                    parent_depth=depth,
                    child_depth=child_depth,
                    rate_change=parameter_value,
                    sigsq=sigsq,
                )
            visit(child, child_depth)

    visit(cloned_root, 0.0)
    return PhyloTree(
        root=cloned_root,
        source_format=tree.source_format,
        rooted=tree.rooted,
    )


def branch_length_rows(
    *,
    original_tree: PhyloTree,
    transformed_tree: PhyloTree,
) -> list[EvolutionaryModeBranchLengthRow]:
    original_rows = tree_branch_lookup(original_tree)
    transformed_rows = tree_branch_lookup(transformed_tree)
    rows: list[EvolutionaryModeBranchLengthRow] = []
    for node_id in sorted(original_rows):
        original = original_rows[node_id]
        transformed = transformed_rows[node_id]
        rows.append(
            EvolutionaryModeBranchLengthRow(
                node=node_id,
                descendant_taxa=list(original["descendant_taxa"]),
                original_branch_length=stable_float(original["branch_length"]),
                transformed_branch_length=stable_float(transformed["branch_length"]),
                parent_depth=stable_float(original["parent_depth"]),
                child_depth=stable_float(original["child_depth"]),
            )
        )
    return rows


def tree_branch_lookup(tree: PhyloTree) -> dict[str, dict[str, object]]:
    rows: dict[str, dict[str, object]] = {}

    def visit(node: TreeNode, depth: float) -> None:
        for child in node.children:
            branch_length = float(child.branch_length or 0.0)
            child_depth = depth + branch_length
            branch_id = node_signature(child)
            rows[branch_id] = {
                "branch_length": branch_length,
                "parent_depth": depth,
                "child_depth": child_depth,
                "descendant_taxa": descendant_taxa(child),
            }
            visit(child, child_depth)

    visit(tree.root, 0.0)
    return rows


def ou_branch_length(
    *,
    parent_depth: float,
    child_depth: float,
    total_depth: float,
    alpha: float,
    sigsq: float,
) -> float:
    if alpha <= 0.0:
        raise ComparativeMethodError("OU alpha must be positive")

    def term(depth: float) -> float:
        return (
            (1.0 / (2.0 * alpha))
            * math.exp(-2.0 * alpha * (total_depth - depth))
            * (1.0 - math.exp(-2.0 * alpha * depth))
        )

    return max(0.0, (term(child_depth) - term(parent_depth)) * sigsq)


def early_burst_branch_length(
    *,
    parent_depth: float,
    child_depth: float,
    rate_change: float,
    sigsq: float,
) -> float:
    if math.isclose(rate_change, 0.0, rel_tol=0.0, abs_tol=1e-12):
        return max(0.0, (child_depth - parent_depth) * sigsq)
    transformed = (
        math.exp(-rate_change * parent_depth) - math.exp(-rate_change * child_depth)
    ) / rate_change
    return max(0.0, transformed * sigsq)


def kappa_branch_length(branch_length: float, *, kappa: float) -> float:
    if kappa < 0.0:
        raise ComparativeMethodError("Pagel kappa must be non-negative")
    if branch_length < 0.0:
        raise ComparativeMethodError(
            "Pagel kappa cannot transform negative branch lengths"
        )
    transformed = math.pow(branch_length, kappa)
    if not math.isfinite(transformed) or transformed < 0.0:
        raise ComparativeMethodError(
            "Pagel kappa produced an invalid transformed branch length"
        )
    return transformed


def delta_transformed_depth(
    depth: float,
    *,
    total_depth: float,
    delta: float,
) -> float:
    if delta < 0.0:
        raise ComparativeMethodError("Pagel delta must be non-negative")
    if depth < 0.0 or total_depth < 0.0:
        raise ComparativeMethodError(
            "Pagel delta cannot transform negative node depths"
        )
    if math.isclose(depth, 0.0, rel_tol=0.0, abs_tol=1e-12):
        return 0.0
    if math.isclose(total_depth, 0.0, rel_tol=0.0, abs_tol=1e-12):
        return 0.0
    proportion = depth / total_depth
    transformed = total_depth * math.pow(proportion, delta)
    if not math.isfinite(transformed) or transformed < 0.0:
        raise ComparativeMethodError(
            "Pagel delta produced an invalid transformed node depth"
        )
    return transformed


def identity_covariance_matrix(size: int) -> list[list[float]]:
    return [
        [1.0 if row_index == column_index else 0.0 for column_index in range(size)]
        for row_index in range(size)
    ]


def reject_negative_transform_branch_lengths(tree: PhyloTree, *, mode: str) -> None:
    message_map = {
        "pagel-lambda": "Pagel lambda cannot transform negative branch lengths",
        "pagel-kappa": "Pagel kappa cannot transform negative branch lengths",
        "pagel-delta": "Pagel delta cannot transform negative branch lengths",
        "early-burst": "Early-burst rescaling cannot transform negative branch lengths",
        "white-noise": "White-noise rescaling cannot transform negative branch lengths",
        "ornstein-uhlenbeck": "OU rescaling cannot transform negative branch lengths",
    }

    def visit(node: TreeNode) -> None:
        for child in node.children:
            branch_length = float(child.branch_length or 0.0)
            if branch_length < 0.0:
                raise ComparativeMethodError(message_map[mode])
            visit(child)

    visit(tree.root)


def clone_node(node: TreeNode) -> TreeNode:
    return TreeNode(
        name=node.name,
        branch_length=node.branch_length,
        children=[clone_node(child) for child in node.children],
    )


def clone_tree(tree: PhyloTree) -> PhyloTree:
    return PhyloTree(
        root=clone_node(tree.root),
        source_format=tree.source_format,
        rooted=tree.rooted,
    )


def max_tip_depth(node: TreeNode, *, depth: float) -> float:
    if node.is_leaf():
        return depth
    return max(
        max_tip_depth(
            child,
            depth=depth + float(child.branch_length or 0.0),
        )
        for child in node.children
    )


def descendant_taxa(node: TreeNode) -> list[str]:
    if node.is_leaf():
        return [node.name] if node.name is not None else []
    taxa: list[str] = []
    for child in node.children:
        taxa.extend(descendant_taxa(child))
    return sorted(taxa)

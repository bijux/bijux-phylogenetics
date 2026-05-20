from __future__ import annotations

import math
from math import sqrt
import random

from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def _iter_tip_trait_values(
    tree: PhyloTree,
    *,
    root_state: float,
    propagate,
) -> dict[str, float]:
    values: dict[str, float] = {}

    def visit(node: TreeNode, state: float) -> None:
        if node.is_leaf():
            if node.name is not None:
                values[node.name] = (
                    round(state, 15) if isinstance(state, float) else state
                )
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            visit(child, propagate(state, branch_length))

    visit(tree.root, root_state)
    return values


def _iter_node_trait_values(
    tree: PhyloTree,
    *,
    root_state,
    propagate,
) -> dict[str, object]:
    from bijux_phylogenetics.ancestral.common import node_signature

    values: dict[str, object] = {}

    def visit(node: TreeNode, state) -> None:
        values[node_signature(node)] = state
        if node.is_leaf():
            return
        for child in node.children:
            branch_length = max(child.branch_length or 0.0, 0.0)
            visit(child, propagate(state, branch_length))

    visit(tree.root, root_state)
    return values


def _tip_values_from_node_map(
    tree: PhyloTree,
    node_values: dict[str, object],
) -> dict[str, object]:
    from bijux_phylogenetics.ancestral.common import node_signature

    return {
        node.name: (
            round(float(node_values[node_signature(node)]), 15)
            if isinstance(node_values[node_signature(node)], float)
            else node_values[node_signature(node)]
        )
        for node in tree.iter_leaves()
        if node.name is not None
    }


def _resolve_brownian_sigma_parameters(
    *,
    sigma: float | None,
    sigma_squared: float | None,
) -> tuple[float, float]:
    if sigma is None and sigma_squared is None:
        return 1.0, 1.0
    if sigma is not None and sigma < 0.0:
        raise ValueError(f"sigma must be nonnegative, got {sigma}")
    if sigma_squared is not None and sigma_squared < 0.0:
        raise ValueError(f"sigma_squared must be nonnegative, got {sigma_squared}")
    if sigma is None:
        resolved_sigma_squared = float(sigma_squared)
        return sqrt(resolved_sigma_squared), resolved_sigma_squared
    resolved_sigma_squared = sigma * sigma
    if sigma_squared is None:
        return sigma, resolved_sigma_squared
    if not math.isclose(
        resolved_sigma_squared,
        sigma_squared,
        rel_tol=1e-12,
        abs_tol=1e-12,
    ):
        raise ValueError(
            "sigma and sigma_squared must describe the same Brownian rate parameter"
        )
    return sigma, float(sigma_squared)


def _simulate_brownian_node_values(
    tree: PhyloTree,
    *,
    root_state: float,
    sigma: float,
    rng: random.Random,
) -> dict[str, float]:
    return _iter_node_trait_values(
        tree,
        root_state=root_state,
        propagate=lambda state, branch_length: (
            state + rng.gauss(0.0, sigma * sqrt(branch_length))
        ),
    )

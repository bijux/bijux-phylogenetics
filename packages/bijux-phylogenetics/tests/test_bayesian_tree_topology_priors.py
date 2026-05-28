from __future__ import annotations

import math

import pytest

from bijux_phylogenetics.bayesian.tree_topology_priors import (
    build_uniform_rooted_tree_topology_prior,
    evaluate_tree_topology_log_prior,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import PhylogeneticsError, UnrootedTreeError


def test_equivalent_child_order_trees_receive_the_same_topology_prior() -> None:
    prior_model = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])
    left_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    children=[
                        TreeNode(name="A"),
                        TreeNode(name="B"),
                    ]
                ),
                TreeNode(
                    children=[
                        TreeNode(name="C"),
                        TreeNode(name="D"),
                    ]
                ),
            ]
        ),
        rooted=True,
    )
    right_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    children=[
                        TreeNode(name="D"),
                        TreeNode(name="C"),
                    ]
                ),
                TreeNode(
                    children=[
                        TreeNode(name="B"),
                        TreeNode(name="A"),
                    ]
                ),
            ]
        ),
        rooted=True,
    )

    left_report = evaluate_tree_topology_log_prior(left_tree, prior_model)
    right_report = evaluate_tree_topology_log_prior(right_tree, prior_model)

    assert left_report.topology_id == right_report.topology_id
    assert math.isclose(
        left_report.log_prior,
        right_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert left_report.topology_count == 15


def test_distinct_labeled_topologies_are_handled_consistently() -> None:
    prior_model = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])
    balanced_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(
                    children=[
                        TreeNode(name="A"),
                        TreeNode(name="B"),
                    ]
                ),
                TreeNode(
                    children=[
                        TreeNode(name="C"),
                        TreeNode(name="D"),
                    ]
                ),
            ]
        ),
        rooted=True,
    )
    pectinate_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A"),
                TreeNode(
                    children=[
                        TreeNode(name="B"),
                        TreeNode(
                            children=[
                                TreeNode(name="C"),
                                TreeNode(name="D"),
                            ]
                        ),
                    ]
                ),
            ]
        ),
        rooted=True,
    )

    balanced_report = evaluate_tree_topology_log_prior(balanced_tree, prior_model)
    pectinate_report = evaluate_tree_topology_log_prior(pectinate_tree, prior_model)

    assert balanced_report.topology_id != pectinate_report.topology_id
    assert math.isclose(
        balanced_report.log_prior,
        pectinate_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_tree_topology_prior_requires_the_exact_prior_taxon_set() -> None:
    prior_model = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])
    mismatched_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A"),
                TreeNode(
                    children=[
                        TreeNode(name="B"),
                        TreeNode(
                            children=[
                                TreeNode(name="C"),
                                TreeNode(name="E"),
                            ]
                        ),
                    ]
                ),
            ]
        ),
        rooted=True,
    )

    with pytest.raises(
        PhylogeneticsError,
        match="tree topology prior requires the exact prior taxon set",
    ):
        evaluate_tree_topology_log_prior(mismatched_tree, prior_model)


def test_tree_topology_prior_requires_a_rooted_tree() -> None:
    prior_model = build_uniform_rooted_tree_topology_prior(["A", "B", "C"])
    unrooted_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A"),
                TreeNode(
                    children=[
                        TreeNode(name="B"),
                        TreeNode(name="C"),
                    ]
                ),
            ]
        ),
        rooted=False,
    )

    with pytest.raises(UnrootedTreeError, match="requires a rooted tree"):
        evaluate_tree_topology_log_prior(unrooted_tree, prior_model)


def test_tree_topology_prior_requires_a_strictly_bifurcating_tree() -> None:
    prior_model = build_uniform_rooted_tree_topology_prior(["A", "B", "C", "D"])
    multifurcating_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A"),
                TreeNode(name="B"),
                TreeNode(
                    children=[
                        TreeNode(name="C"),
                        TreeNode(name="D"),
                    ]
                ),
            ]
        ),
        rooted=True,
    )

    with pytest.raises(
        PhylogeneticsError,
        match="requires a strictly bifurcating tree",
    ):
        evaluate_tree_topology_log_prior(multifurcating_tree, prior_model)

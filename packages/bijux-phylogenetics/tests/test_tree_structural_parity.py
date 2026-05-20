from __future__ import annotations

from bijux_phylogenetics.compare.topology import (
    compare_tree_sets_structurally,
    compare_tree_structurally,
)
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode


def _rooted_example_tree(
    *,
    left_branch: float = 0.2,
    root_label: str = "root",
    left_label: str = "left",
) -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            name=root_label,
            children=[
                TreeNode(
                    name=left_label,
                    branch_length=left_branch,
                    children=[
                        TreeNode(name="A", branch_length=0.1),
                        TreeNode(name="B", branch_length=0.1),
                    ],
                ),
                TreeNode(
                    name="right",
                    branch_length=0.3,
                    children=[
                        TreeNode(name="C", branch_length=0.1),
                        TreeNode(name="D", branch_length=0.1),
                    ],
                ),
            ],
        ),
        rooted=True,
    )


def _rooted_reordered_example_tree() -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            name="root",
            children=[
                TreeNode(
                    name="right",
                    branch_length=0.3,
                    children=[
                        TreeNode(name="D", branch_length=0.1),
                        TreeNode(name="C", branch_length=0.1),
                    ],
                ),
                TreeNode(
                    name="left",
                    branch_length=0.2,
                    children=[
                        TreeNode(name="B", branch_length=0.1),
                        TreeNode(name="A", branch_length=0.1),
                    ],
                ),
            ],
        ),
        rooted=True,
    )


def _unrooted_example_tree() -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.1),
                TreeNode(
                    branch_length=0.2,
                    children=[
                        TreeNode(name="B", branch_length=0.1),
                        TreeNode(
                            branch_length=0.2,
                            children=[
                                TreeNode(name="C", branch_length=0.1),
                                TreeNode(name="D", branch_length=0.1),
                            ],
                        ),
                    ],
                ),
            ],
        ),
        rooted=False,
    )


def _unrooted_reordered_example_tree() -> PhyloTree:
    return PhyloTree(
        root=TreeNode(
            children=[
                TreeNode(
                    branch_length=0.2,
                    children=[
                        TreeNode(
                            branch_length=0.2,
                            children=[
                                TreeNode(name="D", branch_length=0.1),
                                TreeNode(name="C", branch_length=0.1),
                            ],
                        ),
                        TreeNode(name="B", branch_length=0.1),
                    ],
                ),
                TreeNode(name="A", branch_length=0.1),
            ],
        ),
        rooted=False,
    )


def test_compare_tree_structurally_accepts_reordered_rooted_children() -> None:
    report = compare_tree_structurally(
        _rooted_example_tree(),
        _rooted_reordered_example_tree(),
        tolerance=0.0,
    )

    assert report.equivalent is True
    assert report.mismatch_reason is None


def test_compare_tree_structurally_accepts_reordered_unrooted_children() -> None:
    report = compare_tree_structurally(
        _unrooted_example_tree(),
        _unrooted_reordered_example_tree(),
        tolerance=0.0,
    )

    assert report.equivalent is True
    assert report.mismatch_reason is None


def test_compare_tree_structurally_reports_branch_length_mismatch() -> None:
    report = compare_tree_structurally(
        _rooted_example_tree(),
        _rooted_example_tree(left_branch=0.25),
        tolerance=0.0,
    )

    assert report.equivalent is False
    assert report.mismatch_reason is not None
    assert "branch lengths differ" in report.mismatch_reason


def test_compare_tree_structurally_reports_topology_mismatch() -> None:
    observed = PhyloTree(
        root=TreeNode(
            name="root",
            children=[
                TreeNode(
                    name="left",
                    branch_length=0.2,
                    children=[
                        TreeNode(name="A", branch_length=0.1),
                        TreeNode(name="C", branch_length=0.1),
                    ],
                ),
                TreeNode(
                    name="right",
                    branch_length=0.3,
                    children=[
                        TreeNode(name="B", branch_length=0.1),
                        TreeNode(name="D", branch_length=0.1),
                    ],
                ),
            ],
        ),
        rooted=True,
    )

    report = compare_tree_structurally(
        _rooted_example_tree(),
        observed,
        tolerance=0.0,
    )

    assert report.equivalent is False
    assert report.mismatch_reason == "clades or splits differ between trees"


def test_compare_tree_structurally_reports_rootedness_mismatch() -> None:
    report = compare_tree_structurally(
        _rooted_example_tree(),
        _unrooted_example_tree(),
        tolerance=0.0,
    )

    assert report.equivalent is False
    assert (
        report.mismatch_reason
        == "tree rootedness differs: expected True, observed False"
    )


def test_compare_tree_structurally_reports_internal_label_mismatch() -> None:
    report = compare_tree_structurally(
        _rooted_example_tree(),
        _rooted_example_tree(root_label="root", left_label="changed"),
        tolerance=0.0,
    )

    assert report.equivalent is False
    assert report.mismatch_reason is not None
    assert "internal labels differ" in report.mismatch_reason


def test_compare_tree_sets_structurally_reports_tree_specific_mismatch() -> None:
    report = compare_tree_sets_structurally(
        [_rooted_example_tree(), _rooted_reordered_example_tree()],
        [_rooted_example_tree(), _rooted_example_tree(left_branch=0.25)],
        tolerance=0.0,
    )

    assert report.equivalent is False
    assert report.mismatch_reason is not None
    assert report.mismatch_reason.startswith("tree 2:")

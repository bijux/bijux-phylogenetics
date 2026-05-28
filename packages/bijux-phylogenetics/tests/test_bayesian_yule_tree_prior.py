from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_yule_tree_prior,
    evaluate_yule_tree_log_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    NonUltrametricTreeError,
    PhylogeneticsError,
    UnrootedTreeError,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_crown_conditioned_yule_tree_prior_matches_small_ultrametric_fixture() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    prior_model = build_crown_conditioned_yule_tree_prior(speciation_rate=0.5)

    report = evaluate_yule_tree_log_prior(tree, prior_model)

    assert report.family == "crown-conditioned-yule"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.post_root_speciation_count == 2
    assert math.isclose(report.root_age, 3.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.total_branch_length, 9.0, rel_tol=0.0, abs_tol=1e-12)
    assert [row.lineage_count for row in report.interval_rows] == [2, 3, 4]
    assert [row.event_count for row in report.interval_rows] == [1, 1, 0]
    assert math.isclose(
        report.log_prior,
        (2.0 * math.log(0.5)) - (0.5 * 9.0),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_crown_conditioned_yule_tree_prior_handles_two_tip_ultrametric_tree() -> None:
    tree = load_tree(fixture("trees", "example_tree_two_tip_ultrametric.nwk"))
    prior_model = build_crown_conditioned_yule_tree_prior(speciation_rate=2.0)

    report = evaluate_yule_tree_log_prior(tree, prior_model)

    assert report.post_root_speciation_count == 0
    assert len(report.interval_rows) == 1
    assert report.interval_rows[0].event_count == 0
    assert math.isclose(report.total_branch_length, 0.6, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.log_prior, -1.2, rel_tol=0.0, abs_tol=1e-12)


def test_yule_tree_prior_rejects_non_ultrametric_tree() -> None:
    tree = load_tree(fixture("trees", "strict_clock_nonclock_tree_4_taxa.nwk"))
    prior_model = build_crown_conditioned_yule_tree_prior(speciation_rate=0.5)

    with pytest.raises(
        NonUltrametricTreeError,
        match="requires an ultrametric tree",
    ):
        evaluate_yule_tree_log_prior(tree, prior_model)


def test_yule_tree_prior_rejects_unrooted_tree() -> None:
    prior_model = build_crown_conditioned_yule_tree_prior(speciation_rate=1.0)
    tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A", branch_length=1.0),
                TreeNode(name="B", branch_length=1.0),
                TreeNode(name="C", branch_length=1.0),
            ]
        ),
        rooted=False,
    )

    with pytest.raises(UnrootedTreeError, match="requires a rooted tree"):
        evaluate_yule_tree_log_prior(tree, prior_model)


@pytest.mark.parametrize(
    ("tree", "exception_type", "message"),
    [
        (
            PhyloTree(
                TreeNode(
                    children=[
                        TreeNode(name="A", branch_length=1.0),
                        TreeNode(name="B"),
                    ]
                ),
                rooted=True,
            ),
            InvalidBranchLengthError,
            "requires complete branch lengths",
        ),
        (
            PhyloTree(
                TreeNode(
                    children=[
                        TreeNode(name="A", branch_length=-0.1),
                        TreeNode(name="B", branch_length=0.1),
                    ]
                ),
                rooted=True,
            ),
            InvalidBranchLengthError,
            "requires non-negative branch lengths",
        ),
        (
            PhyloTree(
                TreeNode(
                    children=[
                        TreeNode(name="A", branch_length=1.0),
                        TreeNode(name="B", branch_length=1.0),
                        TreeNode(name="C", branch_length=1.0),
                    ]
                ),
                rooted=True,
            ),
            PhylogeneticsError,
            "requires a strictly bifurcating tree",
        ),
    ],
)
def test_yule_tree_prior_rejects_invalid_tree_structures(
    tree: PhyloTree,
    exception_type: type[Exception],
    message: str,
) -> None:
    prior_model = build_crown_conditioned_yule_tree_prior(speciation_rate=1.0)

    with pytest.raises(exception_type, match=message):
        evaluate_yule_tree_log_prior(tree, prior_model)


def test_yule_tree_prior_rejects_invalid_speciation_rate() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires a strictly positive finite speciation rate",
    ):
        build_crown_conditioned_yule_tree_prior(speciation_rate=0.0)

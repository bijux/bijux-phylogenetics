from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.time_tree_priors import (
    SkylineCoalescentEpoch,
    build_skyline_coalescent_tree_prior,
    evaluate_skyline_coalescent_tree_log_prior,
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


def _expected_skyline_coalescent_log_prior() -> float:
    first_interval = math.log(3.0) - (3.0 * 1.0)
    second_interval = -(1.5 * 0.5) - (0.75 * 0.5) + math.log(0.75)
    third_interval = math.log(0.25) - (0.25 * 1.0)
    return first_interval + second_interval + third_interval


def test_skyline_coalescent_prior_matches_segmented_hand_computation() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    prior_model = build_skyline_coalescent_tree_prior(
        [
            SkylineCoalescentEpoch(
                younger_boundary_age=0.0,
                older_boundary_age=1.5,
                effective_population_size=2.0,
            ),
            SkylineCoalescentEpoch(
                younger_boundary_age=1.5,
                older_boundary_age=None,
                effective_population_size=4.0,
            ),
        ]
    )

    report = evaluate_skyline_coalescent_tree_log_prior(tree, prior_model)

    assert report.family == "skyline-coalescent"
    assert report.epoch_count == 2
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert math.isclose(report.root_age, 3.0, rel_tol=0.0, abs_tol=1e-12)
    assert len(report.segment_rows) == 4
    assert [row.coalescent_interval_index for row in report.segment_rows] == [
        1,
        2,
        2,
        3,
    ]
    assert [row.skyline_epoch_index for row in report.segment_rows] == [1, 1, 2, 2]
    assert [row.duration for row in report.segment_rows] == [1.0, 0.5, 0.5, 1.0]
    assert [row.lineage_count for row in report.segment_rows] == [4, 3, 3, 2]
    assert [row.coalescent_event_count for row in report.segment_rows] == [1, 0, 1, 1]
    assert math.isclose(
        report.log_prior,
        _expected_skyline_coalescent_log_prior(),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_skyline_coalescent_prior_changes_when_epoch_population_changes() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    smaller_second_epoch_prior = build_skyline_coalescent_tree_prior(
        [
            SkylineCoalescentEpoch(0.0, 1.5, 2.0),
            SkylineCoalescentEpoch(1.5, None, 4.0),
        ]
    )
    larger_second_epoch_prior = build_skyline_coalescent_tree_prior(
        [
            SkylineCoalescentEpoch(0.0, 1.5, 2.0),
            SkylineCoalescentEpoch(1.5, None, 8.0),
        ]
    )

    smaller_second_epoch_report = evaluate_skyline_coalescent_tree_log_prior(
        tree,
        smaller_second_epoch_prior,
    )
    larger_second_epoch_report = evaluate_skyline_coalescent_tree_log_prior(
        tree,
        larger_second_epoch_prior,
    )

    assert not math.isclose(
        smaller_second_epoch_report.log_prior,
        larger_second_epoch_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_skyline_coalescent_prior_rejects_non_ultrametric_tree() -> None:
    tree = load_tree(fixture("trees", "strict_clock_nonclock_tree_4_taxa.nwk"))
    prior_model = build_skyline_coalescent_tree_prior(
        [
            SkylineCoalescentEpoch(0.0, 1.5, 2.0),
            SkylineCoalescentEpoch(1.5, None, 4.0),
        ]
    )

    with pytest.raises(
        NonUltrametricTreeError,
        match="requires an ultrametric tree",
    ):
        evaluate_skyline_coalescent_tree_log_prior(tree, prior_model)


def test_skyline_coalescent_prior_rejects_unrooted_tree() -> None:
    prior_model = build_skyline_coalescent_tree_prior(
        [SkylineCoalescentEpoch(0.0, None, 2.0)]
    )
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
        evaluate_skyline_coalescent_tree_log_prior(tree, prior_model)


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
def test_skyline_coalescent_prior_rejects_invalid_tree_structures(
    tree: PhyloTree,
    exception_type: type[Exception],
    message: str,
) -> None:
    prior_model = build_skyline_coalescent_tree_prior(
        [SkylineCoalescentEpoch(0.0, None, 2.0)]
    )

    with pytest.raises(exception_type, match=message):
        evaluate_skyline_coalescent_tree_log_prior(tree, prior_model)


@pytest.mark.parametrize(
    ("epochs", "message"),
    [
        (
            [],
            "requires at least one epoch",
        ),
        (
            [SkylineCoalescentEpoch(0.5, None, 2.0)],
            "requires contiguous epochs from age zero",
        ),
        (
            [
                SkylineCoalescentEpoch(0.0, 1.0, 2.0),
                SkylineCoalescentEpoch(1.5, None, 2.0),
            ],
            "requires contiguous epochs from age zero",
        ),
        (
            [
                SkylineCoalescentEpoch(0.0, None, 2.0),
                SkylineCoalescentEpoch(1.0, None, 2.0),
            ],
            "only the final epoch to be open ended",
        ),
        (
            [SkylineCoalescentEpoch(0.0, 1.0, -1.0)],
            "requires strictly positive finite epoch effective population sizes",
        ),
        (
            [SkylineCoalescentEpoch(0.0, 0.0, 2.0)],
            "requires each finite epoch to end after it starts",
        ),
        (
            [SkylineCoalescentEpoch(0.0, 1.0, 2.0)],
            "requires the final epoch to be open ended",
        ),
    ],
)
def test_skyline_coalescent_prior_rejects_invalid_epochs(
    epochs: list[SkylineCoalescentEpoch],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        build_skyline_coalescent_tree_prior(epochs)

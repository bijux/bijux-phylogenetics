from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_constant_population_coalescent_tree_prior,
    evaluate_constant_population_coalescent_tree_log_prior,
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


def _expected_constant_population_coalescent_log_prior(
    *,
    interval_durations: list[float],
    lineage_counts: list[int],
    effective_population_size: float,
) -> float:
    total = 0.0
    for duration, lineage_count in zip(interval_durations, lineage_counts, strict=True):
        coalescent_rate = math.comb(lineage_count, 2) / effective_population_size
        total += math.log(coalescent_rate) - (coalescent_rate * duration)
    return total


def test_constant_population_coalescent_prior_matches_hand_computed_fixture() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    prior_model = build_constant_population_coalescent_tree_prior(
        effective_population_size=2.0,
    )

    report = evaluate_constant_population_coalescent_tree_log_prior(tree, prior_model)

    expected_log_prior = _expected_constant_population_coalescent_log_prior(
        interval_durations=[1.0, 1.0, 1.0],
        lineage_counts=[4, 3, 2],
        effective_population_size=2.0,
    )
    assert report.family == "constant-population-coalescent"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert math.isclose(report.root_age, 3.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(report.total_branch_length, 9.0, rel_tol=0.0, abs_tol=1e-12)
    assert [row.lineage_count for row in report.interval_rows] == [4, 3, 2]
    assert [row.coalescent_event_count for row in report.interval_rows] == [1, 1, 1]
    assert [row.duration for row in report.interval_rows] == [1.0, 1.0, 1.0]
    assert math.isclose(
        report.log_prior,
        expected_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_constant_population_coalescent_prior_handles_two_tip_ultrametric_tree() -> (
    None
):
    tree = load_tree(fixture("trees", "example_tree_two_tip_ultrametric.nwk"))
    prior_model = build_constant_population_coalescent_tree_prior(
        effective_population_size=2.0,
    )

    report = evaluate_constant_population_coalescent_tree_log_prior(tree, prior_model)

    assert len(report.interval_rows) == 1
    assert report.interval_rows[0].lineage_count == 2
    assert report.interval_rows[0].coalescent_event_count == 1
    assert math.isclose(report.interval_rows[0].duration, 0.3, abs_tol=1e-12)
    assert math.isclose(
        report.log_prior,
        math.log(0.5) - (0.5 * 0.3),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_constant_population_coalescent_prior_changes_with_population_size() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    smaller_population_prior = build_constant_population_coalescent_tree_prior(
        effective_population_size=1.0,
    )
    larger_population_prior = build_constant_population_coalescent_tree_prior(
        effective_population_size=4.0,
    )

    smaller_population_report = evaluate_constant_population_coalescent_tree_log_prior(
        tree,
        smaller_population_prior,
    )
    larger_population_report = evaluate_constant_population_coalescent_tree_log_prior(
        tree,
        larger_population_prior,
    )

    assert not math.isclose(
        smaller_population_report.log_prior,
        larger_population_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_constant_population_coalescent_prior_rejects_non_ultrametric_tree() -> None:
    tree = load_tree(fixture("trees", "strict_clock_nonclock_tree_4_taxa.nwk"))
    prior_model = build_constant_population_coalescent_tree_prior(
        effective_population_size=2.0,
    )

    with pytest.raises(
        NonUltrametricTreeError,
        match="requires an ultrametric tree",
    ):
        evaluate_constant_population_coalescent_tree_log_prior(tree, prior_model)


def test_constant_population_coalescent_prior_rejects_unrooted_tree() -> None:
    prior_model = build_constant_population_coalescent_tree_prior(
        effective_population_size=2.0,
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
        evaluate_constant_population_coalescent_tree_log_prior(tree, prior_model)


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
def test_constant_population_coalescent_prior_rejects_invalid_tree_structures(
    tree: PhyloTree,
    exception_type: type[Exception],
    message: str,
) -> None:
    prior_model = build_constant_population_coalescent_tree_prior(
        effective_population_size=2.0,
    )

    with pytest.raises(exception_type, match=message):
        evaluate_constant_population_coalescent_tree_log_prior(tree, prior_model)


def test_constant_population_coalescent_prior_rejects_invalid_population_size() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="requires a strictly positive finite effective population size",
    ):
        build_constant_population_coalescent_tree_prior(
            effective_population_size=0.0,
        )

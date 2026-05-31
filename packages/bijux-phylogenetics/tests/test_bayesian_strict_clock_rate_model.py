from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.clock_models import (
    build_strict_clock_rate_model,
    evaluate_strict_clock_tree_log_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def load_rooted_tree_fixture(name: str) -> PhyloTree:
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree


def scale_tree_branch_lengths(tree: PhyloTree, *, clock_rate: float) -> PhyloTree:
    scaled_tree = tree.copy()
    for _parent, child in scaled_tree.iter_edges():
        child.branch_length = float(child.branch_length or 0.0) * clock_rate
    return scaled_tree


def test_strict_clock_rate_model_matches_one_global_rate_scaled_tree() -> None:
    dated_tree = load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk")
    substitution_tree = scale_tree_branch_lengths(dated_tree, clock_rate=0.5)
    rate_model = build_strict_clock_rate_model(global_clock_rate=0.5)

    report = evaluate_strict_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        rate_model,
    )

    assert report.family == "strict-clock"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.branch_count == 6
    assert report.exact_match_count == 6
    assert report.mismatch_count == 0
    assert report.total_log_prior == 0.0
    assert report.substitution_tree_newick == report.expected_substitution_tree_newick
    row_by_taxa = {tuple(row.descendant_taxa): row for row in report.branch_rows}
    tip_a_row = row_by_taxa[("A",)]
    assert math.isclose(tip_a_row.dated_time_duration, 1.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        tip_a_row.expected_substitution_branch_length,
        0.5,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        tip_a_row.observed_substitution_branch_length,
        0.5,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        tip_a_row.branch_length_residual,
        0.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert tip_a_row.exact_match is True
    assert tip_a_row.log_prior_contribution == 0.0


def test_strict_clock_rate_model_changes_with_changed_global_rate() -> None:
    dated_tree = load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk")
    substitution_tree = scale_tree_branch_lengths(dated_tree, clock_rate=0.5)

    report = evaluate_strict_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_strict_clock_rate_model(global_clock_rate=0.75),
    )

    assert report.total_log_prior == -math.inf
    assert report.exact_match_count == 0
    assert report.mismatch_count == report.branch_count
    row_by_taxa = {tuple(row.descendant_taxa): row for row in report.branch_rows}
    tip_d_row = row_by_taxa[("D",)]
    assert math.isclose(
        tip_d_row.expected_substitution_branch_length,
        2.25,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        tip_d_row.observed_substitution_branch_length,
        1.5,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        tip_d_row.branch_length_residual,
        -0.75,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert tip_d_row.exact_match is False
    assert tip_d_row.log_prior_contribution == -math.inf


def test_strict_clock_rate_model_rejects_mismatched_rooted_topology() -> None:
    dated_tree = load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk")
    substitution_tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.5),
                TreeNode(
                    branch_length=0.5,
                    children=[
                        TreeNode(name="B", branch_length=0.5),
                        TreeNode(
                            branch_length=0.5,
                            children=[
                                TreeNode(name="C", branch_length=0.5),
                                TreeNode(name="D", branch_length=0.5),
                            ],
                        ),
                    ],
                ),
            ]
        ),
        rooted=True,
    )

    with pytest.raises(
        PhylogeneticsError,
        match="requires identical rooted topology",
    ):
        evaluate_strict_clock_tree_log_prior(
            substitution_tree,
            dated_tree,
            build_strict_clock_rate_model(global_clock_rate=0.5),
        )


@pytest.mark.parametrize(
    ("builder_kwargs", "message"),
    [
        (
            {"global_clock_rate": 0.0},
            "requires a strictly positive finite global clock rate",
        ),
        (
            {"global_clock_rate": 0.5, "branch_length_tolerance": -1e-6},
            "requires a non-negative finite branch-length tolerance",
        ),
    ],
)
def test_strict_clock_rate_model_rejects_invalid_parameters(
    builder_kwargs: dict[str, float],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        build_strict_clock_rate_model(**builder_kwargs)


def test_strict_clock_rate_model_rejects_unrooted_and_incomplete_trees() -> None:
    dated_tree = load_rooted_tree_fixture("strict_clock_time_tree_4_taxa.nwk")
    rooted_substitution_tree = scale_tree_branch_lengths(dated_tree, clock_rate=0.5)
    unrooted_substitution_tree = rooted_substitution_tree.copy()
    unrooted_substitution_tree.rooted = False
    incomplete_dated_tree = dated_tree.copy()
    incomplete_dated_tree.root.children[0].branch_length = None

    with pytest.raises(
        UnrootedTreeError, match="requires one rooted substitution tree"
    ):
        evaluate_strict_clock_tree_log_prior(
            unrooted_substitution_tree,
            dated_tree,
            build_strict_clock_rate_model(global_clock_rate=0.5),
        )

    with pytest.raises(
        InvalidBranchLengthError,
        match="requires complete dated branch durations",
    ):
        evaluate_strict_clock_tree_log_prior(
            rooted_substitution_tree,
            incomplete_dated_tree,
            build_strict_clock_rate_model(global_clock_rate=0.5),
        )

from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.branch_length_priors import (
    build_exponential_branch_length_prior,
    build_fixed_branch_length_prior,
    build_gamma_branch_length_prior,
    build_lognormal_branch_length_prior,
    evaluate_tree_branch_length_log_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology.tree import PhyloTree, TreeNode
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_tree_branch_length_prior_report_sums_branch_rows() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    prior_model = build_exponential_branch_length_prior(rate=2.0)

    report = evaluate_tree_branch_length_log_prior(tree, prior_model)

    assert report.family == "exponential"
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.branch_count == 2
    assert report.parameter_values == {"rate": 2.0}
    assert [row.branch_length for row in report.branch_rows] == [0.1, 0.2]
    assert math.isclose(
        report.total_log_prior,
        sum(row.log_density for row in report.branch_rows),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_tree_branch_length_prior_requires_explicit_branch_lengths() -> None:
    tree = PhyloTree(
        TreeNode(
            children=[
                TreeNode(name="A", branch_length=0.1),
                TreeNode(name="B"),
            ]
        )
    )
    prior_model = build_exponential_branch_length_prior(rate=2.0)

    with pytest.raises(
        InvalidBranchLengthError,
        match="requires explicit branch lengths",
    ):
        evaluate_tree_branch_length_log_prior(tree, prior_model)


@pytest.mark.parametrize(
    ("builder", "kwargs", "message"),
    [
        (
            build_exponential_branch_length_prior,
            {"rate": 0.0},
            "requires 'rate' to be positive",
        ),
        (
            build_gamma_branch_length_prior,
            {"shape": -1.0, "scale": 0.5},
            "requires 'shape' to be positive",
        ),
        (
            build_gamma_branch_length_prior,
            {"shape": 2.0, "scale": 0.0},
            "requires 'scale' to be positive",
        ),
        (
            build_lognormal_branch_length_prior,
            {"log_mean": 0.0, "log_standard_deviation": 0.0},
            "requires 'log_standard_deviation' to be positive",
        ),
        (
            build_fixed_branch_length_prior,
            {"fixed_value": -0.1},
            "requires 'fixed_value' to be non-negative",
        ),
    ],
)
def test_branch_length_prior_builders_reject_invalid_parameters(
    builder,
    kwargs: dict[str, float],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        builder(**kwargs)

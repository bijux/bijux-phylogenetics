from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.time_tree_priors import (
    build_crown_conditioned_birth_death_tree_prior,
    build_crown_conditioned_yule_tree_prior,
    evaluate_birth_death_tree_log_prior,
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


def _expected_crown_conditioned_birth_death_log_prior(
    *,
    branching_ages: list[float],
    crown_age: float,
    tip_count: int,
    speciation_rate: float,
    extinction_rate: float,
    sampling_fraction: float,
) -> float:
    def denominator(time_before_present: float) -> float:
        return speciation_rate * sampling_fraction + (
            (speciation_rate * (1.0 - sampling_fraction)) - extinction_rate
        ) * math.exp(-(speciation_rate - extinction_rate) * time_before_present)

    def p1(time_before_present: float) -> float:
        if math.isclose(speciation_rate, extinction_rate, abs_tol=1e-15):
            return (
                sampling_fraction
                / (1.0 + (sampling_fraction * speciation_rate * time_before_present))
                ** 2
            )
        rate_gap = speciation_rate - extinction_rate
        return (
            sampling_fraction
            * (rate_gap**2)
            * math.exp(-(rate_gap * time_before_present))
            / (denominator(time_before_present) ** 2)
        )

    def q(time_before_present: float) -> float:
        if math.isclose(speciation_rate, extinction_rate, abs_tol=1e-15):
            return (sampling_fraction * time_before_present) / (
                1.0 + (sampling_fraction * speciation_rate * time_before_present)
            )
        rate_gap = speciation_rate - extinction_rate
        return (
            sampling_fraction * (1.0 - math.exp(-(rate_gap * time_before_present)))
        ) / denominator(time_before_present)

    crown_normalization = q(crown_age)
    return (-math.log(tip_count - 1)) + sum(
        math.log(p1(branching_age)) - math.log(crown_normalization)
        for branching_age in branching_ages
    )


def test_crown_conditioned_birth_death_tree_prior_matches_small_ultrametric_fixture() -> (
    None
):
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    prior_model = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=1.2,
        extinction_rate=0.4,
        sampling_fraction=0.75,
    )

    report = evaluate_birth_death_tree_log_prior(tree, prior_model)

    expected_log_prior = _expected_crown_conditioned_birth_death_log_prior(
        branching_ages=[2.0, 1.0],
        crown_age=3.0,
        tip_count=4,
        speciation_rate=1.2,
        extinction_rate=0.4,
        sampling_fraction=0.75,
    )
    assert report.family == "crown-conditioned-birth-death"
    assert report.conditioning_mode == "fixed-tip-count-and-crown-age"
    assert report.tip_count == 4
    assert report.internal_node_count == 3
    assert report.post_root_speciation_count == 2
    assert math.isclose(report.root_age, 3.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        report.transformed_birth_rate,
        0.9,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.transformed_death_rate,
        0.1,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [row.branching_age for row in report.branching_rows] == [2.0, 1.0]
    assert all(
        math.isclose(
            row.log_contribution,
            row.log_branching_probability - row.log_crown_normalization,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in report.branching_rows
    )
    assert math.isclose(
        report.log_prior,
        expected_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_birth_death_tree_prior_changes_with_extinction_and_sampling_fraction() -> None:
    tree = load_tree(fixture("trees", "strict_clock_time_tree_4_taxa.nwk"))
    base_prior = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=0.8,
        extinction_rate=0.1,
        sampling_fraction=1.0,
    )
    higher_extinction_prior = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=0.8,
        extinction_rate=0.3,
        sampling_fraction=1.0,
    )
    incomplete_sampling_prior = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=0.8,
        extinction_rate=0.1,
        sampling_fraction=0.6,
    )
    yule_prior = build_crown_conditioned_yule_tree_prior(speciation_rate=0.8)

    base_report = evaluate_birth_death_tree_log_prior(tree, base_prior)
    higher_extinction_report = evaluate_birth_death_tree_log_prior(
        tree,
        higher_extinction_prior,
    )
    incomplete_sampling_report = evaluate_birth_death_tree_log_prior(
        tree,
        incomplete_sampling_prior,
    )
    yule_report = evaluate_yule_tree_log_prior(tree, yule_prior)

    assert not math.isclose(
        base_report.log_prior,
        higher_extinction_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        base_report.log_prior,
        incomplete_sampling_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        incomplete_sampling_report.log_prior,
        yule_report.log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_birth_death_tree_prior_rejects_non_ultrametric_tree() -> None:
    tree = load_tree(fixture("trees", "strict_clock_nonclock_tree_4_taxa.nwk"))
    prior_model = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=0.8,
        extinction_rate=0.2,
        sampling_fraction=0.7,
    )

    with pytest.raises(
        NonUltrametricTreeError,
        match="requires an ultrametric tree",
    ):
        evaluate_birth_death_tree_log_prior(tree, prior_model)


def test_birth_death_tree_prior_rejects_unrooted_tree() -> None:
    prior_model = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=1.0,
        extinction_rate=0.3,
        sampling_fraction=0.9,
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
        evaluate_birth_death_tree_log_prior(tree, prior_model)


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
def test_birth_death_tree_prior_rejects_invalid_tree_structures(
    tree: PhyloTree,
    exception_type: type[Exception],
    message: str,
) -> None:
    prior_model = build_crown_conditioned_birth_death_tree_prior(
        speciation_rate=1.0,
        extinction_rate=0.3,
        sampling_fraction=0.9,
    )

    with pytest.raises(exception_type, match=message):
        evaluate_birth_death_tree_log_prior(tree, prior_model)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        (
            {
                "speciation_rate": 0.0,
                "extinction_rate": 0.1,
                "sampling_fraction": 1.0,
            },
            "requires a strictly positive finite speciation rate",
        ),
        (
            {
                "speciation_rate": 1.0,
                "extinction_rate": -0.1,
                "sampling_fraction": 1.0,
            },
            "requires a non-negative finite extinction rate",
        ),
        (
            {
                "speciation_rate": 1.0,
                "extinction_rate": 0.1,
                "sampling_fraction": 0.0,
            },
            "requires sampling_fraction in \\(0, 1\\]",
        ),
    ],
)
def test_birth_death_tree_prior_rejects_invalid_parameters(
    kwargs: dict[str, float],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        build_crown_conditioned_birth_death_tree_prior(**kwargs)

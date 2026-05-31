from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.bayesian.clock_models import (
    build_relaxed_lognormal_clock_model,
    evaluate_relaxed_lognormal_clock_tree_log_prior,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import (
    InvalidBranchLengthError,
    PhylogeneticsError,
    UnrootedTreeError,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def load_rooted_tree_fixture(name: str):
    tree = load_tree(fixture("trees", name))
    tree.rooted = True
    return tree


def test_relaxed_lognormal_clock_prior_changes_with_rate_variance() -> None:
    substitution_tree = load_rooted_tree_fixture(
        "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
    )
    dated_tree = load_rooted_tree_fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")

    low_variance_report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_relaxed_lognormal_clock_model(
            rate_policy="independent",
            mean_clock_rate=0.2,
            log_standard_deviation=0.25,
        ),
    )
    high_variance_report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_relaxed_lognormal_clock_model(
            rate_policy="independent",
            mean_clock_rate=0.2,
            log_standard_deviation=0.75,
        ),
    )

    assert low_variance_report.family == "relaxed-lognormal"
    assert low_variance_report.rate_policy == "independent"
    assert low_variance_report.branch_count == 6
    assert not math.isclose(
        low_variance_report.total_log_prior,
        high_variance_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        low_variance_report.minimum_branch_rate,
        0.1,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        low_variance_report.maximum_branch_rate,
        0.6,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_relaxed_lognormal_clock_prior_supports_autocorrelated_policy() -> None:
    substitution_tree = load_rooted_tree_fixture(
        "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
    )
    dated_tree = load_rooted_tree_fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")

    independent_report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_relaxed_lognormal_clock_model(
            rate_policy="independent",
            mean_clock_rate=0.1,
            log_standard_deviation=0.5,
        ),
    )
    autocorrelated_report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_relaxed_lognormal_clock_model(
            rate_policy="autocorrelated",
            mean_clock_rate=0.1,
            log_standard_deviation=0.5,
        ),
    )

    independent_row_by_taxa = {
        tuple(row.descendant_taxa): row for row in independent_report.branch_rows
    }
    autocorrelated_row_by_taxa = {
        tuple(row.descendant_taxa): row for row in autocorrelated_report.branch_rows
    }
    a_independent = independent_row_by_taxa[("A",)]
    a_autocorrelated = autocorrelated_row_by_taxa[("A",)]

    assert independent_report.rate_policy == "independent"
    assert autocorrelated_report.rate_policy == "autocorrelated"
    assert a_independent.anchor_branch_id is None
    assert math.isclose(a_independent.anchor_rate, 0.1, rel_tol=0.0, abs_tol=1e-12)
    assert a_autocorrelated.anchor_branch_id is not None
    assert math.isclose(a_autocorrelated.anchor_rate, 0.2, rel_tol=0.0, abs_tol=1e-12)
    assert not math.isclose(
        independent_report.total_log_prior,
        autocorrelated_report.total_log_prior,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    ab_row = autocorrelated_row_by_taxa[("A", "B")]
    assert math.isclose(ab_row.branch_rate, 0.2, rel_tol=0.0, abs_tol=1e-12)
    assert a_autocorrelated.anchor_branch_id == ab_row.branch_id


def test_relaxed_lognormal_clock_prior_reports_expected_branch_fields() -> None:
    substitution_tree = load_rooted_tree_fixture(
        "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
    )
    dated_tree = load_rooted_tree_fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")

    report = evaluate_relaxed_lognormal_clock_tree_log_prior(
        substitution_tree,
        dated_tree,
        build_relaxed_lognormal_clock_model(
            rate_policy="independent",
            mean_clock_rate=0.2,
            log_standard_deviation=0.5,
        ),
    )

    row_by_taxa = {tuple(row.descendant_taxa): row for row in report.branch_rows}
    d_row = row_by_taxa[("D",)]
    assert math.isclose(d_row.dated_time_duration, 8.0, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        d_row.observed_substitution_branch_length,
        0.8,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(d_row.branch_rate, 0.1, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(d_row.anchor_rate, 0.2, rel_tol=0.0, abs_tol=1e-12)
    assert math.isclose(
        d_row.expected_substitution_branch_length,
        1.6,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        d_row.branch_rate_deviation,
        -0.1,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


@pytest.mark.parametrize(
    ("builder_kwargs", "message"),
    [
        (
            {
                "rate_policy": "unsupported",
                "mean_clock_rate": 0.2,
                "log_standard_deviation": 0.5,
            },
            "requires a supported rate policy",
        ),
        (
            {
                "rate_policy": "independent",
                "mean_clock_rate": 0.0,
                "log_standard_deviation": 0.5,
            },
            "requires a strictly positive finite mean clock rate",
        ),
        (
            {
                "rate_policy": "independent",
                "mean_clock_rate": 0.2,
                "log_standard_deviation": 0.0,
            },
            "requires a strictly positive finite log standard deviation",
        ),
    ],
)
def test_relaxed_lognormal_clock_prior_rejects_invalid_parameters(
    builder_kwargs: dict[str, float | str],
    message: str,
) -> None:
    with pytest.raises(PhylogeneticsError, match=message):
        build_relaxed_lognormal_clock_model(**builder_kwargs)


def test_relaxed_lognormal_clock_prior_rejects_unrooted_and_zero_duration_trees() -> (
    None
):
    substitution_tree = load_rooted_tree_fixture(
        "relaxed_rate_summary_substitution_tree_4_taxa.nwk"
    )
    dated_tree = load_rooted_tree_fixture("relaxed_rate_summary_dated_tree_4_taxa.nwk")
    unrooted_substitution_tree = substitution_tree.copy()
    unrooted_substitution_tree.rooted = False
    zero_duration_dated_tree = dated_tree.copy()
    zero_duration_dated_tree.root.children[1].branch_length = 0.0

    with pytest.raises(
        UnrootedTreeError, match="requires one rooted substitution tree"
    ):
        evaluate_relaxed_lognormal_clock_tree_log_prior(
            unrooted_substitution_tree,
            dated_tree,
            build_relaxed_lognormal_clock_model(
                rate_policy="independent",
                mean_clock_rate=0.2,
                log_standard_deviation=0.5,
            ),
        )

    with pytest.raises(
        InvalidBranchLengthError,
        match="requires strictly positive dated branch durations",
    ):
        evaluate_relaxed_lognormal_clock_tree_log_prior(
            substitution_tree,
            zero_duration_dated_tree,
            build_relaxed_lognormal_clock_model(
                rate_policy="independent",
                mean_clock_rate=0.2,
                log_standard_deviation=0.5,
            ),
        )

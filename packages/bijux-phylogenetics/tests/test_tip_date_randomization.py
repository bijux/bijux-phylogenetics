from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.diagnostics.root_to_tip import (
    diagnose_tip_date_randomization,
)
from bijux_phylogenetics.runtime.errors import PhylogeneticsError

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "metadata")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_tip_date_randomization_reports_reproducible_null_distribution_and_p_value() -> (
    None
):
    report = diagnose_tip_date_randomization(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
        permutations=19,
        seed=17,
    )

    assert report.tip_count == 7
    assert report.permutations == 19
    assert report.seed == 17
    assert report.p_value == pytest.approx(0.05, abs=1e-12)
    assert report.permuted_r_squared_at_or_above_observed == 0
    assert report.observed_regression.slope == pytest.approx(
        2.392857142857143,
        abs=1e-12,
    )
    assert report.observed_regression.r_squared == pytest.approx(
        0.6390945330296127,
        abs=1e-12,
    )
    assert report.null_distribution_minimum == pytest.approx(
        0.00014236902050102085,
        abs=1e-18,
    )
    assert report.null_distribution_mean == pytest.approx(
        0.15226741397913907,
        abs=1e-15,
    )
    assert report.null_distribution_maximum == pytest.approx(
        0.547266514806378,
        abs=1e-15,
    )
    assert len(report.permutation_rows) == 19
    assert report.permutation_rows[0].permuted_slope == pytest.approx(
        0.6071428571428571,
        abs=1e-15,
    )
    assert report.permutation_rows[0].permuted_r_squared == pytest.approx(
        0.041144646924829,
        abs=1e-15,
    )
    assert report.permutation_rows[-1].permuted_intercept == pytest.approx(
        0.10714285714285676,
        abs=1e-15,
    )


def test_tip_date_randomization_reuses_seeded_permutation_path() -> None:
    left = diagnose_tip_date_randomization(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
        permutations=11,
        seed=17,
    )
    right = diagnose_tip_date_randomization(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
        permutations=11,
        seed=17,
    )
    different_seed = diagnose_tip_date_randomization(
        fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
        fixture("root_to_tip_regression_dates_7_taxa.tsv"),
        permutations=11,
        seed=18,
    )

    assert left.seed == 17
    assert right.seed == 17
    assert left.permutation_rows == right.permutation_rows
    assert left.p_value == right.p_value
    assert left.null_distribution_mean == right.null_distribution_mean
    assert left.permutation_rows != different_seed.permutation_rows


def test_tip_date_randomization_requires_at_least_two_permutations() -> None:
    with pytest.raises(
        PhylogeneticsError,
        match="tip-date randomization requires at least two permutations",
    ):
        diagnose_tip_date_randomization(
            fixture("root_to_tip_regression_diagnostic_tree_7_taxa.nwk"),
            fixture("root_to_tip_regression_dates_7_taxa.tsv"),
            permutations=1,
        )

from __future__ import annotations

import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import solve_ctmc_stationary_distribution
from bijux_phylogenetics.phylo.likelihood.gtr import gtr_rate_matrix


def test_solve_ctmc_stationary_distribution_matches_two_state_fixture() -> None:
    report = solve_ctmc_stationary_distribution(
        [[-2.0, 2.0], [1.0, -1.0]],
        state_labels=("A", "B"),
    )

    assert report.state_labels == ("A", "B")
    assert report.probabilities == pytest.approx((1.0 / 3.0, 2.0 / 3.0))
    assert report.as_mapping() == pytest.approx({"A": 1.0 / 3.0, "B": 2.0 / 3.0})


def test_solve_ctmc_stationary_distribution_matches_asymmetric_three_state_fixture() -> None:
    report = solve_ctmc_stationary_distribution(
        [
            [-0.7, 0.6, 0.1],
            [1.2, -1.5, 0.3],
            [0.6, 0.9, -1.5],
        ],
        state_labels=("A", "B", "C"),
    )

    assert report.probabilities == pytest.approx((0.6, 0.3, 0.1), abs=1e-12)


def test_solve_ctmc_stationary_distribution_matches_four_state_fixture() -> None:
    rate_matrix = gtr_rate_matrix(
        {
            "AC": 1.5,
            "AG": 1.0,
            "AT": 2.0,
            "CG": 1.25,
            "CT": 5.0 / 3.0,
            "GT": 0.75,
        },
        base_frequencies={"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    )
    report = solve_ctmc_stationary_distribution(
        rate_matrix,
        state_labels=("A", "C", "G", "T"),
    )

    assert report.probabilities == pytest.approx((0.4, 0.1, 0.2, 0.3), abs=1e-12)

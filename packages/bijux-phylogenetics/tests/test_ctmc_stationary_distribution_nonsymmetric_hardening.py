from __future__ import annotations

import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution,
    verify_ctmc_stationary_distribution,
)


def test_solve_ctmc_stationary_distribution_matches_nonsymmetric_three_state_fixture() -> (
    None
):
    report = solve_ctmc_stationary_distribution(
        [
            [-3.0, 2.0, 1.0],
            [4.0, -5.0, 1.0],
            [2.0, 3.0, -5.0],
        ],
        state_labels=("A", "B", "C"),
    )

    assert report.probabilities == pytest.approx(
        (11.0 / 21.0, 13.0 / 42.0, 1.0 / 6.0),
        abs=1e-12,
    )


def test_verify_ctmc_stationary_distribution_accepts_nonsymmetric_mapping() -> None:
    report = verify_ctmc_stationary_distribution(
        [
            [-3.0, 2.0, 1.0],
            [4.0, -5.0, 1.0],
            [2.0, 3.0, -5.0],
        ],
        {
            "A": 11.0 / 21.0,
            "B": 13.0 / 42.0,
            "C": 1.0 / 6.0,
        },
        state_labels=("A", "B", "C"),
    )

    assert report.normalization_error == pytest.approx(0.0, abs=1e-12)
    assert max(abs(value) for value in report.residual_by_state) == pytest.approx(
        0.0,
        abs=1e-12,
    )

from __future__ import annotations

import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution,
)


def test_solve_ctmc_stationary_distribution_accepts_unique_absorbing_fixture() -> (
    None
):
    report = solve_ctmc_stationary_distribution(
        [
            [0.0, 0.0],
            [1.0, -1.0],
        ],
        state_labels=("absorbing", "transient"),
    )

    assert report.probabilities == pytest.approx((1.0, 0.0), abs=1e-12)
    assert report.as_mapping() == pytest.approx(
        {"absorbing": 1.0, "transient": 0.0},
        abs=1e-12,
    )

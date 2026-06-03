from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood.ctmc import (
    solve_ctmc_stationary_distribution,
    verify_ctmc_stationary_distribution,
)
from bijux_phylogenetics.phylo.likelihood.gtr import gtr_rate_matrix
from bijux_phylogenetics.runtime.errors import PhylogeneticsError


def test_solve_ctmc_stationary_distribution_matches_two_state_fixture() -> None:
    report = solve_ctmc_stationary_distribution(
        [[-2.0, 2.0], [1.0, -1.0]],
        state_labels=("A", "B"),
    )

    assert report.state_labels == ("A", "B")
    assert report.probabilities == pytest.approx((1.0 / 3.0, 2.0 / 3.0))
    assert report.as_mapping() == pytest.approx({"A": 1.0 / 3.0, "B": 2.0 / 3.0})


def test_solve_ctmc_stationary_distribution_matches_asymmetric_three_state_fixture() -> (
    None
):
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


def test_verify_ctmc_stationary_distribution_accepts_valid_mapping() -> None:
    report = verify_ctmc_stationary_distribution(
        [[-2.0, 2.0], [1.0, -1.0]],
        {"A": 1.0 / 3.0, "B": 2.0 / 3.0},
        state_labels=("A", "B"),
    )

    assert report.normalization_error == pytest.approx(0.0, abs=1e-12)
    assert max(abs(value) for value in report.residual_by_state) == pytest.approx(
        0.0,
        abs=1e-12,
    )


def test_verify_ctmc_stationary_distribution_rejects_non_stationary_vector() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        verify_ctmc_stationary_distribution(
            [[-2.0, 2.0], [1.0, -1.0]],
            {"A": 0.5, "B": 0.5},
            state_labels=("A", "B"),
        )

    assert error_info.value.code == "ctmc_stationary_distribution_residual_nonzero"


def test_verify_ctmc_stationary_distribution_rejects_negative_probability() -> None:
    with pytest.raises(PhylogeneticsError) as error_info:
        verify_ctmc_stationary_distribution(
            [[-2.0, 2.0], [1.0, -1.0]],
            numpy.array([-0.1, 1.1], dtype=float),
        )

    assert error_info.value.code == "ctmc_stationary_distribution_negative"

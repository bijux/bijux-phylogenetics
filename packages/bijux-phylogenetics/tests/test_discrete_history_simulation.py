from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.simulation import (
    DiscreteHistoryRateRow,
    simulate_discrete_histories,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    for group in ("trees", "alignments", "metadata", "expected"):
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    candidate = FIXTURES / name
    if candidate.exists():
        return candidate
    raise FileNotFoundError(name)


def _binary_rate_rows() -> list[DiscreteHistoryRateRow]:
    return [
        DiscreteHistoryRateRow(source_state="0", target_state="1", rate=0.4),
        DiscreteHistoryRateRow(source_state="1", target_state="0", rate=0.15),
    ]


def test_simulate_discrete_histories_records_transform_truth_and_changes_branch_scale() -> (
    None
):
    tree_path = fixture("example_tree.nwk")
    baseline = simulate_discrete_histories(
        tree_path,
        states=["0", "1"],
        rate_rows=_binary_rate_rows(),
        root_state="0",
        replicates=1,
        seed=7,
    )
    transformed = simulate_discrete_histories(
        tree_path,
        states=["0", "1"],
        rate_rows=_binary_rate_rows(),
        root_state="0",
        transform="kappa",
        transform_parameter_value=0.5,
        replicates=1,
        seed=7,
    )

    baseline_branch_length = sum(
        branch.branch_length for branch in baseline.simulations[0].branch_histories
    )
    transformed_branch_length = sum(
        branch.branch_length for branch in transformed.simulations[0].branch_histories
    )

    assert transformed.transform_name == "kappa"
    assert transformed.transform_parameter_name == "kappa"
    assert transformed.transform_parameter_value == 0.5
    assert transformed.simulations[0].transform_name == "kappa"
    assert transformed.simulations[0].transform_parameter_name == "kappa"
    assert transformed.simulations[0].transform_parameter_value == 0.5
    assert transformed_branch_length != pytest.approx(baseline_branch_length)


def test_simulate_discrete_histories_uses_a_for_early_burst_parameter_name() -> None:
    report = simulate_discrete_histories(
        fixture("example_tree.nwk"),
        states=["0", "1"],
        rate_rows=_binary_rate_rows(),
        root_state="0",
        transform="early-burst",
        transform_parameter_value=1.75,
        replicates=1,
        seed=11,
    )

    assert report.transform_name == "early-burst"
    assert report.transform_parameter_name == "a"
    assert report.transform_parameter_value == 1.75
    assert report.simulations[0].transform_parameter_name == "a"


def test_simulate_discrete_histories_rejects_incomplete_transform_request() -> None:
    with pytest.raises(
        ValueError,
        match="transform_parameter_value is required when a discrete-history transform is supplied",
    ):
        simulate_discrete_histories(
            fixture("example_tree.nwk"),
            states=["0", "1"],
            rate_rows=_binary_rate_rows(),
            root_state="0",
            transform="lambda",
            replicates=1,
            seed=5,
        )

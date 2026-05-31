from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood import (
    tree_log_likelihood as discrete_tree_log_likelihood,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_discrete_likelihood_observation_policies_change_likelihood_on_ambiguity_state() -> (
    None
):
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    rate_matrix = _three_state_rate_matrix()
    root_prior = numpy.full(3, 1.0 / 3.0, dtype=float)
    states_by_taxon = {"A": "0|1", "B": "2"}

    missing_report = discrete_tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=["0", "1", "2"],
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        observation_policy="treat-as-missing",
    )
    ambiguity_report = discrete_tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=["0", "1", "2"],
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        observation_policy="ambiguity-vector",
    )

    assert math.isfinite(missing_report)
    assert math.isfinite(ambiguity_report)
    assert not math.isclose(
        missing_report,
        ambiguity_report,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_likelihood_reject_policy_blocks_ambiguous_states() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))

    with pytest.raises(
        InvalidAlignmentError,
        match="does not allow ambiguous discrete state '0\\|1'",
    ):
        discrete_tree_log_likelihood(
            tree,
            {"A": "0|1", "B": "2"},
            state_order=["0", "1", "2"],
            rate_matrix=_three_state_rate_matrix(),
            root_prior=numpy.full(3, 1.0 / 3.0, dtype=float),
            observation_policy="reject",
        )


def test_discrete_likelihood_treat_as_missing_accepts_missing_token() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))

    log_likelihood = discrete_tree_log_likelihood(
        tree,
        {"A": "?", "B": "2"},
        state_order=["0", "1", "2"],
        rate_matrix=_three_state_rate_matrix(),
        root_prior=numpy.full(3, 1.0 / 3.0, dtype=float),
        observation_policy="treat-as-missing",
    )

    assert math.isfinite(log_likelihood)


def _three_state_rate_matrix() -> numpy.ndarray:
    return numpy.array(
        [
            [-1.0, 0.6, 0.4],
            [0.3, -0.8, 0.5],
            [0.2, 0.2, -0.4],
        ],
        dtype=float,
    )

from __future__ import annotations

import numpy
import pytest

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.finite_state_root_priors import (
    equal_finite_state_root_prior,
    estimate_empirical_finite_state_root_prior,
    fixed_state_finite_state_root_prior,
    resolve_finite_state_root_prior,
    validate_finite_state_root_prior,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

_STATE_ORDER = ("A", "B", "C")
_RECORDS = [
    AlignmentRecord(identifier="left", sequence="ABCA"),
    AlignmentRecord(identifier="right", sequence="AACC"),
]


def test_validate_finite_state_root_prior_normalizes_mapping_in_declared_order() -> None:
    root_prior = validate_finite_state_root_prior(
        {"B": 2.0, "A": 3.0, "C": 5.0},
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
    )

    assert numpy.allclose(root_prior, numpy.array([0.3, 0.2, 0.5], dtype=float))


def test_equal_and_fixed_state_root_priors_follow_declared_order() -> None:
    assert numpy.allclose(
        equal_finite_state_root_prior(_STATE_ORDER),
        numpy.array([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], dtype=float),
    )
    assert numpy.allclose(
        fixed_state_finite_state_root_prior(
            "B",
            state_labels=_STATE_ORDER,
            owner_name="finite-state root prior test",
        ),
        numpy.array([0.0, 1.0, 0.0], dtype=float),
    )


def test_estimate_empirical_finite_state_root_prior_ignores_unknown_tokens() -> None:
    root_prior = estimate_empirical_finite_state_root_prior(
        [
            AlignmentRecord(identifier="left", sequence="AXCA"),
            AlignmentRecord(identifier="right", sequence="?BCC"),
        ],
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
        sequence_tokenizer=tuple,
    )

    assert numpy.allclose(root_prior, numpy.array([2.0, 1.0, 3.0]) / 6.0)


def test_resolve_finite_state_root_prior_supports_all_shared_policies() -> None:
    equal = resolve_finite_state_root_prior(
        _RECORDS,
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
        default_policy="equal",
        root_prior_policy="equal",
        sequence_tokenizer=tuple,
    )
    empirical = resolve_finite_state_root_prior(
        _RECORDS,
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
        default_policy="equal",
        root_prior_policy="empirical",
        sequence_tokenizer=tuple,
    )
    stationary = resolve_finite_state_root_prior(
        _RECORDS,
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
        default_policy="equal",
        root_prior_policy="stationary",
        stationary_frequencies={"A": 0.2, "B": 0.5, "C": 0.3},
        sequence_tokenizer=tuple,
    )
    fixed = resolve_finite_state_root_prior(
        _RECORDS,
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
        default_policy="equal",
        root_prior_policy="fixed-state",
        fixed_root_state="C",
        sequence_tokenizer=tuple,
    )
    provided = resolve_finite_state_root_prior(
        _RECORDS,
        state_labels=_STATE_ORDER,
        owner_name="finite-state root prior test",
        default_policy="equal",
        root_prior_policy="user-supplied",
        root_prior={"A": 0.6, "B": 0.2, "C": 0.2},
        sequence_tokenizer=tuple,
    )

    assert equal.root_prior_source == "equal"
    assert empirical.root_prior_source == "empirical"
    assert stationary.root_prior_source == "stationary"
    assert fixed.root_prior_source == "fixed-state"
    assert provided.root_prior_source == "provided"

    assert numpy.allclose(empirical.root_prior, numpy.array([4.0, 1.0, 3.0]) / 8.0)
    assert numpy.allclose(stationary.root_prior, numpy.array([0.2, 0.5, 0.3]))
    assert numpy.allclose(fixed.root_prior, numpy.array([0.0, 0.0, 1.0]))


def test_resolve_finite_state_root_prior_rejects_missing_policy_inputs() -> None:
    with pytest.raises(ValueError, match="requires 'root_prior'"):
        resolve_finite_state_root_prior(
            _RECORDS,
            state_labels=_STATE_ORDER,
            owner_name="finite-state root prior test",
            default_policy="equal",
            root_prior_policy="provided",
            sequence_tokenizer=tuple,
        )

    with pytest.raises(InvalidAlignmentError, match="exactly the declared finite-state order"):
        validate_finite_state_root_prior(
            {"A": 1.0, "B": 1.0},
            state_labels=_STATE_ORDER,
            owner_name="finite-state root prior test",
        )

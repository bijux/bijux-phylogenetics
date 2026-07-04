from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.optimizer_search import (
    _fit_symmetric_warm_start,
    build_discrete_initial_candidates,
    optimize_single_log_parameter,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.parameter_layout import (
    broadcast_discrete_mk_log_parameters,
    build_discrete_mk_parameter_layout,
)
from bijux_phylogenetics.io.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_discrete_initial_candidates_stay_deterministic_for_unordered_ard() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    state_order = ["north", "south", "west"]
    states_by_taxon = {"A": "north", "B": "south"}
    allowed_transition_pairs = {
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 2),
        (2, 0),
        (2, 1),
    }

    first = build_discrete_initial_candidates(
        tree,
        states_by_taxon,
        state_order=state_order,
        parameter_count=6,
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode="equal",
        ascertainment_policy="none",
    )
    second = build_discrete_initial_candidates(
        tree,
        states_by_taxon,
        state_order=state_order,
        parameter_count=6,
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode="equal",
        ascertainment_policy="none",
    )

    assert len(first) == len(second)
    assert len(first) <= 10
    for (first_scale, first_seed), (second_scale, second_seed) in zip(
        first,
        second,
        strict=True,
    ):
        assert math.isclose(first_scale, second_scale, rel_tol=0.0, abs_tol=1e-12)
        assert numpy.allclose(first_seed, second_seed, rtol=0.0, atol=1e-12)


def test_symmetric_warm_start_skips_unsupported_directional_surface() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    state_order = ["north", "south", "west"]
    states_by_taxon = {"A": "north", "B": "south"}
    allowed_transition_pairs = {
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 2),
        (2, 1),
    }
    equal_layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model="equal-rates",
        state_ordering="unordered",
        allowed_transition_pairs=allowed_transition_pairs,
    )
    target_layout = build_discrete_mk_parameter_layout(
        state_order=state_order,
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs=allowed_transition_pairs,
    )
    equal_run = optimize_single_log_parameter(
        tree,
        states_by_taxon,
        state_order=state_order,
        model="equal-rates",
        state_ordering="unordered",
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode="equal",
        ascertainment_policy="none",
    )
    equal_seed = numpy.array(
        broadcast_discrete_mk_log_parameters(
            equal_layout,
            target_layout,
            tuple(float(value) for value in equal_run.log_parameters),
        ),
        dtype=float,
    )

    symmetric_seed = _fit_symmetric_warm_start(
        tree,
        states_by_taxon,
        state_order=state_order,
        state_ordering="unordered",
        allowed_transition_pairs=allowed_transition_pairs,
        root_prior_mode="equal",
        ascertainment_policy="none",
        equal_seed=equal_seed,
        target_layout=target_layout,
    )

    assert symmetric_seed is None

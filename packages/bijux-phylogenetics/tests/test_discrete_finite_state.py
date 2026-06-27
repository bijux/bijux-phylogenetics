from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    transition_probability_matrix,
    tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    compile_finite_state_likelihood_context,
    evaluate_finite_state_site_log_likelihood,
)
from tests.support.two_tip_likelihood_reference import (
    two_tip_log_likelihood_for_labels,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_unordered_discrete_surface_matches_native_mk_kernel() -> None:
    state_order = ["north", "south", "west"]
    states_by_taxon = {"A": "north", "B": "south"}
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.8], dtype=float)),
        state_order=state_order,
        model="equal-rates",
        state_ordering="unordered",
        allowed_transition_pairs={
            (0, 1),
            (0, 2),
            (1, 0),
            (1, 2),
            (2, 0),
            (2, 1),
        },
    )

    _assert_discrete_surface_matches_native_kernel(
        state_order=state_order,
        states_by_taxon=states_by_taxon,
        rate_matrix=rate_matrix,
        model_name="unordered discrete Mk",
    )
    _assert_discrete_surface_matches_two_tip_known_answer(
        state_order=state_order,
        states_by_taxon=states_by_taxon,
        rate_matrix=rate_matrix,
        model_name="unordered discrete Mk",
    )


def test_ordered_discrete_surface_matches_native_mk_kernel() -> None:
    state_order = ["north", "south", "island"]
    states_by_taxon = {"A": "north", "B": "island"}
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.6], dtype=float)),
        state_order=state_order,
        model="equal-rates",
        state_ordering="ordered",
        allowed_transition_pairs={
            (0, 1),
            (1, 0),
            (1, 2),
            (2, 1),
        },
    )

    _assert_discrete_surface_matches_native_kernel(
        state_order=state_order,
        states_by_taxon=states_by_taxon,
        rate_matrix=rate_matrix,
        model_name="ordered discrete Mk",
    )
    _assert_discrete_surface_matches_two_tip_known_answer(
        state_order=state_order,
        states_by_taxon=states_by_taxon,
        rate_matrix=rate_matrix,
        model_name="ordered discrete Mk",
    )


def test_irreversible_discrete_surface_matches_native_mk_kernel() -> None:
    state_order = ["present", "absent"]
    states_by_taxon = {"A": "present", "B": "absent"}
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.4], dtype=float)),
        state_order=state_order,
        model="all-rates-different",
        state_ordering="unordered",
        allowed_transition_pairs={(0, 1)},
    )

    _assert_discrete_surface_matches_native_kernel(
        state_order=state_order,
        states_by_taxon=states_by_taxon,
        rate_matrix=rate_matrix,
        model_name="irreversible discrete Mk",
    )
    _assert_discrete_surface_matches_two_tip_known_answer(
        state_order=state_order,
        states_by_taxon=states_by_taxon,
        rate_matrix=rate_matrix,
        model_name="irreversible discrete Mk",
    )


def _assert_discrete_surface_matches_native_kernel(
    *,
    state_order: list[str],
    states_by_taxon: dict[str, str],
    rate_matrix: numpy.ndarray,
    model_name: str,
) -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name=model_name,
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_site_log_likelihood(
        context,
        tuple(states_by_taxon[taxon] for taxon in ["A", "B"]),
        taxon_order=["A", "B"],
    )
    native_log_likelihood = tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )

    assert math.isclose(
        protocol_log_likelihood,
        native_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def _assert_discrete_surface_matches_two_tip_known_answer(
    *,
    state_order: list[str],
    states_by_taxon: dict[str, str],
    rate_matrix: numpy.ndarray,
    model_name: str,
) -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name=model_name,
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    expected_log_likelihood = two_tip_log_likelihood_for_labels(
        observed_pairs=((states_by_taxon["A"], states_by_taxon["B"]),),
        state_index={state: index for index, state in enumerate(state_order)},
        root_prior=root_prior,
        left_transition=transition_probability_matrix(rate_matrix, 1.0),
        right_transition=transition_probability_matrix(rate_matrix, 1.0),
    )

    protocol_log_likelihood = evaluate_finite_state_site_log_likelihood(
        context,
        tuple(states_by_taxon[taxon] for taxon in ["A", "B"]),
        taxon_order=["A", "B"],
    )
    native_log_likelihood = tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )

    assert math.isclose(
        protocol_log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        native_log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

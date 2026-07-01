from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
    evaluate_finite_state_tree_likelihood_with_invariant_mixture,
    evaluate_finite_state_tree_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_finite_state_invariant_mixture_keeps_invariant_mass_inside_likelihood() -> None:
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
    state_order = ["north", "south", "west"]
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=3,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("north", "north"),
                weight=3,
                site_positions=[1, 2, 3],
            )
        ],
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix_from_log_parameters(
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
        ),
        root_prior=numpy.full(len(state_order), 1.0 / len(state_order), dtype=float),
        model_name="unordered discrete Mk",
    )

    fixed_rate_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compressed_patterns,
        surface=surface,
    )
    invariant_report = evaluate_finite_state_tree_likelihood_with_invariant_mixture(
        tree,
        compressed_patterns,
        surface=surface,
        invariant_proportion=0.35,
    )

    assert invariant_report.model_name == "unordered discrete Mk"
    assert invariant_report.invariant_proportion == 0.35
    assert math.isclose(
        invariant_report.log_likelihood,
        sum(row.log_likelihood for row in invariant_report.site_likelihoods),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert invariant_report.log_likelihood > fixed_rate_log_likelihood

    first_row = invariant_report.site_likelihoods[0]
    assert math.isclose(
        first_row.invariant_component_likelihood,
        1.0 / 3.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert first_row.variable_component_likelihood > 0.0
    assert math.isclose(
        first_row.mixture_likelihood,
        (0.35 * first_row.invariant_component_likelihood)
        + (0.65 * first_row.variable_component_likelihood),
        rel_tol=0.0,
        abs_tol=1e-12,
    )

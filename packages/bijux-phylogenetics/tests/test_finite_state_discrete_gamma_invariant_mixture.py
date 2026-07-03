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
    evaluate_finite_state_tree_likelihood_with_discrete_gamma,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
)
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    return FIXTURES / "trees" / name


def test_finite_state_gamma_invariant_mixture_reports_components_and_boundary_hits() -> (
    None
):
    tree = load_tree(fixture("felsenstein_two_tip_tree.nwk"))
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
        state_order=["north", "south", "west"],
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8], dtype=float)),
            state_order=["north", "south", "west"],
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
        root_prior=numpy.full(3, 1.0 / 3.0, dtype=float),
        model_name="unordered discrete Mk",
    )

    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.8,
        category_count=4,
    )
    combined_report = (
        evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            compressed_patterns,
            surface=surface,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            lower_invariant_proportion_bound=0.35,
            upper_invariant_proportion_bound=0.9,
            emit_boundary_warnings=True,
        )
    )

    assert combined_report.alpha == 0.8
    assert combined_report.category_count == 4
    assert combined_report.invariant_proportion == 0.35
    assert combined_report.hit_lower_invariant_proportion_boundary is True
    assert combined_report.hit_upper_invariant_proportion_boundary is False
    assert combined_report.boundary_warnings == [
        "invariant proportion hit lower search boundary"
    ]
    assert not math.isclose(
        combined_report.log_likelihood,
        gamma_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    first_row = combined_report.site_likelihoods[0]
    assert len(first_row.category_likelihoods) == 4
    assert math.isclose(
        first_row.invariant_component_likelihood,
        1.0 / 3.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first_row.variable_component_likelihood,
        sum(
            category.weight * category_likelihood
            for category, category_likelihood in zip(
                combined_report.category_rates,
                first_row.category_likelihoods,
                strict=True,
            )
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first_row.mixture_likelihood,
        (0.35 * first_row.invariant_component_likelihood)
        + (0.65 * first_row.variable_component_likelihood),
        rel_tol=0.0,
        abs_tol=1e-12,
    )

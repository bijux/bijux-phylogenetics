from __future__ import annotations

from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
    build_discrete_mk_finite_state_surface,
    evaluate_discrete_mk_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
)

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_TREE = FIXTURES / "trees" / "felsenstein_two_tip_tree.nwk"


@pytest.mark.parametrize(
    ("model_name", "state_order", "rate_matrix", "records"),
    [
        (
            "unordered discrete Mk",
            ("north", "south", "west"),
            rate_matrix_from_log_parameters(
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
            [
                AlignmentRecord(identifier="A", sequence="north"),
                AlignmentRecord(identifier="B", sequence="north"),
            ],
        ),
        (
            "ordered discrete Mk",
            ("north", "south", "island"),
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.8, 0.35, 1.1, 0.7], dtype=float)),
                state_order=["north", "south", "island"],
                model="all-rates-different",
                state_ordering="ordered",
                allowed_transition_pairs={(0, 1), (1, 0), (1, 2), (2, 1)},
            ),
            [
                AlignmentRecord(identifier="A", sequence="north"),
                AlignmentRecord(identifier="B", sequence="north"),
            ],
        ),
        (
            "irreversible discrete Mk",
            ("present", "lost"),
            rate_matrix_from_log_parameters(
                numpy.log(numpy.array([0.8], dtype=float)),
                state_order=["present", "lost"],
                model="all-rates-different",
                state_ordering="ordered",
                allowed_transition_pairs={(0, 1)},
            ),
            [
                AlignmentRecord(identifier="A", sequence="present"),
                AlignmentRecord(identifier="B", sequence="present"),
            ],
        ),
    ],
)
def test_discrete_mk_discrete_gamma_invariant_mixture_matches_shared_surface(
    model_name: str,
    state_order: tuple[str, ...],
    rate_matrix: numpy.ndarray,
    records: list[AlignmentRecord],
) -> None:
    tree = load_tree(FIXTURE_TREE)

    report = evaluate_discrete_mk_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        state_order=state_order,
        rate_matrix=rate_matrix,
        model_name=model_name,
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        records=records,
        model_name=model_name,
    )
    shared_report = (
        evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            CompressedAlignmentSitePatterns(
                source_path=None,
                taxon_order=["A", "B"],
                alignment_length=1,
                pattern_count=1,
                patterns=[
                    AlignmentSitePattern(
                        pattern_id="pattern-1",
                        states=(records[0].sequence, records[1].sequence),
                        weight=1,
                        site_positions=[1],
                    )
                ],
            ),
            surface=surface,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
        )
    )

    assert report.model_name == model_name
    assert report.alpha == 0.8
    assert report.category_count == 4
    assert report.invariant_proportion == 0.35
    assert report.root_prior_source == surface.root_prior_source
    assert report.normalization_policy == surface.normalization_policy
    assert report.log_likelihood == shared_report.log_likelihood
    assert report.site_likelihoods[0].category_likelihoods == (
        shared_report.site_likelihoods[0].category_likelihoods
    )
    assert report.site_likelihoods[0].invariant_component_likelihood == (
        shared_report.site_likelihoods[0].invariant_component_likelihood
    )


def test_discrete_mk_discrete_gamma_invariant_mixture_emits_boundary_warnings() -> (
    None
):
    tree = load_tree(FIXTURE_TREE)
    records = [
        AlignmentRecord(identifier="A", sequence="north"),
        AlignmentRecord(identifier="B", sequence="north"),
    ]
    rate_matrix = rate_matrix_from_log_parameters(
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
    )

    report = evaluate_discrete_mk_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        state_order=("north", "south", "west"),
        rate_matrix=rate_matrix,
        model_name="unordered discrete Mk",
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        lower_invariant_proportion_bound=0.35,
        upper_invariant_proportion_bound=0.5,
        emit_boundary_warnings=True,
    )

    assert report.hit_lower_invariant_proportion_boundary is True
    assert report.hit_upper_invariant_proportion_boundary is False
    assert report.boundary_warnings == [
        "invariant proportion hit lower search boundary"
    ]

from __future__ import annotations

import math
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
    evaluate_discrete_mk_tree_likelihood_with_invariant_mixture,
    evaluate_finite_state_tree_log_likelihood,
)

FIXTURES = Path(__file__).parent / "fixtures"


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
def test_discrete_mk_invariant_mixture_runs_for_unordered_ordered_and_irreversible_surfaces(
    model_name: str,
    state_order: tuple[str, ...],
    rate_matrix: numpy.ndarray,
    records: list[AlignmentRecord],
) -> None:
    tree = load_tree(FIXTURES / "trees" / "felsenstein_two_tip_tree.nwk")
    report = evaluate_discrete_mk_tree_likelihood_with_invariant_mixture(
        tree,
        records,
        state_order=state_order,
        rate_matrix=rate_matrix,
        invariant_proportion=0.35,
        model_name=model_name,
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        records=records,
        model_name=model_name,
    )
    fixed_rate_log_likelihood = evaluate_finite_state_tree_log_likelihood(
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
    )

    assert report.model_name == model_name
    assert report.invariant_proportion == 0.35
    assert report.root_prior_source == surface.root_prior_source
    assert report.normalization_policy == surface.normalization_policy
    assert report.site_likelihoods[0].invariant_component_likelihood > 0.0
    assert not math.isclose(
        report.log_likelihood,
        fixed_rate_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

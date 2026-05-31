from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_empirical_protein_invariant_zero_matches_fixed_rate_surface() -> None:
    fixed_rate_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
    )
    zero_invariant_report = (
        evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments",
                "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            matrix_label="compact-polar",
            invariant_proportion=0.0,
        )
    )

    assert zero_invariant_report.invariant_proportion == 0.0
    assert math.isclose(
        zero_invariant_report.log_likelihood,
        fixed_rate_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert all(
        math.isclose(
            row.mixture_likelihood,
            row.variable_component_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in zero_invariant_report.site_likelihoods
    )
    assert all(
        math.isclose(
            row.log_likelihood,
            math.log(row.variable_component_likelihood),
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in zero_invariant_report.site_likelihoods
    )


def test_empirical_protein_gamma_invariant_zero_matches_gamma_surface() -> None:
    gamma_report = evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        alpha=0.8,
        category_count=4,
    )
    gamma_invariant_report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments",
                "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            matrix_label="compact-polar",
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.0,
        )
    )

    assert gamma_invariant_report.invariant_proportion == 0.0
    assert math.isclose(
        gamma_invariant_report.log_likelihood,
        gamma_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [row.category_likelihoods for row in gamma_invariant_report.site_likelihoods] == [
        row.category_likelihoods for row in gamma_report.site_likelihoods
    ]
    assert all(
        math.isclose(
            row.variable_component_likelihood,
            row.mixture_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in gamma_invariant_report.site_likelihoods
    )


def _compact_polar_rate_matrix() -> numpy.ndarray:
    return _build_empirical_rate_matrix(
        boosted_pairs={
            ("A", "C"): 0.45,
            ("C", "D"): 0.35,
            ("D", "E"): 0.55,
            ("A", "E"): 0.20,
        }
    )


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = _protein_state_index()
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _build_empirical_rate_matrix(
    *,
    boosted_pairs: dict[tuple[str, str], float],
) -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = _protein_state_index()
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _protein_state_order() -> tuple[str, ...]:
    return (
        "A",
        "C",
        "D",
        "E",
        "F",
        "G",
        "H",
        "I",
        "K",
        "L",
        "M",
        "N",
        "P",
        "Q",
        "R",
        "S",
        "T",
        "V",
        "W",
        "Y",
    )


def _protein_state_index() -> dict[str, int]:
    return {state: index for index, state in enumerate(_protein_state_order())}

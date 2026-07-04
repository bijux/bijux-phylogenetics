from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_discrete_gamma_empirical_protein_report_emits_categories_and_site_rows() -> (
    None
):
    report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=0.8,
            category_count=4,
            matrix_label="compact-polar",
        )
    )

    assert report.alpha == 0.8
    assert report.category_count == 4
    assert report.matrix_label == "compact-polar"
    assert report.root_prior_source == "provided"
    assert len(report.category_rates) == 4
    assert len(report.site_likelihoods) == 4
    assert [row.site_position for row in report.site_likelihoods] == [1, 2, 3, 4]
    assert all(len(row.category_likelihoods) == 4 for row in report.site_likelihoods)
    assert math.isclose(
        sum(category.weight for category in report.category_rates),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum(category.weight * category.rate for category in report.category_rates),
        1.0,
        rel_tol=0.0,
        abs_tol=1e-9,
    )
    assert math.isclose(
        report.log_likelihood,
        sum(row.log_likelihood for row in report.site_likelihoods),
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    for row in report.site_likelihoods:
        expected_mixture = sum(
            category.weight * category_likelihood
            for category, category_likelihood in zip(
                report.category_rates,
                row.category_likelihoods,
                strict=True,
            )
        )
        assert math.isclose(
            row.mixture_likelihood,
            expected_mixture,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        assert math.isclose(
            row.log_likelihood,
            math.log(row.mixture_likelihood),
            rel_tol=0.0,
            abs_tol=1e-12,
        )


def test_discrete_gamma_alpha_changes_empirical_protein_likelihood() -> None:
    low_alpha_report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=0.5,
            category_count=4,
            matrix_label="compact-polar",
        )
    )
    high_alpha_report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=8.0,
            category_count=4,
            matrix_label="compact-polar",
        )
    )

    assert not math.isclose(
        low_alpha_report.log_likelihood,
        high_alpha_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [
        round(category.rate, 12) for category in low_alpha_report.category_rates
    ] != [round(category.rate, 12) for category in high_alpha_report.category_rates]


def test_discrete_gamma_high_alpha_nearly_recovers_fixed_rate_empirical_likelihood() -> (
    None
):
    fixed_rate_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
    )
    gamma_report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=250.0,
            category_count=4,
            matrix_label="compact-polar",
        )
    )

    assert math.isclose(
        fixed_rate_report.log_likelihood,
        gamma_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=2e-3,
    )


def _compact_polar_rate_matrix():
    import numpy

    state_order = _protein_state_order()
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in {
        ("A", "C"): 0.45,
        ("C", "D"): 0.35,
        ("D", "E"): 0.55,
        ("A", "E"): 0.20,
    }.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(rate_matrix[row_index, :].sum())
    return rate_matrix


def _biased_root_prior():
    import numpy

    prior = numpy.full(20, 0.02, dtype=float)
    state_index = {state: index for index, state in enumerate(_protein_state_order())}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


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

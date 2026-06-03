from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment,
    evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_discrete_gamma_invariant_mixture_reports_alpha_invariant_and_total_likelihood() -> (
    None
):
    report = optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        alpha=0.8,
        category_count=4,
        matrix_label="compact-polar",
        initial_invariant_proportion=0.1,
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 20
    assert report.pattern_count == 9
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.state_count == 20
    assert report.matrix_label == "compact-polar"
    assert report.root_prior_source == "provided"
    assert report.gap_policy == "treat-as-missing"
    assert report.missing_policy == "treat-as-missing"
    assert report.alpha == 0.8
    assert report.category_count == 4
    assert [round(category.rate, 12) for category in report.category_rates] == [
        0.09825292034,
        0.442946818611,
        1.038681202727,
        2.420119058323,
    ]
    assert report.initial_invariant_proportion == 0.1
    assert math.isclose(
        report.invariant_proportion,
        0.10743071401974319,
        rel_tol=0.0,
        abs_tol=1e-7,
    )
    assert math.isclose(
        report.initial_log_likelihood,
        -58.343094924776175,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.log_likelihood,
        -58.342920438693604,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.function_evaluation_count == 46
    assert report.converged is True
    assert report.hit_lower_invariant_proportion_boundary is False
    assert report.hit_upper_invariant_proportion_boundary is False
    assert report.boundary_warnings == []
    assert len(report.site_likelihoods) == report.site_count


def test_discrete_gamma_invariant_mixture_rows_keep_gamma_and_invariant_components_active() -> (
    None
):
    report = evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        matrix_label="compact-polar",
    )

    assert report.alpha == 0.8
    assert report.invariant_proportion == 0.35
    assert report.boundary_warnings == []
    assert math.isclose(
        report.log_likelihood,
        -58.56155736229263,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.log_likelihood,
        sum(row.log_likelihood for row in report.site_likelihoods),
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    first_row = report.site_likelihoods[0]
    assert [round(value, 12) for value in first_row.category_likelihoods] == [
        0.194265051806,
        0.175656945764,
        0.148345599524,
        0.102714546823,
    ]
    assert math.isclose(
        first_row.invariant_component_likelihood,
        0.2,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first_row.variable_component_likelihood,
        0.15524553597929605,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        first_row.mixture_likelihood,
        0.1709095983865424,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    expected_gamma_mixture = sum(
        category.weight * category_likelihood
        for category, category_likelihood in zip(
            report.category_rates,
            first_row.category_likelihoods,
            strict=True,
        )
    )
    assert math.isclose(
        first_row.variable_component_likelihood,
        expected_gamma_mixture,
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


def test_discrete_gamma_invariant_mixture_boundary_warning_reports_bound_hit() -> None:
    report = optimize_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        alpha=0.8,
        category_count=4,
        matrix_label="compact-polar",
        initial_invariant_proportion=0.31,
        lower_invariant_proportion_bound=0.3,
        upper_invariant_proportion_bound=0.5,
    )

    assert math.isclose(
        report.invariant_proportion,
        0.3000000003338482,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.hit_lower_invariant_proportion_boundary is True
    assert report.hit_upper_invariant_proportion_boundary is False
    assert report.boundary_warnings == [
        "invariant proportion hit lower search boundary"
    ]


def test_discrete_gamma_invariant_mixture_is_not_a_standalone_report_merge() -> None:
    combined_report = evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        matrix_label="compact-polar",
    )
    gamma_report = (
        evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment(
            fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
            fixture(
                "alignments",
                "empirical_protein_invariant_mixture_alignment_2_taxa.fasta",
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            root_prior=_biased_root_prior(),
            alpha=0.8,
            category_count=4,
            matrix_label="compact-polar",
        )
    )
    invariant_report = evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        invariant_proportion=0.35,
        matrix_label="compact-polar",
    )

    assert math.isclose(
        gamma_report.log_likelihood,
        -58.37757694225371,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        invariant_report.log_likelihood,
        -58.78116686054642,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        combined_report.log_likelihood,
        gamma_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        combined_report.log_likelihood,
        invariant_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
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

from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
    optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_empirical_protein_invariant_mixture_estimates_nonzero_proportion() -> None:
    report = optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
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
    assert report.initial_invariant_proportion == 0.1
    assert math.isclose(
        report.invariant_proportion,
        0.21159721908532134,
        rel_tol=0.0,
        abs_tol=1e-7,
    )
    assert math.isclose(
        report.initial_log_likelihood,
        -58.74292047974127,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.log_likelihood,
        -58.69533356244757,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.log_likelihood > report.initial_log_likelihood
    assert report.function_evaluation_count == 46
    assert len(report.site_likelihoods) == report.site_count
    assert [row.site_position for row in report.site_likelihoods] == list(range(1, 21))


def test_empirical_protein_invariant_mixture_rows_emit_real_component_activity() -> (
    None
):
    report = evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        invariant_proportion=0.35,
    )

    assert report.invariant_proportion == 0.35
    assert report.initial_invariant_proportion == 0.35
    assert report.function_evaluation_count == 1
    assert report.converged is True
    assert math.isclose(
        report.log_likelihood,
        -58.78116686054642,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.log_likelihood,
        sum(row.log_likelihood for row in report.site_likelihoods),
        rel_tol=0.0,
        abs_tol=1e-12,
    )

    invariant_rows = report.site_likelihoods[:16]
    variable_rows = report.site_likelihoods[16:]
    assert all(row.invariant_component_likelihood > 0.0 for row in invariant_rows)
    assert all(row.invariant_component_likelihood == 0.0 for row in variable_rows)
    assert math.isclose(
        invariant_rows[0].invariant_component_likelihood,
        0.2,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        invariant_rows[0].variable_component_likelihood,
        0.14995334420240528,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        invariant_rows[0].mixture_likelihood,
        0.1674696737315634,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    for row in report.site_likelihoods:
        expected_mixture_likelihood = (0.35 * row.invariant_component_likelihood) + (
            0.65 * row.variable_component_likelihood
        )
        assert math.isclose(
            row.mixture_likelihood,
            expected_mixture_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        assert math.isclose(
            row.log_likelihood,
            math.log(row.mixture_likelihood),
            rel_tol=0.0,
            abs_tol=1e-12,
        )


def test_empirical_protein_invariant_mixture_improves_over_zero_invariant_case() -> (
    None
):
    zero_invariant_report = evaluate_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        invariant_proportion=0.0,
    )
    optimized_report = optimize_empirical_protein_tree_likelihood_with_invariant_mixture_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture(
            "alignments", "empirical_protein_invariant_mixture_alignment_2_taxa.fasta"
        ),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        initial_invariant_proportion=0.1,
    )

    assert zero_invariant_report.invariant_proportion == 0.0
    assert math.isclose(
        zero_invariant_report.log_likelihood,
        -58.85817779561978,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert (
        optimized_report.invariant_proportion
        > zero_invariant_report.invariant_proportion
    )
    assert optimized_report.log_likelihood > zero_invariant_report.log_likelihood


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

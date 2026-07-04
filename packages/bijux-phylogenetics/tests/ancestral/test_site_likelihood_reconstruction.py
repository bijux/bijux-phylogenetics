from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.comparative.discrete_mk import fit_discrete_mk_model
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_with_discrete_gamma_from_alignment,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_nucleotide_site_log_likelihoods_from_alignment,
    evaluate_protein_poisson_tree_likelihood_from_alignment,
    sum_site_log_likelihood_rows,
    sum_weighted_site_pattern_log_likelihood_rows,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_jc69_site_rows_reconstruct_total_likelihood() -> None:
    tree_path = fixture("trees", "likelihood_site_pattern_tree_4_taxa.nwk")
    alignment_path = fixture("alignments", "jc69_site_pattern_alignment.fasta")
    report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        tree_path,
        alignment_path,
        model_name="jc69",
    )
    direct = evaluate_jc69_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
    )

    assert math.isclose(report.log_likelihood, direct.log_likelihood, abs_tol=1e-12)
    assert math.isclose(
        sum_site_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum_weighted_site_pattern_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )


def test_hky85_site_rows_reconstruct_total_likelihood() -> None:
    tree_path = fixture("trees", "hky85_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        tree_path,
        alignment_path,
        model_name="hky85",
        kappa=4.0,
        base_frequencies=base_frequencies,
    )
    direct = evaluate_hky85_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        kappa=4.0,
        base_frequencies=base_frequencies,
    )

    assert math.isclose(report.log_likelihood, direct.log_likelihood, abs_tol=1e-12)
    assert math.isclose(
        sum_site_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum_weighted_site_pattern_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )


def test_gtr_site_rows_reconstruct_total_likelihood() -> None:
    tree_path = fixture("trees", "gtr_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
    base_frequencies = {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3}
    exchangeabilities = {
        "AC": 1.0,
        "AG": 4.5,
        "AT": 0.8,
        "CG": 1.6,
        "CT": 2.4,
        "GT": 3.1,
    }
    report = evaluate_nucleotide_site_log_likelihoods_from_alignment(
        tree_path,
        alignment_path,
        model_name="gtr",
        exchangeabilities=exchangeabilities,
        base_frequencies=base_frequencies,
    )
    direct = evaluate_gtr_tree_likelihood(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        exchangeabilities=exchangeabilities,
        base_frequencies=base_frequencies,
    )

    assert math.isclose(report.log_likelihood, direct.log_likelihood, abs_tol=1e-12)
    assert math.isclose(
        sum_site_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum_weighted_site_pattern_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )


def test_protein_poisson_site_rows_reconstruct_total_likelihood() -> None:
    report = evaluate_protein_poisson_tree_likelihood_from_alignment(
        fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"),
    )

    assert math.isclose(
        sum_site_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum_weighted_site_pattern_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )


def test_empirical_protein_gamma_site_rows_reconstruct_total_likelihood() -> None:
    report = (
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

    assert math.isclose(
        sum_site_log_likelihood_rows(report.site_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum_weighted_site_pattern_log_likelihood_rows(report.site_likelihoods),
        report.log_likelihood,
        abs_tol=1e-12,
    )


def test_discrete_mk_pattern_rows_reconstruct_total_likelihood() -> None:
    report = fit_discrete_mk_model(
        fixture("trees", "example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("metadata", "example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="binary_state",
        taxon_column="taxon",
        model="equal-rates",
    )

    assert math.isclose(
        math.fsum(
            row.pattern_weight * row.log_likelihood
            for row in report.pattern_likelihood_rows
        ),
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_mk_variable_only_pattern_rows_reconstruct_total_likelihood() -> None:
    report = fit_discrete_mk_model(
        fixture("trees", "example_tree_phytools_ultrametric_twenty_four_taxa.nwk"),
        fixture("metadata", "example_traits_phytools_signal_twenty_four_taxa.tsv"),
        trait="binary_state",
        taxon_column="taxon",
        model="equal-rates",
        ascertainment_policy="lewis-variable-only",
    )

    assert report.ascertainment_conditioning_log_probability is not None
    assert math.isclose(
        math.fsum(
            row.pattern_weight * row.log_likelihood
            for row in report.pattern_likelihood_rows
        ),
        report.log_likelihood,
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
    states = tuple("ACDEFGHIKLMNPQRSTVWY")
    index = {state: position for position, state in enumerate(states)}
    matrix = numpy.full((len(states), len(states)), 0.03, dtype=float)
    numpy.fill_diagonal(matrix, 0.0)
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = index[left_state]
        right_index = index[right_state]
        matrix[left_index, right_index] = rate
        matrix[right_index, left_index] = rate
    row_sums = matrix.sum(axis=1)
    numpy.fill_diagonal(matrix, -row_sums)
    return matrix


def _protein_state_index() -> dict[str, int]:
    return {state: index for index, state in enumerate("ACDEFGHIKLMNPQRSTVWY")}

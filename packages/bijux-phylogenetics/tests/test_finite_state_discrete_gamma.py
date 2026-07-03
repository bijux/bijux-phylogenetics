from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    build_codon_finite_state_surface,
    build_discrete_mk_finite_state_surface,
    build_protein_poisson_finite_state_surface,
    build_selected_nucleotide_finite_state_surface,
    compress_alignment_site_patterns_from_records,
    compress_codon_site_patterns_from_records,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma,
    evaluate_finite_state_tree_log_likelihood,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_finite_state_discrete_gamma_reports_shared_jc69_categories() -> None:
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="JC69",
        owner_name="finite-state shared gamma proof",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)

    report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.8,
        category_count=4,
    )
    fixed_rate_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compressed_patterns,
        surface=surface,
    )
    high_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=250.0,
        category_count=4,
    )

    assert report.model_name == "JC69"
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.state_count == 4
    assert report.category_count == 4
    assert len(report.category_rates) == 4
    assert len(report.site_likelihoods) == 4
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
    assert math.isclose(
        high_alpha_report.log_likelihood,
        fixed_rate_log_likelihood,
        rel_tol=0.0,
        abs_tol=2e-3,
    )


def test_finite_state_discrete_gamma_alpha_changes_codon_likelihood() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    surface = build_codon_finite_state_surface(records)
    compressed_patterns = compress_codon_site_patterns_from_records(records)

    low_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.5,
        category_count=4,
    )
    high_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=8.0,
        category_count=4,
    )

    assert low_alpha_report.model_name == "codon CTMC"
    assert low_alpha_report.state_count == 61
    assert not math.isclose(
        low_alpha_report.log_likelihood,
        high_alpha_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert [
        round(category.rate, 12) for category in low_alpha_report.category_rates
    ] != [round(category.rate, 12) for category in high_alpha_report.category_rates]


def test_finite_state_discrete_gamma_alpha_changes_protein_poisson_likelihood() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    surface = build_protein_poisson_finite_state_surface(records)
    compressed_patterns = compress_alignment_site_patterns_from_records(records)

    low_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.5,
        category_count=4,
    )
    high_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=8.0,
        category_count=4,
    )

    assert low_alpha_report.model_name == "protein Poisson"
    assert low_alpha_report.state_count == 20
    assert not math.isclose(
        low_alpha_report.log_likelihood,
        high_alpha_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_finite_state_discrete_gamma_alpha_changes_discrete_mk_likelihood() -> None:
    tree = load_tree(fixture("trees", "example_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="0"),
        AlignmentRecord(identifier="B", sequence="0"),
        AlignmentRecord(identifier="C", sequence="1"),
        AlignmentRecord(identifier="D", sequence="1"),
    ]
    rate_matrix = numpy.asarray(((-1.0, 1.0), (1.0, -1.0)), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=("0", "1"),
        rate_matrix=rate_matrix,
        records=records,
        root_prior_policy="equal",
        model_name="discrete Mk",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        records,
        source_path=fixture(
            "metadata", "example_traits_discrete_mk_variable_only_four_taxa.tsv"
        ),
        site_tokenizer=lambda sequence: (sequence,),
    )

    low_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.5,
        category_count=4,
    )
    high_alpha_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=8.0,
        category_count=4,
    )

    assert low_alpha_report.model_name == "discrete Mk"
    assert low_alpha_report.state_count == 2
    assert len(low_alpha_report.site_likelihoods) == 1
    assert not math.isclose(
        low_alpha_report.log_likelihood,
        high_alpha_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

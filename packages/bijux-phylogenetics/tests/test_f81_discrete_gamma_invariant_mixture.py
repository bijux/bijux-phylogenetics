from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_selected_nucleotide_finite_state_surface,
    compress_alignment_site_patterns_from_records,
    evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
)

FIXTURES = Path(__file__).parent / "fixtures"

BASE_FREQUENCIES = {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15}


def test_f81_discrete_gamma_invariant_mixture_matches_owned_surface_expectations() -> (
    None
):
    tree = load_tree(FIXTURES / "trees" / "f81_likelihood_tree_2_taxa.nwk")
    records = load_fasta_alignment(
        FIXTURES / "alignments" / "f81_likelihood_alignment_2_taxa.fasta"
    )
    report = evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        base_frequencies=BASE_FREQUENCIES,
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="f81",
        owner_name="F81 likelihood",
        base_frequencies=BASE_FREQUENCIES,
    )
    shared_report = (
        evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            compress_alignment_site_patterns_from_records(records),
            surface=surface,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
        )
    )

    assert report.model_name == surface.model_name
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


def test_f81_discrete_gamma_invariant_mixture_emits_boundary_warnings() -> None:
    tree = load_tree(FIXTURES / "trees" / "f81_likelihood_tree_2_taxa.nwk")
    records = load_fasta_alignment(
        FIXTURES / "alignments" / "f81_likelihood_alignment_2_taxa.fasta"
    )

    report = evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        tree,
        records,
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        base_frequencies=BASE_FREQUENCIES,
        lower_invariant_proportion_bound=0.35,
        upper_invariant_proportion_bound=0.5,
        emit_boundary_warnings=True,
    )

    assert report.hit_lower_invariant_proportion_boundary is True
    assert report.hit_upper_invariant_proportion_boundary is False
    assert report.boundary_warnings == [
        "invariant proportion hit lower search boundary"
    ]


def test_f81_discrete_gamma_invariant_mixture_from_alignment_matches_direct_call() -> (
    None
):
    tree_path = FIXTURES / "trees" / "f81_likelihood_tree_2_taxa.nwk"
    alignment_path = FIXTURES / "alignments" / "f81_likelihood_alignment_2_taxa.fasta"

    from_alignment = (
        evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
            tree_path,
            alignment_path,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            base_frequencies=BASE_FREQUENCIES,
        )
    )
    direct = evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        base_frequencies=BASE_FREQUENCIES,
    )

    assert from_alignment.log_likelihood == direct.log_likelihood
    assert from_alignment.site_likelihoods[0].mixture_likelihood == (
        direct.site_likelihoods[0].mixture_likelihood
    )

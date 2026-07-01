from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_selected_nucleotide_finite_state_surface,
    compress_alignment_site_patterns_from_records,
    evaluate_finite_state_tree_log_likelihood,
    evaluate_hky85_tree_likelihood_with_invariant_mixture,
    evaluate_hky85_tree_likelihood_with_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"

BASE_FREQUENCIES = {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15}


def test_hky85_invariant_mixture_matches_owned_surface_expectations() -> None:
    tree = load_tree(FIXTURES / "trees" / "hky85_likelihood_tree_2_taxa.nwk")
    records = load_fasta_alignment(
        FIXTURES / "alignments" / "hky85_likelihood_alignment_2_taxa.fasta"
    )
    report = evaluate_hky85_tree_likelihood_with_invariant_mixture(
        tree,
        records,
        kappa=4.0,
        invariant_proportion=0.35,
        base_frequencies=BASE_FREQUENCIES,
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="hky85",
        owner_name="HKY85 likelihood",
        kappa=4.0,
        base_frequencies=BASE_FREQUENCIES,
    )
    fixed_rate_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compress_alignment_site_patterns_from_records(records),
        surface=surface,
    )

    assert report.model_name == surface.model_name
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


def test_hky85_invariant_mixture_from_alignment_matches_direct_call() -> None:
    tree_path = FIXTURES / "trees" / "hky85_likelihood_tree_2_taxa.nwk"
    alignment_path = (
        FIXTURES / "alignments" / "hky85_likelihood_alignment_2_taxa.fasta"
    )

    from_alignment = (
        evaluate_hky85_tree_likelihood_with_invariant_mixture_from_alignment(
            tree_path,
            alignment_path,
            kappa=4.0,
            invariant_proportion=0.35,
            base_frequencies=BASE_FREQUENCIES,
        )
    )
    direct = evaluate_hky85_tree_likelihood_with_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        kappa=4.0,
        invariant_proportion=0.35,
        base_frequencies=BASE_FREQUENCIES,
    )

    assert from_alignment.log_likelihood == direct.log_likelihood
    assert from_alignment.site_likelihoods[0].mixture_likelihood == (
        direct.site_likelihoods[0].mixture_likelihood
    )

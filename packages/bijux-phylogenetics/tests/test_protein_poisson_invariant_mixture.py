from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    build_protein_poisson_finite_state_surface,
    compress_alignment_site_patterns_from_records,
    evaluate_finite_state_tree_log_likelihood,
    evaluate_protein_poisson_tree_likelihood_with_invariant_mixture,
    evaluate_protein_poisson_tree_likelihood_with_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def test_protein_poisson_invariant_mixture_matches_shared_surface_expectations() -> (
    None
):
    tree = load_tree(FIXTURES / "trees" / "protein_poisson_likelihood_tree_2_taxa.nwk")
    records = load_fasta_alignment(
        FIXTURES / "alignments" / "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )
    report = evaluate_protein_poisson_tree_likelihood_with_invariant_mixture(
        tree,
        records,
        invariant_proportion=0.35,
    )
    surface = build_protein_poisson_finite_state_surface(records)
    fixed_rate_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compress_alignment_site_patterns_from_records(records),
        surface=surface,
    )

    assert report.model_name == "protein Poisson"
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


def test_protein_poisson_invariant_mixture_from_alignment_matches_direct_call() -> None:
    tree_path = FIXTURES / "trees" / "protein_poisson_likelihood_tree_2_taxa.nwk"
    alignment_path = (
        FIXTURES / "alignments" / "protein_poisson_likelihood_alignment_2_taxa.fasta"
    )

    from_alignment = (
        evaluate_protein_poisson_tree_likelihood_with_invariant_mixture_from_alignment(
            tree_path,
            alignment_path,
            invariant_proportion=0.35,
        )
    )
    direct = evaluate_protein_poisson_tree_likelihood_with_invariant_mixture(
        load_tree(tree_path),
        load_fasta_alignment(alignment_path),
        invariant_proportion=0.35,
    )

    assert from_alignment.log_likelihood == direct.log_likelihood
    assert from_alignment.site_likelihoods[0].mixture_likelihood == (
        direct.site_likelihoods[0].mixture_likelihood
    )

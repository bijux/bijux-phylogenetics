from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_gtr_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_gtr_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_hky85_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_hky85_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_jc69_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_jc69_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_k80_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_k80_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
    evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"

BASE_FREQUENCIES = {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15}
EXCHANGEABILITIES = {
    ("A", "C"): 0.4,
    ("A", "G"): 2.0,
    ("A", "T"): 0.6,
    ("C", "G"): 0.5,
    ("C", "T"): 1.6,
    ("G", "T"): 0.7,
}


@pytest.mark.parametrize(
    (
        "model_name",
        "owner_name",
        "tree_name",
        "alignment_name",
        "owned_wrapper",
        "owned_from_alignment_wrapper",
        "wrapper_kwargs",
    ),
    [
        (
            "jc69",
            "JC69 likelihood",
            "jc69_likelihood_tree_2_taxa.nwk",
            "jc69_likelihood_alignment_2_taxa.fasta",
            evaluate_jc69_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
            evaluate_jc69_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
            {},
        ),
        (
            "k80",
            "K80 likelihood",
            "jc69_likelihood_tree_2_taxa.nwk",
            "jc69_likelihood_alignment_2_taxa.fasta",
            evaluate_k80_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
            evaluate_k80_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
            {"kappa": 3.2},
        ),
        (
            "f81",
            "F81 likelihood",
            "f81_likelihood_tree_2_taxa.nwk",
            "f81_likelihood_alignment_2_taxa.fasta",
            evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
            evaluate_f81_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
            {"base_frequencies": BASE_FREQUENCIES},
        ),
        (
            "hky85",
            "HKY85 likelihood",
            "hky85_likelihood_tree_2_taxa.nwk",
            "hky85_likelihood_alignment_2_taxa.fasta",
            evaluate_hky85_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
            evaluate_hky85_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
            {"kappa": 4.0, "base_frequencies": BASE_FREQUENCIES},
        ),
        (
            "gtr",
            "GTR likelihood",
            "gtr_likelihood_tree_2_taxa.nwk",
            "gtr_likelihood_alignment_2_taxa.fasta",
            evaluate_gtr_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
            evaluate_gtr_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
            {
                "base_frequencies": BASE_FREQUENCIES,
                "exchangeabilities": EXCHANGEABILITIES,
            },
        ),
    ],
)
def test_owned_nucleotide_discrete_gamma_invariant_wrappers_match_selected_surface(
    model_name: str,
    owner_name: str,
    tree_name: str,
    alignment_name: str,
    owned_wrapper,
    owned_from_alignment_wrapper,
    wrapper_kwargs: dict[str, object],
) -> None:
    tree_path = FIXTURES / "trees" / tree_name
    alignment_path = FIXTURES / "alignments" / alignment_name
    tree = load_tree(tree_path)
    records = load_fasta_alignment(alignment_path)

    owned_report = owned_wrapper(
        tree,
        records,
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        **wrapper_kwargs,
    )
    selected_report = (
        evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            records,
            model_name=model_name,
            owner_name=owner_name,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            **wrapper_kwargs,
        )
    )

    assert owned_report.log_likelihood == selected_report.log_likelihood
    assert owned_report.root_prior_source == selected_report.root_prior_source
    assert owned_report.normalization_policy == selected_report.normalization_policy
    assert owned_report.site_likelihoods[0].category_likelihoods == (
        selected_report.site_likelihoods[0].category_likelihoods
    )
    assert owned_report.site_likelihoods[0].invariant_component_likelihood == (
        selected_report.site_likelihoods[0].invariant_component_likelihood
    )

    owned_from_alignment_report = owned_from_alignment_wrapper(
        tree_path,
        alignment_path,
        alpha=0.8,
        category_count=4,
        invariant_proportion=0.35,
        **wrapper_kwargs,
    )
    selected_from_alignment_report = (
        evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
            tree_path,
            alignment_path,
            model_name=model_name,
            owner_name=owner_name,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            **wrapper_kwargs,
        )
    )

    assert owned_from_alignment_report.log_likelihood == (
        selected_from_alignment_report.log_likelihood
    )
    assert owned_from_alignment_report.site_likelihoods[0].mixture_likelihood == (
        selected_from_alignment_report.site_likelihoods[0].mixture_likelihood
    )

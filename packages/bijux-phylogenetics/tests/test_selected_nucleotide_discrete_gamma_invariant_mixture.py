from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    build_selected_nucleotide_finite_state_surface,
    compress_alignment_site_patterns_from_records,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
    evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize(
    ("model_name", "kwargs"),
    [
        ("jc69", {}),
        ("k80", {"kappa": 3.2}),
        (
            "f81",
            {"base_frequencies": {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15}},
        ),
        (
            "hky85",
            {
                "kappa": 4.0,
                "base_frequencies": {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
            },
        ),
        (
            "gtr",
            {
                "base_frequencies": {"A": 0.55, "C": 0.15, "G": 0.15, "T": 0.15},
                "exchangeabilities": {
                    ("A", "C"): 0.4,
                    ("A", "G"): 2.0,
                    ("A", "T"): 0.6,
                    ("C", "G"): 0.5,
                    ("C", "T"): 1.6,
                    ("G", "T"): 0.7,
                },
            },
        ),
    ],
)
def test_selected_nucleotide_discrete_gamma_invariant_mixture_matches_shared_surface(
    model_name: str,
    kwargs: dict[str, object],
) -> None:
    tree = load_tree(FIXTURES / "trees" / "jc69_likelihood_tree_2_taxa.nwk")
    records = [
        AlignmentRecord(identifier="A", sequence="AAAA"),
        AlignmentRecord(identifier="B", sequence="AAAA"),
    ]

    report = (
        evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            records,
            model_name=model_name,
            owner_name="selected nucleotide gamma-invariant test",
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            **kwargs,
        )
    )
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name=model_name,
        owner_name="selected nucleotide gamma-invariant test",
        **kwargs,
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


def test_selected_nucleotide_discrete_gamma_invariant_mixture_emits_boundary_warnings() -> (
    None
):
    tree = load_tree(FIXTURES / "trees" / "jc69_likelihood_tree_2_taxa.nwk")
    records = [
        AlignmentRecord(identifier="A", sequence="AAAA"),
        AlignmentRecord(identifier="B", sequence="AAAA"),
    ]

    report = (
        evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            records,
            model_name="jc69",
            owner_name="selected nucleotide gamma-invariant test",
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
            lower_invariant_proportion_bound=0.35,
            upper_invariant_proportion_bound=0.5,
            emit_boundary_warnings=True,
        )
    )

    assert report.hit_lower_invariant_proportion_boundary is True
    assert report.hit_upper_invariant_proportion_boundary is False
    assert report.boundary_warnings == [
        "invariant proportion hit lower search boundary"
    ]


def test_selected_nucleotide_discrete_gamma_invariant_mixture_from_alignment_matches_direct_call() -> (
    None
):
    tree_path = FIXTURES / "trees" / "jc69_likelihood_tree_2_taxa.nwk"
    alignment_path = FIXTURES / "alignments" / "jc69_likelihood_alignment_2_taxa.fasta"

    from_alignment = (
        evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture_from_alignment(
            tree_path,
            alignment_path,
            model_name="jc69",
            owner_name="selected nucleotide gamma-invariant test",
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
        )
    )
    direct = (
        evaluate_selected_nucleotide_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            load_tree(tree_path),
            load_fasta_alignment(alignment_path),
            model_name="jc69",
            owner_name="selected nucleotide gamma-invariant test",
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
        )
    )

    assert from_alignment.log_likelihood == direct.log_likelihood
    assert from_alignment.site_likelihoods[0].mixture_likelihood == (
        direct.site_likelihoods[0].mixture_likelihood
    )

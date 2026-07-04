from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma,
    evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    compress_codon_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
    build_codon_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    build_discrete_mk_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_named_empirical_protein_finite_state_surface,
    build_protein_poisson_finite_state_surface,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("family_name", "case_factory"),
    [
        ("nucleotide", lambda: _nucleotide_case()),
        ("protein-poisson", lambda: _protein_poisson_case()),
        ("protein-empirical", lambda: _protein_empirical_case()),
        ("codon", lambda: _codon_case()),
        ("discrete-unordered", lambda: _discrete_unordered_case()),
    ],
)
def test_gamma_invariant_mixture_changes_likelihood_for_every_eligible_finite_state_family(
    family_name: str,
    case_factory,
) -> None:
    tree, compressed_patterns, surface = case_factory()

    gamma_report = evaluate_finite_state_tree_likelihood_with_discrete_gamma(
        tree,
        compressed_patterns,
        surface=surface,
        alpha=0.8,
        category_count=4,
    )
    combined_report = (
        evaluate_finite_state_tree_likelihood_with_discrete_gamma_and_invariant_mixture(
            tree,
            compressed_patterns,
            surface=surface,
            alpha=0.8,
            category_count=4,
            invariant_proportion=0.35,
        )
    )

    assert not math.isclose(
        combined_report.log_likelihood,
        gamma_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), family_name
    assert combined_report.site_likelihoods
    assert combined_report.site_likelihoods[0].invariant_component_likelihood > 0.0
    assert combined_report.site_likelihoods[0].variable_component_likelihood > 0.0


def _nucleotide_case():
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="AAAA"),
        AlignmentRecord(identifier="B", sequence="AAAA"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="jc69",
        owner_name="nucleotide gamma-invariant mixture test",
    )
    return tree, compressed_patterns, surface


def _protein_poisson_case():
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="AAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAA"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(records)
    return tree, compressed_patterns, surface


def _protein_empirical_case():
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="AAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAA"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_named_empirical_protein_finite_state_surface(
        records,
        model_name="WAG",
    )
    return tree, compressed_patterns, surface


def _codon_case():
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="AAAAAA"),
        AlignmentRecord(identifier="B", sequence="AAAAAA"),
    ]
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    surface = build_codon_finite_state_surface(records)
    return tree, compressed_patterns, surface


def _discrete_unordered_case():
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=3,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("north", "north"),
                weight=3,
                site_positions=[1, 2, 3],
            )
        ],
    )
    surface = build_discrete_mk_finite_state_surface(
        state_order=["north", "south", "west"],
        rate_matrix=rate_matrix_from_log_parameters(
            numpy.log(numpy.array([0.8], dtype=float)),
            state_order=["north", "south", "west"],
            model="equal-rates",
            state_ordering="unordered",
            allowed_transition_pairs={
                (0, 1),
                (0, 2),
                (1, 0),
                (1, 2),
                (2, 0),
                (2, 1),
            },
        ),
        root_prior=numpy.full(3, 1.0 / 3.0, dtype=float),
        model_name="unordered discrete Mk",
    )
    return tree, compressed_patterns, surface

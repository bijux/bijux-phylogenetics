from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    sum_site_log_likelihood_rows,
    sum_weighted_site_pattern_log_likelihood_rows,
)
from bijux_phylogenetics.phylo.likelihood.codon import (
    compress_codon_site_patterns_from_records,
    evaluate_codon_ctmc_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
    build_codon_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import resolve_codon_state_space
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    compile_finite_state_likelihood_context,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_codon_surface_matches_existing_tree_report() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    surface = build_codon_finite_state_surface(records)
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_codon_ctmc_tree_likelihood(tree, records)

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_surface_respects_ambiguity_vector_policy() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_ambiguity_2_taxa.fasta")
    )
    compressed_patterns = compress_codon_site_patterns_from_records(
        records,
        observation_policy="ambiguity-vector",
    )
    surface = build_codon_finite_state_surface(
        records,
        observation_policy="ambiguity-vector",
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        observation_policy="ambiguity-vector",
    )

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_surface_respects_custom_frequency_surface() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    state_space = resolve_codon_state_space()
    biased_frequencies = dict.fromkeys(state_space.state_order, 1.0)
    biased_frequencies["AAA"] = 20.0
    biased_frequencies["AAG"] = 15.0
    biased_frequencies["GCT"] = 0.5
    biased_frequencies["GCC"] = 0.5
    surface = build_codon_finite_state_surface(
        records,
        codon_frequencies=biased_frequencies,
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        codon_frequencies=biased_frequencies,
    )

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_surface_preserves_partial_missing_and_gap_constraints() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="AT?A-G"),
        AlignmentRecord(identifier="B", sequence="ATGAAG"),
    ]
    compressed_patterns = compress_codon_site_patterns_from_records(
        records,
        observation_policy="ambiguity-vector",
    )
    surface = build_codon_finite_state_surface(
        records,
        observation_policy="ambiguity-vector",
    )
    context = compile_finite_state_likelihood_context(tree, surface)

    protocol_log_likelihood = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )
    family_report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        observation_policy="ambiguity-vector",
    )
    missing_report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        observation_policy="treat-as-missing",
    )

    assert math.isclose(
        protocol_log_likelihood,
        family_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert not math.isclose(
        family_report.log_likelihood,
        missing_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_codon_surface_normalizes_rna_records_for_empirical_root_prior() -> None:
    records = [
        AlignmentRecord(identifier="A", sequence="augaaa"),
        AlignmentRecord(identifier="B", sequence="augaag"),
    ]

    surface = build_codon_finite_state_surface(
        records,
        root_prior_policy="empirical",
    )

    assert surface.root_prior_source == "empirical"
    assert surface.state_count == 61
    assert math.isclose(float(numpy.sum(surface.root_prior)), 1.0, abs_tol=1e-12)
    assert math.isclose(
        float(surface.root_prior[surface.state_labels.index("ATG")]),
        0.5,
        abs_tol=1e-12,
    )
    assert math.isclose(
        float(surface.root_prior[surface.state_labels.index("AAA")]),
        0.25,
        abs_tol=1e-12,
    )
    assert math.isclose(
        float(surface.root_prior[surface.state_labels.index("AAG")]),
        0.25,
        abs_tol=1e-12,
    )


def test_codon_tree_report_expands_site_rows_in_original_order() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_site_pattern_alignment_2_taxa.fasta")
    )

    report = evaluate_codon_ctmc_tree_likelihood(tree, records)

    assert report.expansion_policy == "expanded-site-rows"
    assert report.site_count == 6
    assert report.pattern_count == 3
    assert [row.site_position for row in report.site_log_likelihoods] == [
        1,
        2,
        3,
        4,
        5,
        6,
    ]
    assert [row.pattern_id for row in report.site_log_likelihoods] == [
        "pattern-1",
        "pattern-2",
        "pattern-1",
        "pattern-2",
        "pattern-3",
        "pattern-3",
    ]
    assert [row.pattern_weight for row in report.site_log_likelihoods] == [
        2,
        2,
        2,
        2,
        2,
        2,
    ]
    assert [row.site_states for row in report.site_log_likelihoods] == [
        ("AAA", "AAG"),
        ("GCT", "GCC"),
        ("AAA", "AAG"),
        ("GCT", "GCC"),
        ("TTT", "TTT"),
        ("TTT", "TTT"),
    ]
    assert math.isclose(
        sum_site_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        sum_weighted_site_pattern_log_likelihood_rows(report.site_log_likelihoods),
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

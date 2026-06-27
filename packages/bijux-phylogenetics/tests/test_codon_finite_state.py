from __future__ import annotations

import math
from pathlib import Path

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
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

from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.likelihood import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
    compress_alignment_site_patterns_from_records,
    evaluate_codon_ctmc_tree_likelihood,
    evaluate_empirical_protein_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_protein_poisson_tree_likelihood,
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
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    evaluate_finite_state_tree_log_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.nucleotide_finite_state import (
    build_selected_nucleotide_finite_state_surface,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    build_empirical_protein_finite_state_surface,
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
        ("discrete-ordered", lambda: _discrete_ordered_case()),
        ("discrete-irreversible", lambda: _discrete_irreversible_case()),
    ],
)
def test_one_finite_state_evaluator_spans_every_goal_one_model_family(
    family_name: str,
    case_factory,
) -> None:
    tree, compressed_patterns, surface, expected_log_likelihood = case_factory()

    observed_log_likelihood = evaluate_finite_state_tree_log_likelihood(
        tree,
        compressed_patterns,
        surface=surface,
    )

    assert math.isclose(
        observed_log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), family_name


def _nucleotide_case():
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_selected_nucleotide_finite_state_surface(
        records,
        model_name="jc69",
        owner_name="family unification test",
    )
    expected = evaluate_jc69_tree_likelihood(tree, records).log_likelihood
    return tree, compressed_patterns, surface, expected


def _protein_poisson_case():
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    surface = build_protein_poisson_finite_state_surface(records)
    expected = evaluate_protein_poisson_tree_likelihood(tree, records).log_likelihood
    return tree, compressed_patterns, surface, expected


def _protein_empirical_case():
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(records)
    rate_matrix = _compact_polar_rate_matrix()
    root_prior = _biased_root_prior()
    surface = build_empirical_protein_finite_state_surface(
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    )
    expected = evaluate_empirical_protein_tree_likelihood(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    ).log_likelihood
    return tree, compressed_patterns, surface, expected


def _codon_case():
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    compressed_patterns = compress_codon_site_patterns_from_records(records)
    surface = build_codon_finite_state_surface(records)
    expected = evaluate_codon_ctmc_tree_likelihood(tree, records).log_likelihood
    return tree, compressed_patterns, surface, expected


def _discrete_unordered_case():
    state_order = ["north", "south", "west"]
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    states_by_taxon = {"A": "north", "B": "south"}
    compressed_patterns = _single_site_patterns(states_by_taxon)
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.8], dtype=float)),
        state_order=state_order,
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
    )
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="unordered discrete Mk",
    )
    expected = tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    return tree, compressed_patterns, surface, expected


def _discrete_ordered_case():
    state_order = ["north", "south", "island"]
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    states_by_taxon = {"A": "north", "B": "island"}
    compressed_patterns = _single_site_patterns(states_by_taxon)
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.6], dtype=float)),
        state_order=state_order,
        model="equal-rates",
        state_ordering="ordered",
        allowed_transition_pairs={
            (0, 1),
            (1, 0),
            (1, 2),
            (2, 1),
        },
    )
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="ordered discrete Mk",
    )
    expected = tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    return tree, compressed_patterns, surface, expected


def _discrete_irreversible_case():
    state_order = ["present", "absent"]
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    states_by_taxon = {"A": "present", "B": "absent"}
    compressed_patterns = _single_site_patterns(states_by_taxon)
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.4], dtype=float)),
        state_order=state_order,
        model="irreversible",
        state_ordering="unordered",
        allowed_transition_pairs={(0, 1)},
    )
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    surface = build_discrete_mk_finite_state_surface(
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="irreversible discrete Mk",
    )
    expected = tree_log_likelihood(
        tree,
        states_by_taxon,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )
    return tree, compressed_patterns, surface, expected


def _single_site_patterns(states_by_taxon: dict[str, str]) -> CompressedAlignmentSitePatterns:
    return CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=1,
        pattern_count=1,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=(states_by_taxon["A"], states_by_taxon["B"]),
                weight=1,
                site_positions=[1],
            )
        ],
    )


def _compact_polar_rate_matrix() -> numpy.ndarray:
    rate_matrix = numpy.full((20, 20), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    state_index = _protein_state_index()
    for (left_state, right_state), rate in {
        ("A", "C"): 0.45,
        ("C", "D"): 0.35,
        ("D", "E"): 0.55,
        ("A", "E"): 0.20,
    }.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = _protein_state_index()
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _protein_state_index() -> dict[str, int]:
    return {
        state: index
        for index, state in enumerate(
            (
                "A",
                "C",
                "D",
                "E",
                "F",
                "G",
                "H",
                "I",
                "K",
                "L",
                "M",
                "N",
                "P",
                "Q",
                "R",
                "S",
                "T",
                "V",
                "W",
                "Y",
            )
        )
    }

from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    compile_codon_likelihood_analysis,
    compile_discrete_mk_likelihood_analysis,
    compile_empirical_protein_likelihood_analysis,
    compile_protein_poisson_likelihood_analysis,
    compile_selected_nucleotide_likelihood_analysis,
    evaluate_compiled_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_site_log_likelihood_rows,
    validate_compiled_finite_state_site_log_likelihood_reconstruction,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("family_name", "case_factory"),
    [
        ("nucleotide", lambda: _compiled_nucleotide_case()),
        ("protein-poisson", lambda: _compiled_protein_poisson_case()),
        ("protein-empirical", lambda: _compiled_protein_empirical_case()),
        ("codon", lambda: _compiled_codon_case()),
        ("discrete-unordered", lambda: _compiled_discrete_case()),
    ],
)
def test_compiled_site_rows_reconstruct_owned_family_totals(
    family_name: str,
    case_factory,
) -> None:
    analysis, expected_site_positions = case_factory()

    rows, row_total = evaluate_compiled_finite_state_site_log_likelihood_rows(analysis)
    direct_total = evaluate_compiled_finite_state_likelihood_analysis(analysis)

    validate_compiled_finite_state_site_log_likelihood_reconstruction(
        rows,
        expected_total_log_likelihood=row_total,
        expected_site_count=analysis.alignment_length,
        expected_pattern_count=analysis.pattern_count,
        owner_name=f"{family_name} compiled site reconstruction",
    )
    assert [row.source_site_position for row in rows] == expected_site_positions
    assert math.isclose(
        sum(row.site_weight * row.log_likelihood for row in rows),
        row_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        row_total,
        direct_total,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), family_name


def _compiled_nucleotide_case():
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_selected_nucleotide_likelihood_analysis(
        tree,
        records,
        model_name="jc69",
        owner_name="compiled site reconstruction test",
    )
    return analysis, [1, 2, 3, 4]


def _compiled_protein_poisson_case():
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_protein_poisson_likelihood_analysis(tree, records)
    return analysis, [1, 2, 3, 4]


def _compiled_protein_empirical_case():
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_empirical_protein_likelihood_analysis(
        tree,
        records,
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
    )
    return analysis, [1, 2, 3, 4]


def _compiled_codon_case():
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_codon_likelihood_analysis(tree, records)
    return analysis, [1, 2]


def _compiled_discrete_case():
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="north"),
        AlignmentRecord(identifier="B", sequence="south"),
    ]
    state_order = ["north", "south", "west"]
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
    analysis = compile_discrete_mk_likelihood_analysis(
        tree,
        records,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="unordered discrete Mk",
    )
    return analysis, [1]


def _compact_polar_rate_matrix() -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = {state: index for index, state in enumerate(state_order)}
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
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
    state_index = {state: index for index, state in enumerate(_protein_state_order())}
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _protein_state_order() -> tuple[str, ...]:
    return (
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

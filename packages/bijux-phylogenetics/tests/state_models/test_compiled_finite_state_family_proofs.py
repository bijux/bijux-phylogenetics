from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    transition_probability_matrix as discrete_transition_probability_matrix,
)
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
    evaluate_empirical_protein_tree_likelihood,
    evaluate_gtr_tree_likelihood,
    evaluate_hky85_tree_likelihood,
    evaluate_jc69_tree_likelihood,
    evaluate_protein_poisson_tree_likelihood,
    transition_probability_matrix,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import (
    build_equal_rate_codon_ctmc_rate_matrix,
)
from tests.support.two_tip_likelihood_reference import (
    two_tip_log_likelihood_for_labels,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


@pytest.mark.parametrize(
    ("family_name", "analysis_factory"),
    [
        ("nucleotide-jc69", lambda: _compiled_nucleotide_jc69_case()),
        ("nucleotide-hky85", lambda: _compiled_nucleotide_hky85_case()),
        ("nucleotide-gtr", lambda: _compiled_nucleotide_gtr_case()),
        ("protein-poisson", lambda: _compiled_protein_poisson_case()),
        ("protein-empirical", lambda: _compiled_protein_empirical_case()),
        ("codon", lambda: _compiled_codon_case()),
        ("discrete-unordered", lambda: _compiled_discrete_unordered_case()),
        ("discrete-ordered", lambda: _compiled_discrete_ordered_case()),
        ("discrete-irreversible", lambda: _compiled_discrete_irreversible_case()),
    ],
)
def test_compiled_finite_state_analysis_proves_every_goal_two_family(
    family_name: str,
    analysis_factory,
) -> None:
    analysis, expected_log_likelihood = analysis_factory()

    first_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(analysis)
    second_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(analysis)

    assert math.isclose(
        first_log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), family_name
    assert math.isclose(
        second_log_likelihood,
        expected_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    ), family_name


def _compiled_nucleotide_jc69_case():
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "jc69_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_selected_nucleotide_likelihood_analysis(
        tree,
        records,
        model_name="jc69",
        owner_name="compiled finite-state family proof",
    )
    expected = evaluate_jc69_tree_likelihood(tree, records).log_likelihood
    return analysis, expected


def _compiled_nucleotide_hky85_case():
    tree = load_tree(fixture("trees", "hky85_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "hky85_likelihood_alignment_2_taxa.fasta")
    )
    kwargs = {
        "kappa": 4.0,
        "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
    }
    analysis = compile_selected_nucleotide_likelihood_analysis(
        tree,
        records,
        model_name="hky85",
        owner_name="compiled finite-state family proof",
        **kwargs,
    )
    expected = evaluate_hky85_tree_likelihood(tree, records, **kwargs).log_likelihood
    return analysis, expected


def _compiled_nucleotide_gtr_case():
    tree = load_tree(fixture("trees", "gtr_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "gtr_likelihood_alignment_2_taxa.fasta")
    )
    kwargs = {
        "base_frequencies": {"A": 0.4, "C": 0.1, "G": 0.2, "T": 0.3},
        "exchangeabilities": {
            "AC": 1.0,
            "AG": 4.5,
            "AT": 0.8,
            "CG": 1.6,
            "CT": 2.4,
            "GT": 3.1,
        },
    }
    analysis = compile_selected_nucleotide_likelihood_analysis(
        tree,
        records,
        model_name="gtr",
        owner_name="compiled finite-state family proof",
        **kwargs,
    )
    expected = evaluate_gtr_tree_likelihood(tree, records, **kwargs).log_likelihood
    return analysis, expected


def _compiled_protein_poisson_case():
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_protein_poisson_likelihood_analysis(tree, records)
    expected = evaluate_protein_poisson_tree_likelihood(tree, records).log_likelihood
    return analysis, expected


def _compiled_protein_empirical_case():
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    rate_matrix = _compact_polar_rate_matrix()
    root_prior = _biased_root_prior()
    analysis = compile_empirical_protein_likelihood_analysis(
        tree,
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
    return analysis, expected


def _compiled_codon_case():
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_codon_likelihood_analysis(tree, records)
    state_space, rate_matrix, root_prior, _frequency_source = (
        build_equal_rate_codon_ctmc_rate_matrix()
    )
    expected = two_tip_log_likelihood_for_labels(
        observed_pairs=(("AAA", "AAG"), ("GCT", "GCC")),
        state_index=state_space.state_index,
        root_prior=root_prior,
        left_transition=transition_probability_matrix(rate_matrix, 0.1),
        right_transition=transition_probability_matrix(rate_matrix, 0.2),
    )
    return analysis, expected


def _compiled_discrete_unordered_case():
    state_order = ["north", "south", "west"]
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="north"),
        AlignmentRecord(identifier="B", sequence="south"),
    ]
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
    expected = two_tip_log_likelihood_for_labels(
        observed_pairs=(("north", "south"),),
        state_index={state: index for index, state in enumerate(state_order)},
        root_prior=root_prior,
        left_transition=discrete_transition_probability_matrix(rate_matrix, 1.0),
        right_transition=discrete_transition_probability_matrix(rate_matrix, 1.0),
    )
    return analysis, expected


def _compiled_discrete_ordered_case():
    state_order = ["north", "south", "island"]
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="north"),
        AlignmentRecord(identifier="B", sequence="island"),
    ]
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
    analysis = compile_discrete_mk_likelihood_analysis(
        tree,
        records,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="ordered discrete Mk",
    )
    expected = two_tip_log_likelihood_for_labels(
        observed_pairs=(("north", "island"),),
        state_index={state: index for index, state in enumerate(state_order)},
        root_prior=root_prior,
        left_transition=discrete_transition_probability_matrix(rate_matrix, 1.0),
        right_transition=discrete_transition_probability_matrix(rate_matrix, 1.0),
    )
    return analysis, expected


def _compiled_discrete_irreversible_case():
    state_order = ["present", "absent"]
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="present"),
        AlignmentRecord(identifier="B", sequence="absent"),
    ]
    rate_matrix = rate_matrix_from_log_parameters(
        numpy.log(numpy.array([0.4], dtype=float)),
        state_order=state_order,
        model="irreversible",
        state_ordering="unordered",
        allowed_transition_pairs={(0, 1)},
    )
    root_prior = numpy.full(len(state_order), 1.0 / len(state_order), dtype=float)
    analysis = compile_discrete_mk_likelihood_analysis(
        tree,
        records,
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        model_name="irreversible discrete Mk",
    )
    expected = two_tip_log_likelihood_for_labels(
        observed_pairs=(("present", "absent"),),
        state_index={state: index for index, state in enumerate(state_order)},
        root_prior=root_prior,
        left_transition=discrete_transition_probability_matrix(rate_matrix, 1.0),
        right_transition=discrete_transition_probability_matrix(rate_matrix, 1.0),
    )
    return analysis, expected


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

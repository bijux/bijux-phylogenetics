from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.ancestral.discrete.likelihood.likelihood_math import (
    tree_log_likelihood,
)
from bijux_phylogenetics.ancestral.discrete.likelihood.rate_matrix import (
    rate_matrix_from_log_parameters,
)
from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood.codon import (
    evaluate_codon_ctmc_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.codon_finite_state import (
    compile_codon_likelihood_analysis,
)
from bijux_phylogenetics.phylo.likelihood.codon_states import resolve_codon_state_space
from bijux_phylogenetics.phylo.likelihood.discrete_finite_state import (
    compile_discrete_mk_likelihood_analysis,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    evaluate_compiled_finite_state_likelihood_analysis,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compiled_codon_analysis_matches_existing_tree_report() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_codon_likelihood_analysis(tree, records)

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = evaluate_codon_ctmc_tree_likelihood(tree, records)

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_discrete_analysis_matches_native_kernel() -> None:
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
        site_weight_by_position={1: 6.0},
    )

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    native_log_likelihood = tree_log_likelihood(
        tree,
        {"A": "north", "B": "south"},
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )

    assert math.isclose(
        compiled_log_likelihood,
        native_log_likelihood * 6.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_discrete_analysis_accepts_variable_length_trait_labels() -> None:
    tree = load_tree(fixture("trees", "felsenstein_two_tip_tree.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="north"),
        AlignmentRecord(identifier="B", sequence="island"),
    ]
    state_order = ["north", "south", "island"]
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

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    native_log_likelihood = tree_log_likelihood(
        tree,
        {"A": "north", "B": "island"},
        state_order=state_order,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
    )

    assert math.isclose(
        compiled_log_likelihood,
        native_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_codon_analysis_respects_custom_frequency_surface() -> None:
    tree = load_tree(fixture("trees", "codon_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta")
    )
    state_space = resolve_codon_state_space()
    biased_frequencies = dict.fromkeys(state_space.state_order, 1.0)
    biased_frequencies["AAA"] = 20.0
    biased_frequencies["AAG"] = 15.0
    biased_frequencies["GCT"] = 0.5
    biased_frequencies["GCC"] = 0.5
    analysis = compile_codon_likelihood_analysis(
        tree,
        records,
        codon_frequencies=biased_frequencies,
    )

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = evaluate_codon_ctmc_tree_likelihood(
        tree,
        records,
        codon_frequencies=biased_frequencies,
    )

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

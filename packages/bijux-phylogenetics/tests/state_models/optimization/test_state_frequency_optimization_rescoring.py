from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
    FiniteStateLikelihoodSurface,
    build_protein_equal_input_finite_state_surface,
    compile_finite_state_likelihood_context,
    evaluate_codon_ctmc_tree_likelihood_from_alignment,
    evaluate_finite_state_tree_log_likelihood_from_patterns,
    optimize_codon_ctmc_state_frequencies_from_alignment,
    optimize_discrete_equal_input_state_frequencies,
    optimize_protein_equal_input_state_frequencies,
)
from bijux_phylogenetics.phylo.likelihood.discrete_observation_policies import (
    resolve_discrete_observation_leaf_vector,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.protein import (
    normalize_protein_likelihood_records,
)
from bijux_phylogenetics.phylo.likelihood.pruning import transition_probability_matrix
from bijux_phylogenetics.phylo.likelihood.state_frequency_optimization import (
    build_equal_input_rate_matrix,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_codon_state_frequency_optimization_rescores_optimized_solution() -> None:
    report = optimize_codon_ctmc_state_frequencies_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
        max_coordinate_passes=4,
    )
    optimized_frequencies = {
        row.state_label: row.optimized_value for row in report.state_rows
    }

    rescored = evaluate_codon_ctmc_tree_likelihood_from_alignment(
        fixture("trees", "codon_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "codon_likelihood_alignment_2_taxa.fasta"),
        codon_frequencies=optimized_frequencies,
    )

    assert math.isclose(
        report.optimized_log_likelihood,
        rescored.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_protein_equal_input_state_frequency_optimization_rescores_optimized_solution() -> (
    None
):
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta")
    )
    report = optimize_protein_equal_input_state_frequencies(
        tree,
        records,
        max_coordinate_passes=4,
    )
    optimized_frequencies = {
        row.state_label: row.optimized_value for row in report.state_rows
    }
    normalized_records = normalize_protein_likelihood_records(
        records,
        model_name="protein equal-input",
        gap_policy="treat-as-missing",
        missing_policy="treat-as-missing",
        ambiguity_policy="reject",
    )
    compressed_patterns = compress_alignment_site_patterns_from_records(
        normalized_records
    )
    surface = build_protein_equal_input_finite_state_surface(
        normalized_records,
        stationary_frequencies=optimized_frequencies,
        gap_policy="treat-as-missing",
        missing_policy="treat-as-missing",
        ambiguity_policy="reject",
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    rescored = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )

    assert math.isclose(
        report.optimized_log_likelihood,
        rescored,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_discrete_equal_input_state_frequency_optimization_rescores_optimized_solution() -> (
    None
):
    tree = load_tree(fixture("trees", "jc69_likelihood_tree_2_taxa.nwk"))
    records = [
        AlignmentRecord(identifier="A", sequence="0"),
        AlignmentRecord(identifier="B", sequence="2"),
    ]
    state_order = ("0", "1", "2", "3")
    report = optimize_discrete_equal_input_state_frequencies(
        tree,
        records,
        state_order=state_order,
        max_coordinate_passes=4,
    )
    optimized_frequencies = {
        row.state_label: row.optimized_value for row in report.state_rows
    }
    compressed_patterns = compress_alignment_site_patterns_from_records(
        records,
        site_tokenizer=lambda sequence: (sequence,),
    )
    rate_matrix = build_equal_input_rate_matrix(
        state_labels=state_order,
        stationary_frequencies=optimized_frequencies,
        owner_name="finite-state equal-input rescore",
    )
    stationary_vector = numpy.asarray(
        [optimized_frequencies[state] for state in state_order],
        dtype=float,
    )
    surface = FiniteStateLikelihoodSurface(
        model_name="finite-state equal-input",
        state_labels=state_order,
        root_prior=stationary_vector,
        root_prior_source="stationary",
        observation_policy="reject",
        normalization_policy=CALLER_SUPPLIED_RATE_SCALE_NORMALIZATION_POLICY,
        transition_rate_matrix=rate_matrix,
        transition_matrix_for_branch_length=lambda branch_length: (
            transition_probability_matrix(rate_matrix, branch_length)
        ),
        leaf_likelihood_vector_for_state=lambda taxon_name, observed_state: (
            resolve_discrete_observation_leaf_vector(
                observed_state,
                model_name="finite-state equal-input",
                node_name=taxon_name,
                state_order=list(state_order),
                observation_policy="reject",
            )
        ),
    )
    context = compile_finite_state_likelihood_context(tree, surface)
    rescored = evaluate_finite_state_tree_log_likelihood_from_patterns(
        context,
        compressed_patterns,
    )

    assert math.isclose(
        report.optimized_log_likelihood,
        rescored,
        rel_tol=0.0,
        abs_tol=1e-12,
    )

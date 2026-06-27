from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.io.fasta.core import load_fasta_alignment
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.alignment.partitions import (
    LocusPartition,
    LocusSegment,
)
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood,
    EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY,
    evaluate_protein_poisson_tree_likelihood,
)
from bijux_phylogenetics.phylo.likelihood.finite_state import (
    evaluate_compiled_finite_state_likelihood_analysis,
    evaluate_compiled_finite_state_site_log_likelihood_rows,
)
from bijux_phylogenetics.phylo.likelihood.protein_finite_state import (
    compile_empirical_protein_likelihood_analysis,
    compile_protein_poisson_likelihood_analysis,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_compiled_protein_poisson_analysis_matches_existing_tree_report() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta")
    )
    analysis = compile_protein_poisson_likelihood_analysis(tree, records)

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = evaluate_protein_poisson_tree_likelihood(tree, records)

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert analysis.alignment_length == 4
    assert analysis.pattern_count == 4
    assert analysis.partition_count == 1
    assert math.isclose(analysis.site_weight_total, 4.0, rel_tol=0.0, abs_tol=1e-12)
    assert analysis.context.surface.model_name == "protein Poisson"
    assert analysis.context.surface.state_count == 20
    assert analysis.context.surface.root_prior_source == "stationary"
    assert (
        analysis.context.surface.normalization_policy
        == EXPECTED_SUBSTITUTION_RATE_NORMALIZATION_POLICY
    )


def test_compiled_protein_poisson_site_rows_reconstruct_governed_report() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture(
            "alignments",
            "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta",
        )
    )
    analysis = compile_protein_poisson_likelihood_analysis(
        tree,
        records,
        ambiguity_policy="ambiguity-vector",
    )
    report = evaluate_protein_poisson_tree_likelihood(
        tree,
        records,
        ambiguity_policy="ambiguity-vector",
    )
    compiled_rows, weighted_total = evaluate_compiled_finite_state_site_log_likelihood_rows(
        analysis
    )

    assert len(compiled_rows) == report.site_count
    assert [row.source_site_position for row in compiled_rows] == [1, 2]
    assert all(row.partition_name == "full-alignment" for row in compiled_rows)
    assert all(row.partition_data_type is None for row in compiled_rows)
    assert math.isclose(
        sum(row.log_likelihood * row.site_weight for row in compiled_rows),
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        weighted_total,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_empirical_protein_analysis_matches_existing_tree_report() -> None:
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

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = evaluate_empirical_protein_tree_likelihood(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    )

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_protein_poisson_analysis_matches_ambiguity_vector_report() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture(
            "alignments",
            "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta",
        )
    )
    analysis = compile_protein_poisson_likelihood_analysis(
        tree,
        records,
        ambiguity_policy="ambiguity-vector",
    )

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = evaluate_protein_poisson_tree_likelihood(
        tree,
        records,
        ambiguity_policy="ambiguity-vector",
    )

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_empirical_protein_analysis_matches_ambiguity_vector_report() -> None:
    tree = load_tree(fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture(
            "alignments",
            "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta",
        )
    )
    rate_matrix = _compact_polar_rate_matrix()
    root_prior = _biased_root_prior()
    analysis = compile_empirical_protein_likelihood_analysis(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
        ambiguity_policy="ambiguity-vector",
    )

    compiled_log_likelihood = evaluate_compiled_finite_state_likelihood_analysis(
        analysis
    )
    report = evaluate_empirical_protein_tree_likelihood(
        tree,
        records,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
        ambiguity_policy="ambiguity-vector",
    )

    assert math.isclose(
        compiled_log_likelihood,
        report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_compiled_protein_analysis_keeps_policy_weights_and_partitions() -> None:
    tree = load_tree(fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"))
    records = load_fasta_alignment(
        fixture(
            "alignments",
            "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta",
        )
    )
    analysis = compile_protein_poisson_likelihood_analysis(
        tree,
        records,
        ambiguity_policy="ambiguity-vector",
        partition_map=(
            LocusPartition(
                name="front",
                segments=(LocusSegment(start=1, end=1),),
                total_sites=1,
                data_type="PROTEIN",
            ),
            LocusPartition(
                name="back",
                segments=(LocusSegment(start=2, end=2),),
                total_sites=1,
                data_type="PROTEIN",
            ),
        ),
        site_weight_by_position={2: 2.0},
    )

    assert analysis.partition_count == 2
    assert (
        analysis.missing_data_policy
        == "gap=treat-as-missing;missing=treat-as-missing;ambiguity=ambiguity-vector"
    )
    assert math.isclose(analysis.site_weight_total, 3.0, abs_tol=1e-12)


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

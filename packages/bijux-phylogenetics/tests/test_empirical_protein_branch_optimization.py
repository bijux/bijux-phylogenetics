from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood import (
    optimize_empirical_protein_branch_lengths_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_empirical_protein_branch_optimization_improves_bad_fixed_rate_start_tree() -> (
    None
):
    report = optimize_empirical_protein_branch_lengths_from_alignment(
        fixture("trees", "empirical_protein_branch_optimization_start_tree_2_taxa.nwk"),
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"),
        rate_matrix=_compact_polar_rate_matrix(),
        likelihood_model="fixed-rate",
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
        max_coordinate_passes=4,
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.branch_count == 2
    assert report.initial_tree_newick == "(A:1.2,B:1.4);"
    assert report.state_count == 20
    assert report.matrix_label == "compact-polar"
    assert report.root_prior_source == "provided"
    assert report.gap_policy == "treat-as-missing"
    assert report.missing_policy == "treat-as-missing"
    assert report.likelihood_model == "fixed-rate"
    assert report.alpha is None
    assert report.invariant_proportion is None
    assert math.isclose(
        report.initial_log_likelihood,
        -15.859022906023862,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.optimized_log_likelihood,
        -13.201024634746688,
        rel_tol=0.0,
        abs_tol=1e-8,
    )
    assert report.optimized_log_likelihood > report.initial_log_likelihood
    assert report.optimization_pass_count == 4
    assert report.function_evaluation_count == 401
    assert report.converged is False
    assert report.lower_branch_length_bound == 0.0
    assert report.upper_branch_length_bound == 5.0
    assert len(report.branches) == report.branch_count
    assert all(row.optimized_branch_length >= 0.0 for row in report.branches)
    assert any(
        not math.isclose(
            row.initial_branch_length,
            row.optimized_branch_length,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        for row in report.branches
    )
    assert report.branches[0].branch_id == "root:clade:A|B/taxon:A"
    assert report.branches[0].child_name == "A"
    assert report.branches[0].descendant_taxa == ["A"]
    assert math.isclose(
        report.branches[0].initial_branch_length,
        1.2,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.branches[0].optimized_branch_length,
        0.0751962334061901,
        rel_tol=0.0,
        abs_tol=1e-7,
    )
    assert report.branches[1].branch_id == "root:clade:A|B/taxon:B"
    assert report.branches[1].child_name == "B"
    assert report.branches[1].descendant_taxa == ["B"]
    assert math.isclose(
        report.branches[1].initial_branch_length,
        1.4,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.branches[1].optimized_branch_length,
        0.451186820215953,
        rel_tol=0.0,
        abs_tol=1e-7,
    )


def test_empirical_protein_branch_optimization_rejects_irrelevant_alpha_for_fixed_rate() -> (
    None
):
    with pytest.raises(ValueError, match="does not accept alpha"):
        optimize_empirical_protein_branch_lengths_from_alignment(
            fixture(
                "trees", "empirical_protein_branch_optimization_start_tree_2_taxa.nwk"
            ),
            fixture(
                "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
            ),
            rate_matrix=_compact_polar_rate_matrix(),
            likelihood_model="fixed-rate",
            alpha=0.8,
            root_prior=_biased_root_prior(),
        )


def _compact_polar_rate_matrix() -> numpy.ndarray:
    return _build_empirical_rate_matrix(
        boosted_pairs={
            ("A", "C"): 0.45,
            ("C", "D"): 0.35,
            ("D", "E"): 0.55,
            ("A", "E"): 0.20,
        }
    )


def _biased_root_prior() -> numpy.ndarray:
    prior = numpy.full(20, 0.02, dtype=float)
    state_index = _protein_state_index()
    prior[state_index["A"]] = 0.18
    prior[state_index["C"]] = 0.10
    prior[state_index["D"]] = 0.14
    prior[state_index["E"]] = 0.12
    prior[state_index["F"]] = 0.06
    return prior / float(prior.sum())


def _build_empirical_rate_matrix(
    *,
    boosted_pairs: dict[tuple[str, str], float],
) -> numpy.ndarray:
    state_order = _protein_state_order()
    state_index = _protein_state_index()
    rate_matrix = numpy.full((len(state_order), len(state_order)), 0.02, dtype=float)
    numpy.fill_diagonal(rate_matrix, 0.0)
    for (left_state, right_state), rate in boosted_pairs.items():
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        rate_matrix[left_index, right_index] = rate
        rate_matrix[right_index, left_index] = rate
    for row_index in range(rate_matrix.shape[0]):
        rate_matrix[row_index, row_index] = -float(numpy.sum(rate_matrix[row_index, :]))
    return rate_matrix


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


def _protein_state_index() -> dict[str, int]:
    return {state: index for index, state in enumerate(_protein_state_order())}

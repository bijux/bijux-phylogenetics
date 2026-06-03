from __future__ import annotations

import math
from pathlib import Path

import numpy
import pytest

from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_tree_likelihood_from_alignment,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_empirical_protein_matrix_fixed_tree_likelihood_matches_independent_fixture() -> (
    None
):
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
    )
    rate_matrix = _compact_polar_rate_matrix()
    root_prior = _biased_root_prior()

    report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        rate_matrix=rate_matrix,
        root_prior=root_prior,
        matrix_label="compact-polar",
    )

    assert report.taxa == ["A", "B"]
    assert report.site_count == 4
    assert report.pattern_count == 4
    assert report.compression_used is True
    assert report.tree_newick == "(A:0.1,B:0.2);"
    assert report.state_count == 20
    assert report.matrix_label == "compact-polar"
    assert report.root_prior_source == "provided"
    assert report.gap_policy == "treat-as-missing"
    assert report.missing_policy == "treat-as-missing"
    assert report.ambiguity_policy == "reject"
    assert math.isclose(
        report.log_likelihood,
        _expected_two_tip_log_likelihood(
            rate_matrix=rate_matrix,
            root_prior=root_prior,
            observed_pairs=(("A", "A"), ("C", "D"), ("D", "C"), ("E", "E")),
        ),
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_empirical_protein_matrix_argument_changes_likelihood() -> None:
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
    )
    compact_polar_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="compact-polar",
    )
    aromatic_shift_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        rate_matrix=_aromatic_shift_rate_matrix(),
        root_prior=_biased_root_prior(),
        matrix_label="aromatic-shift",
    )

    assert compact_polar_report.matrix_label == "compact-polar"
    assert aromatic_shift_report.matrix_label == "aromatic-shift"
    assert not math.isclose(
        compact_polar_report.log_likelihood,
        aromatic_shift_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_empirical_protein_matrix_validation_rejects_non_20x20_input() -> None:
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
    )

    with pytest.raises(
        InvalidAlignmentError, match="requires one 20x20 empirical protein rate matrix"
    ):
        evaluate_empirical_protein_tree_likelihood_from_alignment(
            tree_path,
            alignment_path,
            rate_matrix=numpy.eye(19, dtype=float),
        )


def test_empirical_protein_matrix_validation_rejects_nonzero_row_sums() -> None:
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"
    )
    invalid_matrix = _compact_polar_rate_matrix()
    invalid_matrix[0, 0] += 0.25

    with pytest.raises(InvalidAlignmentError, match="rows must sum to zero"):
        evaluate_empirical_protein_tree_likelihood_from_alignment(
            tree_path,
            alignment_path,
            rate_matrix=invalid_matrix,
        )


def test_empirical_protein_matrix_ambiguity_policies_change_likelihood() -> None:
    tree_path = fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk")
    alignment_path = fixture(
        "alignments", "protein_poisson_likelihood_alignment_ambiguity_2_taxa.fasta"
    )

    missing_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        ambiguity_policy="treat-as-missing",
    )
    ambiguity_report = evaluate_empirical_protein_tree_likelihood_from_alignment(
        tree_path,
        alignment_path,
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        ambiguity_policy="ambiguity-vector",
    )

    assert missing_report.ambiguity_policy == "treat-as-missing"
    assert ambiguity_report.ambiguity_policy == "ambiguity-vector"
    assert math.isfinite(missing_report.log_likelihood)
    assert math.isfinite(ambiguity_report.log_likelihood)
    assert not math.isclose(
        missing_report.log_likelihood,
        ambiguity_report.log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
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


def _aromatic_shift_rate_matrix() -> numpy.ndarray:
    return _build_empirical_rate_matrix(
        boosted_pairs={
            ("A", "D"): 0.50,
            ("C", "E"): 0.40,
            ("D", "F"): 0.45,
            ("E", "Y"): 0.30,
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


def _expected_two_tip_log_likelihood(
    *,
    rate_matrix: numpy.ndarray,
    root_prior: numpy.ndarray,
    observed_pairs: tuple[tuple[str, str], ...],
) -> float:
    left_transition = _transition_probability_matrix(rate_matrix, 0.1)
    right_transition = _transition_probability_matrix(rate_matrix, 0.2)
    state_index = _protein_state_index()
    probability = 1.0
    for left_state, right_state in observed_pairs:
        left_index = state_index[left_state]
        right_index = state_index[right_state]
        site_probability = 0.0
        for root_index, root_probability in enumerate(root_prior):
            site_probability += (
                float(root_probability)
                * float(left_transition[root_index, left_index])
                * float(right_transition[root_index, right_index])
            )
        probability *= site_probability
    return math.log(probability)


def _transition_probability_matrix(
    rate_matrix: numpy.ndarray,
    branch_length: float,
) -> numpy.ndarray:
    eigenvalues, eigenvectors = numpy.linalg.eig(rate_matrix)
    inverse_vectors = numpy.linalg.inv(eigenvectors)
    diagonal = numpy.diag(numpy.exp(eigenvalues * branch_length))
    transition = eigenvectors @ diagonal @ inverse_vectors
    transition = numpy.real_if_close(transition, tol=1000).astype(float)
    transition[transition < 0.0] = 0.0
    row_sums = transition.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0.0] = 1.0
    return transition / row_sums


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

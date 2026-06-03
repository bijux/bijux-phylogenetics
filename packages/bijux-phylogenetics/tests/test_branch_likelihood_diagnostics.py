from __future__ import annotations

import math
from pathlib import Path

import numpy

from bijux_phylogenetics.phylo import likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_empirical_protein_branch_likelihood_diagnostics_from_alignment,
    evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment,
    evaluate_protein_poisson_branch_likelihood_diagnostics_from_alignment,
)

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_branch_likelihood_diagnostics() -> None:
    assert (
        likelihood_api.evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment
        is evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment
    )
    assert (
        likelihood_api.evaluate_protein_poisson_branch_likelihood_diagnostics_from_alignment
        is evaluate_protein_poisson_branch_likelihood_diagnostics_from_alignment
    )
    assert (
        likelihood_api.evaluate_empirical_protein_branch_likelihood_diagnostics_from_alignment
        is evaluate_empirical_protein_branch_likelihood_diagnostics_from_alignment
    )


def test_nucleotide_branch_likelihood_diagnostics_report_collapse_delta_per_branch() -> (
    None
):
    report = evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment(
        fixture("trees", "strict_clock_nonclock_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
        model_name="jc69",
    )

    assert report.model_name == "JC69"
    assert report.taxa == ["A", "B", "C", "D"]
    assert report.site_count == 12
    assert report.pattern_count == 2
    assert report.branch_count == 6
    assert len(report.branch_diagnostics) == 6
    assert math.isfinite(report.baseline_log_likelihood)
    for row in report.branch_diagnostics:
        assert row.branch_id
        assert math.isfinite(row.branch_length)
        assert math.isfinite(row.contribution_proxy)
        assert math.isclose(
            row.contribution_proxy,
            report.baseline_log_likelihood - row.collapsed_branch_log_likelihood,
            rel_tol=0.0,
            abs_tol=1e-12,
        )


def test_nucleotide_branch_likelihood_diagnostics_flag_zero_length_branches(
    tmp_path: Path,
) -> None:
    tree_path = tmp_path / "zero-branch-tree.nwk"
    tree_path.write_text(
        "(((A:0.0,B:0.5):0.4,C:0.6):0.3,D:0.8);",
        encoding="utf-8",
    )

    report = evaluate_nucleotide_branch_likelihood_diagnostics_from_alignment(
        tree_path,
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
        model_name="jc69",
    )

    row_by_descendants = {
        tuple(row.descendant_taxa): row for row in report.branch_diagnostics
    }
    assert "zero-branch-length" in row_by_descendants[("A",)].warning_flags


def test_protein_poisson_branch_likelihood_diagnostics_report_rows_are_finite() -> None:
    report = evaluate_protein_poisson_branch_likelihood_diagnostics_from_alignment(
        fixture("trees", "protein_poisson_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "protein_poisson_likelihood_alignment_2_taxa.fasta"),
    )

    assert report.model_name == "protein Poisson"
    assert report.branch_count == 2
    assert len(report.branch_diagnostics) == 2
    assert all(
        math.isfinite(row.contribution_proxy) for row in report.branch_diagnostics
    )


def test_empirical_protein_branch_likelihood_diagnostics_support_fixed_rate_surface() -> (
    None
):
    report = evaluate_empirical_protein_branch_likelihood_diagnostics_from_alignment(
        fixture("trees", "empirical_protein_likelihood_tree_2_taxa.nwk"),
        fixture("alignments", "empirical_protein_likelihood_alignment_2_taxa.fasta"),
        rate_matrix=_compact_polar_rate_matrix(),
        root_prior=_biased_root_prior(),
        likelihood_model="fixed-rate",
        matrix_label="compact-polar",
    )

    assert report.model_name == "empirical protein matrix"
    assert report.branch_count == 2
    assert len(report.branch_diagnostics) == 2
    assert all(
        math.isfinite(row.contribution_proxy) for row in report.branch_diagnostics
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

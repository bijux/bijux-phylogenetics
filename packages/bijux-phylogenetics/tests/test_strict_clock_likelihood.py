from __future__ import annotations

import math
from pathlib import Path

import pytest

import bijux_phylogenetics.phylo.likelihood as likelihood_api
from bijux_phylogenetics.phylo.likelihood import (
    evaluate_jc69_tree_likelihood_from_alignment,
    fit_strict_clock_likelihood_from_alignment,
)
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(group: str, name: str) -> Path:
    return FIXTURES / group / name


def test_package_likelihood_gateway_exports_strict_clock_surface() -> None:
    assert (
        likelihood_api.fit_strict_clock_likelihood_from_alignment
        is fit_strict_clock_likelihood_from_alignment
    )


def test_strict_clock_likelihood_fits_one_global_rate_on_time_tree_fixture() -> None:
    report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )

    assert report.model_name == "JC69"
    assert report.taxa == ["A", "B", "C", "D"]
    assert report.site_count == 12
    assert report.pattern_count == 2
    assert report.branch_count == 6
    assert report.compression_used is True
    assert report.time_tree_newick == "(((A:1,B:1):1,C:2):1,D:3);"
    assert report.scaled_tree_newick == (
        "(((A:0.521738366943483,B:0.521738366943483):0.521738366943483,"
        "C:1.04347673388697):0.521738366943483,D:1.56521510083045);"
    )
    assert math.isclose(
        report.initial_clock_rate,
        1.0,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.optimized_clock_rate,
        0.5217383669434831,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.initial_log_likelihood,
        -65.42488338407935,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert math.isclose(
        report.optimized_log_likelihood,
        -64.38812241070909,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.parameter_count == 1
    assert math.isclose(
        report.aic,
        130.77624482141817,
        rel_tol=0.0,
        abs_tol=1e-12,
    )
    assert report.function_evaluation_count == 51
    assert report.converged is True
    assert report.lower_clock_rate_bound == 1e-6
    assert report.upper_clock_rate_bound == 5.0

    rows_by_descendant_taxa = {
        tuple(row.descendant_taxa): row for row in report.branch_rows
    }
    expected_rows = {
        ("A",): (1.0, 0.5217383669434831),
        ("B",): (1.0, 0.5217383669434831),
        ("A", "B"): (1.0, 0.5217383669434831),
        ("C",): (2.0, 1.0434767338869662),
        ("A", "B", "C"): (1.0, 0.5217383669434831),
        ("D",): (3.0, 1.5652151008304491),
    }
    assert set(rows_by_descendant_taxa) == set(expected_rows)
    for descendant_taxa, (
        expected_time_duration,
        expected_branch_length,
    ) in expected_rows.items():
        row = rows_by_descendant_taxa[descendant_taxa]
        assert math.isclose(
            row.time_duration,
            expected_time_duration,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        assert math.isclose(
            row.optimized_branch_length,
            expected_branch_length,
            rel_tol=0.0,
            abs_tol=1e-12,
        )
        assert math.isclose(
            row.optimized_clock_rate,
            report.optimized_clock_rate,
            rel_tol=0.0,
            abs_tol=1e-12,
        )


def test_strict_clock_likelihood_differs_from_nonclock_likelihood_on_nonclock_fixture() -> (
    None
):
    strict_clock_report = fit_strict_clock_likelihood_from_alignment(
        fixture("trees", "strict_clock_time_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )
    free_branch_report = evaluate_jc69_tree_likelihood_from_alignment(
        fixture("trees", "strict_clock_nonclock_tree_4_taxa.nwk"),
        fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
    )

    assert (
        free_branch_report.log_likelihood > strict_clock_report.optimized_log_likelihood
    )
    assert not math.isclose(
        free_branch_report.log_likelihood,
        strict_clock_report.optimized_log_likelihood,
        rel_tol=0.0,
        abs_tol=1e-12,
    )


def test_strict_clock_likelihood_requires_explicit_branch_lengths() -> None:
    with pytest.raises(InvalidBranchLengthError) as error:
        fit_strict_clock_likelihood_from_alignment(
            fixture("trees", "example_tree_no_lengths.nwk"),
            fixture("alignments", "strict_clock_likelihood_alignment_4_taxa.fasta"),
        )

    assert "requires explicit branch lengths" in str(error.value)

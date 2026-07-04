from __future__ import annotations

import math
from pathlib import Path

import pytest

from bijux_phylogenetics.comparative.covariance import (
    summarize_brownian_covariance,
    summarize_brownian_covariance_from_tree,
    write_brownian_covariance_long_table,
    write_brownian_covariance_matrix_table,
)
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import ComparativeMethodError

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_summarize_brownian_covariance_matches_rooted_ultrametric_matrix() -> None:
    report = summarize_brownian_covariance(fixture("example_tree.nwk"))

    assert report.taxa == ["A", "B", "C", "D"]
    assert report.tree_is_rooted is True
    assert report.tree_is_ultrametric is True
    assert report.matrix_dimension == 4
    assert report.matrix_rank == 4
    assert report.singular is False
    assert report.near_singular is False
    matrix = {
        (row.left_taxon, row.right_taxon): row.shared_ancestry_covariance
        for row in report.rows
    }
    assert math.isclose(matrix[("A", "A")], 0.3, abs_tol=1e-12)
    assert math.isclose(matrix[("A", "B")], 0.2, abs_tol=1e-12)
    assert math.isclose(matrix[("C", "D")], 0.1, abs_tol=1e-12)
    assert math.isclose(matrix[("B", "D")], 0.0, abs_tol=1e-12)


def test_summarize_brownian_covariance_from_tree_matches_path_surface() -> None:
    tree = load_tree(fixture("example_tree.nwk"))

    direct_report = summarize_brownian_covariance_from_tree(
        tree,
        tree_path=fixture("example_tree.nwk"),
    )
    path_report = summarize_brownian_covariance(fixture("example_tree.nwk"))

    assert direct_report.tree_path == fixture("example_tree.nwk")
    assert direct_report.taxa == path_report.taxa
    assert direct_report.tree_is_rooted == path_report.tree_is_rooted
    assert direct_report.tree_is_ultrametric == path_report.tree_is_ultrametric
    assert direct_report.matrix_rank == path_report.matrix_rank
    assert direct_report.condition_number == path_report.condition_number
    assert direct_report.rows == path_report.rows


def test_summarize_brownian_covariance_handles_non_ultrametric_tip_order() -> None:
    report = summarize_brownian_covariance(
        fixture("example_tree_ladderized.nwk"),
        taxa=["D", "C", "B", "A"],
    )

    assert report.taxa == ["D", "C", "B", "A"]
    assert report.tree_is_rooted is True
    assert report.tree_is_ultrametric is False
    assert math.isclose(report.minimum_root_to_tip_depth, 0.1, abs_tol=1e-12)
    assert math.isclose(report.maximum_root_to_tip_depth, 0.3, abs_tol=1e-12)
    matrix = {
        (row.left_taxon, row.right_taxon): row.shared_ancestry_covariance
        for row in report.rows
    }
    assert math.isclose(matrix[("D", "D")], 0.1, abs_tol=1e-12)
    assert math.isclose(matrix[("C", "B")], 0.1, abs_tol=1e-12)
    assert math.isclose(matrix[("A", "D")], 0.0, abs_tol=1e-12)
    assert math.isclose(report.condition_number, 5.56155281280883, abs_tol=1e-12)


def test_summarize_brownian_covariance_handles_unrooted_branch_length_tree() -> None:
    report = summarize_brownian_covariance(fixture("example_tree_unrooted.nwk"))

    assert report.tree_is_rooted is False
    assert report.tree_is_ultrametric is False
    assert report.matrix_dimension == 4
    matrix = {
        (row.left_taxon, row.right_taxon): row.shared_ancestry_covariance
        for row in report.rows
    }
    assert math.isclose(matrix[("A", "A")], 0.1, abs_tol=1e-12)
    assert math.isclose(matrix[("B", "B")], 0.2, abs_tol=1e-12)
    assert math.isclose(matrix[("A", "B")], 0.0, abs_tol=1e-12)
    assert math.isclose(matrix[("C", "D")], 0.0, abs_tol=1e-12)


def test_summarize_brownian_covariance_reports_singular_zero_length_matrix() -> None:
    report = summarize_brownian_covariance(fixture("example_tree_zero_lengths.nwk"))

    assert report.singular is True
    assert report.near_singular is True
    assert report.positive_definite is False
    assert report.raw_log_determinant is None
    assert math.isinf(report.condition_number)


def test_summarize_brownian_covariance_rejects_missing_branch_lengths() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        summarize_brownian_covariance(
            fixture("example_tree_branch_lengths_missing.nwk")
        )

    assert error.value.details["failure_reason"] == (
        "brownian_covariance_branch_lengths_incomplete"
    )


def test_summarize_brownian_covariance_rejects_negative_branch_lengths() -> None:
    with pytest.raises(ComparativeMethodError) as error:
        summarize_brownian_covariance(fixture("example_tree_negative_length.nwk"))

    assert error.value.details["failure_reason"] == (
        "brownian_covariance_negative_branch_lengths"
    )


def test_write_brownian_covariance_tables_preserve_explicit_tip_order(
    tmp_path: Path,
) -> None:
    report = summarize_brownian_covariance(
        fixture("example_tree.nwk"),
        taxa=["D", "C", "B", "A"],
    )
    long_path = tmp_path / "brownian-covariance-long.tsv"
    matrix_path = tmp_path / "brownian-covariance-matrix.tsv"

    write_brownian_covariance_long_table(long_path, report)
    write_brownian_covariance_matrix_table(matrix_path, report)

    long_rows = long_path.read_text(encoding="utf-8").splitlines()
    matrix_rows = matrix_path.read_text(encoding="utf-8").splitlines()
    assert long_rows[0].startswith(
        "left_taxon\tright_taxon\tis_diagonal\tshared_ancestry_covariance"
    )
    assert long_rows[1].startswith("D\tD\ttrue\t0.3")
    assert matrix_rows[0] == "taxon\tD\tC\tB\tA"
    assert matrix_rows[1] == "D\t0.3\t0.1\t0\t0"

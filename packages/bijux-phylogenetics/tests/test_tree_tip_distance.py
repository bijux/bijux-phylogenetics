from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.phylo.topology.tip_distances import (
    compute_tree_tip_distance_matrix,
    write_tree_tip_distance_long_table,
    write_tree_tip_distance_matrix,
)
from bijux_phylogenetics.runtime.errors import InvalidBranchLengthError


def fixture(name: str) -> Path:
    return Path(__file__).resolve().parent / "fixtures" / "trees" / name


def test_compute_tree_tip_distance_matrix_reports_rooted_branch_lengths() -> None:
    report = compute_tree_tip_distance_matrix(fixture("example_tree.nwk"))

    assert report.identifiers == ["A", "B", "C", "D"]
    assert report.rooted is True
    assert report.complete_branch_lengths is True
    assert report.branch_length_count == 6
    assert report.expected_branch_length_count == 6
    assert report.pair_count == 16
    assert report.diagonal_zero is True
    assert report.symmetric is True
    assert report.matrix == [
        [0.0, 0.2, 0.6, 0.6],
        [0.2, 0.0, 0.6, 0.6],
        [0.6, 0.6, 0.0, 0.4],
        [0.6, 0.6, 0.4, 0.0],
    ]


def test_compute_tree_tip_distance_matrix_reports_unrooted_branch_lengths() -> None:
    report = compute_tree_tip_distance_matrix(fixture("example_tree_unrooted.nwk"))

    assert report.identifiers == ["A", "B", "C", "D"]
    assert report.rooted is False
    assert report.complete_branch_lengths is True
    assert report.matrix == [
        [0.0, 0.3, 0.4, 0.5],
        [0.3, 0.0, 0.5, 0.6],
        [0.4, 0.5, 0.0, 0.7],
        [0.5, 0.6, 0.7, 0.0],
    ]


def test_compute_tree_tip_distance_matrix_rejects_missing_branch_lengths_by_default() -> (
    None
):
    with pytest.raises(
        InvalidBranchLengthError,
        match="tree requires complete branch lengths for tip-distance calculations",
    ):
        compute_tree_tip_distance_matrix(fixture("example_tree_partial_lengths.nwk"))


def test_compute_tree_tip_distance_matrix_can_use_explicit_unit_length_policy() -> None:
    report = compute_tree_tip_distance_matrix(
        fixture("example_tree_partial_lengths.nwk"),
        missing_branch_length_policy="unit-length",
    )

    assert report.complete_branch_lengths is False
    assert report.branch_length_count == 3
    assert report.expected_branch_length_count == 4
    assert report.matrix == [
        [0.0, 1.1, 0.8],
        [1.1, 0.0, 1.7],
        [0.8, 1.7, 0.0],
    ]


def test_tree_tip_distance_writers_emit_matrix_and_long_form_outputs(
    tmp_path: Path,
) -> None:
    report = compute_tree_tip_distance_matrix(fixture("example_tree.nwk"))
    matrix_path = tmp_path / "cophenetic-matrix.tsv"
    long_path = tmp_path / "cophenetic-long.tsv"

    write_tree_tip_distance_matrix(matrix_path, report)
    write_tree_tip_distance_long_table(long_path, report)

    assert matrix_path.read_text(encoding="utf-8").splitlines() == [
        "taxon\tA\tB\tC\tD",
        "A\t0\t0.2\t0.6\t0.6",
        "B\t0.2\t0\t0.6\t0.6",
        "C\t0.6\t0.6\t0\t0.4",
        "D\t0.6\t0.6\t0.4\t0",
    ]
    assert long_path.read_text(encoding="utf-8").splitlines()[:5] == [
        "left_identifier\tright_identifier\tdistance",
        "A\tA\t0",
        "A\tB\t0.2",
        "A\tC\t0.6",
        "A\tD\t0.6",
    ]

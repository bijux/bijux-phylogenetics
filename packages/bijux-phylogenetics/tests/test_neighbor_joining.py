from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_distance_matrix_fixture,
)
from bijux_phylogenetics.distance import build_tree_from_imported_distance_matrix
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError


def _write_distance_matrix(path: Path, rows: list[str]) -> Path:
    path.write_text(
        "left_identifier\tright_identifier\tdistance\tcomparable_sites\n"
        + "\n".join(rows)
        + "\n",
        encoding="utf-8",
    )
    return path


def test_build_tree_from_imported_distance_matrix_matches_analytical_nj_fixture() -> (
    None
):
    tree, report = build_tree_from_imported_distance_matrix(
        get_shared_distance_matrix_fixture("analytical_three_taxon").path,
        method="neighbor-joining",
    )
    assert dumps_newick(tree) == "(A:0,B:0.125,C:0.5)Inner1;"
    assert tree.rooted is False
    assert report.method == "neighbor-joining"
    assert report.taxon_count == 3


def test_build_tree_from_imported_distance_matrix_handles_tied_joins_deterministically() -> (
    None
):
    tree, _report = build_tree_from_imported_distance_matrix(
        get_shared_distance_matrix_fixture("ultrametric_four_taxon").path,
        method="neighbor-joining",
    )
    assert dumps_newick(tree) == "((A:1,B:1)Inner1:4,C:1,D:1)Inner2;"


def test_build_tree_from_imported_distance_matrix_matches_nonultrametric_branch_lengths() -> (
    None
):
    tree, _report = build_tree_from_imported_distance_matrix(
        get_shared_distance_matrix_fixture("nonultrametric_four_taxon").path,
        method="neighbor-joining",
    )
    assert dumps_newick(tree) == "((A:1,B:1)Inner1:4.5,C:0.5,D:1.5)Inner2;"


def test_build_tree_from_imported_distance_matrix_rejects_nonzero_diagonal(
    tmp_path: Path,
) -> None:
    matrix_path = _write_distance_matrix(
        tmp_path / "distance-matrix.tsv",
        [
            "A\tA\t0.1\t10",
            "A\tB\t1\t10",
            "B\tA\t1\t10",
            "B\tB\t0\t10",
        ],
    )
    try:
        build_tree_from_imported_distance_matrix(
            matrix_path,
            method="neighbor-joining",
        )
    except InvalidDistanceMatrixError as error:
        assert error.code == "invalid_distance_matrix_error"
        assert "nonzero diagonal" in str(error)
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidDistanceMatrixError")


def test_build_tree_from_imported_distance_matrix_rejects_negative_distance(
    tmp_path: Path,
) -> None:
    matrix_path = _write_distance_matrix(
        tmp_path / "distance-matrix.tsv",
        [
            "A\tA\t0\t10",
            "A\tB\t-1\t10",
            "B\tA\t-1\t10",
            "B\tB\t0\t10",
        ],
    )
    try:
        build_tree_from_imported_distance_matrix(
            matrix_path,
            method="neighbor-joining",
        )
    except InvalidDistanceMatrixError as error:
        assert error.code == "invalid_distance_matrix_error"
        assert "negative distances" in str(error)
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidDistanceMatrixError")

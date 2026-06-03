from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    fit_ordinary_least_squares_tree,
    fit_ordinary_least_squares_tree_from_imported_distance_matrix,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree

FIXTURES = Path(__file__).parent / "fixtures"


def fixture(name: str) -> Path:
    for group in ("metadata", "trees"):
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def distance_lookup(name: str) -> tuple[list[str], dict[tuple[str, str], float]]:
    matrix_path = fixture(name)
    validation = validate_imported_distance_matrix(matrix_path)
    entries = load_imported_distance_matrix(matrix_path)
    lookup = {(identifier, identifier): 0.0 for identifier in validation.identifiers}
    for entry in entries:
        lookup[(entry.left_identifier, entry.right_identifier)] = entry.distance
    return validation.identifiers, lookup


def test_fit_ordinary_least_squares_tree_reports_full_matrix_diagnostics() -> None:
    tree = load_tree(fixture("example_tree_minimum_evolution_five_taxon.nwk"))
    identifiers, lookup = distance_lookup(
        "example_distance_matrix_minimum_evolution_five_taxon.tsv"
    )

    fitted_tree, report = fit_ordinary_least_squares_tree(
        tree,
        identifiers,
        lookup,
    )

    assert dumps_newick(fitted_tree) == "((A:1,B:2):3,C:4,(D:5,E:6):7);"
    assert report.residual_sum_squares == 0.0
    assert report.matrix_rank == 7
    assert report.condition_number == 4.08308918215
    assert report.negative_branch_count == 0
    assert report.fitted_distance_matrix == [
        [0.0, 3.0, 8.0, 16.0, 17.0],
        [3.0, 0.0, 9.0, 17.0, 18.0],
        [8.0, 9.0, 0.0, 16.0, 17.0],
        [16.0, 17.0, 16.0, 0.0, 11.0],
        [17.0, 18.0, 17.0, 11.0, 0.0],
    ]
    assert report.residual_matrix == [
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -0.0, -0.0],
        [0.0, 0.0, -0.0, 0.0, 0.0],
        [0.0, 0.0, -0.0, 0.0, 0.0],
    ]
    assert [row.fitted_branch_length for row in report.branch_fits] == [
        3.0,
        1.0,
        2.0,
        4.0,
        7.0,
        5.0,
        6.0,
    ]


def test_fit_ordinary_least_squares_tree_reports_negative_branches() -> None:
    report = fit_ordinary_least_squares_tree_from_imported_distance_matrix(
        fixture(
            "example_distance_matrix_ordinary_least_squares_negative_branch_five_taxon.tsv"
        ),
        fixture("example_tree_minimum_evolution_five_taxon.nwk"),
    )[1]

    assert report.residual_sum_squares == 44.333333333333
    assert report.matrix_rank == 7
    assert report.condition_number == 4.08308918215
    assert report.negative_branch_count == 1
    assert [row.fitted_branch_length for row in report.branch_fits] == [
        -2.75,
        2.0,
        5.0,
        6.75,
        8.25,
        2.666666666667,
        7.333333333333,
    ]

from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import (
    fit_nonnegative_least_squares_tree,
    fit_nonnegative_least_squares_tree_from_imported_distance_matrix,
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


def test_fit_nonnegative_least_squares_tree_preserves_positive_analytical_fixture() -> (
    None
):
    tree = load_tree(fixture("example_tree_minimum_evolution_five_taxon.nwk"))
    identifiers, lookup = distance_lookup(
        "example_distance_matrix_minimum_evolution_five_taxon.tsv"
    )

    fitted_tree, report = fit_nonnegative_least_squares_tree(
        tree,
        identifiers,
        lookup,
    )

    assert dumps_newick(fitted_tree) == "((A:1,B:2):3,C:4,(D:5,E:6):7);"
    assert report.residual_sum_squares == 0.0
    assert report.active_constraint_count == 0
    assert report.active_constraints == []
    assert [row.fitted_branch_length for row in report.branch_fits] == [
        3.0,
        1.0,
        2.0,
        4.0,
        7.0,
        5.0,
        6.0,
    ]


def test_fit_nonnegative_least_squares_tree_constrains_negative_ols_branch() -> None:
    fitted_tree, report = (
        fit_nonnegative_least_squares_tree_from_imported_distance_matrix(
            fixture(
                "example_distance_matrix_ordinary_least_squares_negative_branch_five_taxon.tsv"
            ),
            fixture("example_tree_minimum_evolution_five_taxon.nwk"),
        )
    )

    assert dumps_newick(fitted_tree) == (
        "((A:0.777777777778,B:3.777777777778):0,C:6.444444444444,(D:2.666666666667,E:7.333333333333):7.333333333333);"
    )
    assert report.residual_sum_squares == 57.777777777778
    assert report.condition_number == 4.08308918215
    assert report.active_constraint_count == 1
    assert report.active_constraints[0].descendant_taxa == ["A", "B"]
    assert [row.fitted_branch_length for row in report.branch_fits] == [
        0.0,
        0.777777777778,
        3.777777777778,
        6.444444444444,
        7.333333333333,
        2.666666666667,
        7.333333333333,
    ]

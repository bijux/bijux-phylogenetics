from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    fit_minimum_evolution_tree,
    fit_minimum_evolution_tree_from_imported_distance_matrix,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

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


def test_fit_minimum_evolution_tree_recovers_analytical_branch_lengths() -> None:
    tree_path = fixture("example_tree_minimum_evolution_five_taxon.nwk")
    identifiers, lookup = distance_lookup(
        "example_distance_matrix_minimum_evolution_five_taxon.tsv"
    )
    tree = load_tree(tree_path)

    fitted_tree, report = fit_minimum_evolution_tree(
        tree,
        identifiers,
        lookup,
    )

    assert tree.total_branch_length() == 63.0
    assert dumps_newick(fitted_tree) == "((A:1,B:2):3,C:4,(D:5,E:6):7);"
    assert report.minimum_evolution_score == 28.0
    assert report.total_fitted_branch_length == 28.0
    assert report.branch_count == 7
    assert report.pair_count == 10
    assert report.negative_branch_count == 0
    assert [row.fitted_branch_length for row in report.branch_fits] == [
        3.0,
        1.0,
        2.0,
        4.0,
        7.0,
        5.0,
        6.0,
    ]


def test_fit_minimum_evolution_tree_from_imported_distance_matrix_supports_public_surface() -> (
    None
):
    fitted_tree, report = fit_minimum_evolution_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv"),
        fixture("example_tree_minimum_evolution_five_taxon.nwk"),
    )

    assert dumps_newick(fitted_tree) == "((A:1,B:2):3,C:4,(D:5,E:6):7);"
    assert report.taxa == ["A", "B", "C", "D", "E"]
    assert report.minimum_evolution_score == 28.0


def test_fit_minimum_evolution_tree_rejects_taxon_mismatch(tmp_path: Path) -> None:
    tree_path = tmp_path / "mismatch-tree.nwk"
    tree_path.write_text("((A:1,B:1):1,C:1,(D:1,F:1):1);\n", encoding="utf-8")

    with pytest.raises(InvalidDistanceMatrixError, match="taxa do not match"):
        fit_minimum_evolution_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_minimum_evolution_five_taxon.tsv"),
            tree_path,
        )

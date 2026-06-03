from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    fit_fitch_margoliash_tree,
    fit_fitch_margoliash_tree_from_imported_distance_matrix,
    fit_minimum_evolution_tree_from_imported_distance_matrix,
    load_imported_distance_matrix,
    validate_imported_distance_matrix,
)
from bijux_phylogenetics.io.newick import dumps_newick, write_newick
from bijux_phylogenetics.io.trees import load_tree
from bijux_phylogenetics.phylo.topology import compute_tree_tip_distance_matrix
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


def weighted_residual_sum_squares(
    tree_path: Path,
    matrix_name: str,
    *,
    weighting_power: float,
) -> float:
    identifiers, lookup = distance_lookup(matrix_name)
    tip_distances = compute_tree_tip_distance_matrix(tree_path)
    predicted = {
        (row.left_identifier, row.right_identifier): row.distance
        for row in tip_distances.pairs
    }
    total = 0.0
    for left_index, left_identifier in enumerate(identifiers):
        for right_identifier in identifiers[left_index + 1 :]:
            observed = lookup[(left_identifier, right_identifier)]
            residual = observed - predicted[(left_identifier, right_identifier)]
            total += (observed ** (-weighting_power)) * (residual**2)
    return round(total, 12)


def test_fit_fitch_margoliash_tree_matches_weighted_fixture() -> None:
    tree = load_tree(fixture("example_tree_minimum_evolution_five_taxon.nwk"))
    identifiers, lookup = distance_lookup(
        "example_distance_matrix_fitch_margoliash_five_taxon.tsv"
    )

    fitted_tree, report = fit_fitch_margoliash_tree(
        tree,
        identifiers,
        lookup,
    )

    assert dumps_newick(fitted_tree) == (
        "((A:1.186387219954,B:0.813612780046):3.461532690732,C:1.538467309268,"
        "(D:4.498139020967,E:8.501860979033):6.564266585263);"
    )
    assert report.weighting_power == 2.0
    assert report.residual_sum_squares == 4.562310999055
    assert report.weighted_residual_sum_squares == 0.020472720566
    assert report.matrix_rank == 7
    assert report.negative_branch_count == 0
    assert [row.fitted_branch_length for row in report.branch_fits] == [
        3.461532690732,
        1.186387219954,
        0.813612780046,
        1.538467309268,
        6.564266585263,
        4.498139020967,
        8.501860979033,
    ]


def test_fit_fitch_margoliash_tree_beats_ols_under_weighted_objective(
    tmp_path: Path,
) -> None:
    matrix_name = "example_distance_matrix_fitch_margoliash_five_taxon.tsv"
    matrix_path = fixture(matrix_name)
    tree_path = fixture("example_tree_minimum_evolution_five_taxon.nwk")

    fm_tree, fm_report = fit_fitch_margoliash_tree_from_imported_distance_matrix(
        matrix_path,
        tree_path,
    )
    ols_tree, _ols_report = fit_minimum_evolution_tree_from_imported_distance_matrix(
        matrix_path,
        tree_path,
    )
    fm_tree_path = write_newick(tmp_path / "fm-tree.nwk", fm_tree)
    ols_tree_path = write_newick(tmp_path / "ols-tree.nwk", ols_tree)

    assert dumps_newick(fm_tree) != dumps_newick(ols_tree)
    assert (
        weighted_residual_sum_squares(
            fm_tree_path,
            matrix_name,
            weighting_power=2.0,
        )
        == fm_report.weighted_residual_sum_squares
    )
    assert (
        weighted_residual_sum_squares(
            ols_tree_path,
            matrix_name,
            weighting_power=2.0,
        )
        == 0.027825546887
    )
    assert fm_report.weighted_residual_sum_squares < weighted_residual_sum_squares(
        ols_tree_path,
        matrix_name,
        weighting_power=2.0,
    )


def test_fit_fitch_margoliash_tree_rejects_zero_off_diagonal_distance(
    tmp_path: Path,
) -> None:
    matrix_path = tmp_path / "zero-distance.tsv"
    matrix_path.write_text(
        "left_identifier\tright_identifier\tdistance\tcomparable_sites\n"
        "A\tA\t0\t10\n"
        "A\tB\t0\t10\n"
        "B\tA\t0\t10\n"
        "B\tB\t0\t10\n",
        encoding="utf-8",
    )
    tree_path = tmp_path / "two-tip-tree.nwk"
    tree_path.write_text("(A:1,B:1);\n", encoding="utf-8")

    with pytest.raises(
        InvalidDistanceMatrixError,
        match="strictly positive finite off-diagonal observed distances",
    ):
        fit_fitch_margoliash_tree_from_imported_distance_matrix(
            matrix_path,
            tree_path,
        )

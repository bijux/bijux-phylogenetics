from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    build_tree_from_imported_distance_matrix,
    build_upgma_tree,
    build_wpgma_tree,
    load_imported_distance_matrix,
)
import bijux_phylogenetics.distance.average_linkage as average_linkage_module
import bijux_phylogenetics.distance.upgma as upgma_module
from bijux_phylogenetics.io.newick import dumps_newick
from bijux_phylogenetics.runtime.errors import InvalidDistanceMatrixError

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


def _lookup_from_imported_matrix(
    path: Path,
) -> tuple[list[str], dict[tuple[str, str], float]]:
    entries = load_imported_distance_matrix(path)
    identifiers = sorted(
        {entry.left_identifier for entry in entries}
        | {entry.right_identifier for entry in entries}
    )
    lookup = {(identifier, identifier): 0.0 for identifier in identifiers}
    for entry in entries:
        lookup[(entry.left_identifier, entry.right_identifier)] = entry.distance
    return identifiers, lookup


def test_build_wpgma_tree_matches_uneven_cluster_fixture_and_differs_from_upgma() -> (
    None
):
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_wpgma_uneven_cluster.tsv")
    )
    upgma_tree, upgma_report = build_upgma_tree(identifiers, lookup)
    wpgma_tree, wpgma_report = build_wpgma_tree(identifiers, lookup)
    assert dumps_newick(upgma_tree) == (
        "((((A:0.5,D:0.5)Inner1:2.25,E:2.75)Inner2:0.583333333333333,C:3.33333333333333)Inner3:0.541666666666667,B:3.875)Inner4;"
    )
    assert dumps_newick(wpgma_tree) == (
        "((((A:0.5,D:0.5)Inner1:2.25,E:2.75)Inner2:0.75,C:3.5)Inner3:0.5625,B:4.0625)Inner4;"
    )
    assert dumps_newick(wpgma_tree) != dumps_newick(upgma_tree)
    assert [
        (row.left_cluster, row.right_cluster, row.pair_distance, row.merge_height)
        for row in upgma_report.merge_history
    ] == [
        ("A", "D", 1.0, 0.5),
        ("A|D", "E", 5.5, 2.75),
        ("A|D|E", "C", 6.666666666666667, 3.3333333333333335),
        ("A|C|D|E", "B", 7.75, 3.875),
    ]
    assert [
        (row.left_cluster, row.right_cluster, row.pair_distance, row.merge_height)
        for row in wpgma_report.merge_history
    ] == [
        ("A", "D", 1.0, 0.5),
        ("A|D", "E", 5.5, 2.75),
        ("A|D|E", "C", 7.0, 3.5),
        ("A|C|D|E", "B", 8.125, 4.0625),
    ]
    assert [row.height for row in upgma_report.cluster_heights] == [
        0.5,
        2.75,
        3.3333333333333335,
        3.875,
    ]
    assert [row.height for row in wpgma_report.cluster_heights] == [
        0.5,
        2.75,
        3.5,
        4.0625,
    ]
    assert wpgma_report.ultrametric_compatible is False
    assert wpgma_report.assumption_warnings == [
        "pairwise distances are not ultrametric, so wpgma's strict clock-like assumption is violated"
    ]


def test_build_wpgma_tree_does_not_delegate_to_upgma(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_wpgma_uneven_cluster.tsv")
    )

    def _boom(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("wpgma must not call build_upgma_tree")

    monkeypatch.setattr(upgma_module, "build_upgma_tree", _boom)
    tree, report = build_wpgma_tree(identifiers, lookup)
    assert dumps_newick(tree).endswith(";")
    assert report.merge_history[-1].resulting_cluster == "A|B|C|D|E"


def test_build_tree_from_imported_distance_matrix_supports_wpgma() -> None:
    tree, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_wpgma_uneven_cluster.tsv"),
        method="wpgma",
    )
    assert dumps_newick(tree) == (
        "((((A:0.5,D:0.5)Inner1:2.25,E:2.75)Inner2:0.75,C:3.5)Inner3:0.5625,B:4.0625)Inner4;"
    )
    assert report.method == "wpgma"
    assert report.method_policy.method == "wpgma"
    assert report.assumptions.ultrametric_compatible is False


def test_build_tree_from_imported_distance_matrix_rejects_asymmetric_wpgma_input() -> (
    None
):
    with pytest.raises(InvalidDistanceMatrixError) as error:
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_asymmetric.tsv"),
            method="wpgma",
        )
    assert error.value.code == "invalid_distance_matrix_error"
    assert "asymmetric directional entries" in str(error.value)


def test_build_wpgma_tree_routes_through_shared_agglomerative_engine(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_wpgma_uneven_cluster.tsv")
    )

    def _boom(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("wpgma must call the shared agglomerative engine")

    monkeypatch.setattr(
        average_linkage_module,
        "build_agglomerative_clustering_tree",
        _boom,
    )
    with pytest.raises(AssertionError, match="shared agglomerative engine"):
        build_wpgma_tree(identifiers, lookup)

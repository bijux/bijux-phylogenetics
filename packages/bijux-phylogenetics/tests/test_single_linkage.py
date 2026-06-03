from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    build_single_linkage_tree,
    build_tree_from_imported_distance_matrix,
    build_upgma_tree,
    build_wpgma_tree,
    load_imported_distance_matrix,
)
import bijux_phylogenetics.distance.extremal_linkage as extremal_linkage_module
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


def test_build_single_linkage_tree_matches_chaining_fixture() -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_single_linkage_chain.tsv")
    )
    tree, report = build_single_linkage_tree(identifiers, lookup)
    assert dumps_newick(tree) == (
        "((((A:0.5,B:0.5)Inner1:0.5,C:1)Inner2:0.5,D:1.5)Inner3:0.5,E:2)Inner4;"
    )
    assert [
        (
            row.merge_index,
            row.left_cluster,
            row.right_cluster,
            row.pair_distance,
            row.merge_height,
            row.resulting_cluster,
        )
        for row in report.merge_history
    ] == [
        (1, "A", "B", 1.0, 0.5, "A|B"),
        (2, "A|B", "C", 2.0, 1.0, "A|B|C"),
        (3, "A|B|C", "D", 3.0, 1.5, "A|B|C|D"),
        (4, "A|B|C|D", "E", 4.0, 2.0, "A|B|C|D|E"),
    ]
    assert [(row.cluster, row.height) for row in report.cluster_heights] == [
        ("A|B", 0.5),
        ("A|B|C", 1.0),
        ("A|B|C|D", 1.5),
        ("A|B|C|D|E", 2.0),
    ]


def test_build_single_linkage_tree_differs_from_average_linkage_on_chain_fixture() -> (
    None
):
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_single_linkage_chain.tsv")
    )
    single_tree, single_report = build_single_linkage_tree(identifiers, lookup)
    upgma_tree, upgma_report = build_upgma_tree(identifiers, lookup)
    wpgma_tree, wpgma_report = build_wpgma_tree(identifiers, lookup)
    assert dumps_newick(single_tree) != dumps_newick(upgma_tree)
    assert dumps_newick(single_tree) != dumps_newick(wpgma_tree)
    assert [row.pair_distance for row in single_report.merge_history] == [
        1.0,
        2.0,
        3.0,
        4.0,
    ]
    assert [row.pair_distance for row in upgma_report.merge_history] != [
        1.0,
        2.0,
        3.0,
        4.0,
    ]
    assert [row.pair_distance for row in wpgma_report.merge_history] != [
        1.0,
        2.0,
        3.0,
        4.0,
    ]


def test_build_single_linkage_tree_routes_through_shared_agglomerative_engine(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_single_linkage_chain.tsv")
    )

    def _boom(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("single-linkage must call the shared agglomerative engine")

    monkeypatch.setattr(
        extremal_linkage_module,
        "build_agglomerative_clustering_tree",
        _boom,
    )
    with pytest.raises(AssertionError, match="shared agglomerative engine"):
        build_single_linkage_tree(identifiers, lookup)


def test_build_tree_from_imported_distance_matrix_supports_single_linkage() -> None:
    tree, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_single_linkage_chain.tsv"),
        method="single-linkage",
    )
    assert dumps_newick(tree) == (
        "((((A:0.5,B:0.5)Inner1:0.5,C:1)Inner2:0.5,D:1.5)Inner3:0.5,E:2)Inner4;"
    )
    assert report.method == "single-linkage"
    assert report.method_policy.method == "single-linkage"


def test_build_tree_from_imported_distance_matrix_rejects_asymmetric_single_linkage_input() -> (
    None
):
    with pytest.raises(InvalidDistanceMatrixError) as error:
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_asymmetric.tsv"),
            method="single-linkage",
        )
    assert error.value.code == "invalid_distance_matrix_error"
    assert "asymmetric directional entries" in str(error.value)

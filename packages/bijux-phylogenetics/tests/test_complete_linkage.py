from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    build_complete_linkage_tree,
    build_single_linkage_tree,
    build_tree_from_imported_distance_matrix,
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


def test_build_complete_linkage_tree_matches_compact_cluster_fixture() -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv")
    )
    tree, report = build_complete_linkage_tree(identifiers, lookup)
    assert dumps_newick(tree) == (
        "(((A:1,D:1)Inner2:2.5,E:3.5)Inner3:2,(B:1,C:1)Inner1:4.5)Inner4;"
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
        (1, "B", "C", 2.0, 1.0, "B|C"),
        (2, "A", "D", 2.0, 1.0, "A|D"),
        (3, "A|D", "E", 7.0, 3.5, "A|D|E"),
        (4, "A|D|E", "B|C", 11.0, 5.5, "A|B|C|D|E"),
    ]
    assert [(row.cluster, row.height) for row in report.cluster_heights] == [
        ("B|C", 1.0),
        ("A|D", 1.0),
        ("A|D|E", 3.5),
        ("A|B|C|D|E", 5.5),
    ]


def test_build_complete_linkage_tree_differs_from_single_linkage_on_compact_cluster_fixture() -> (
    None
):
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv")
    )
    complete_tree, complete_report = build_complete_linkage_tree(identifiers, lookup)
    single_tree, single_report = build_single_linkage_tree(identifiers, lookup)
    assert dumps_newick(complete_tree) == (
        "(((A:1,D:1)Inner2:2.5,E:3.5)Inner3:2,(B:1,C:1)Inner1:4.5)Inner4;"
    )
    assert dumps_newick(single_tree) == (
        "(((A:1,D:1)Inner2:0.5,(B:1,C:1)Inner1:0.5)Inner3:1,E:2.5)Inner4;"
    )
    assert [row.pair_distance for row in complete_report.merge_history] == [
        2.0,
        2.0,
        7.0,
        11.0,
    ]
    assert [row.pair_distance for row in single_report.merge_history] == [
        2.0,
        2.0,
        3.0,
        5.0,
    ]


def test_build_complete_linkage_tree_routes_through_shared_agglomerative_engine(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv")
    )

    def _boom(*_args: object, **_kwargs: object) -> object:
        raise AssertionError(
            "complete-linkage must call the shared agglomerative engine"
        )

    monkeypatch.setattr(
        extremal_linkage_module,
        "build_agglomerative_clustering_tree",
        _boom,
    )
    with pytest.raises(AssertionError, match="shared agglomerative engine"):
        build_complete_linkage_tree(identifiers, lookup)


def test_build_tree_from_imported_distance_matrix_supports_complete_linkage() -> None:
    tree, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv"),
        method="complete-linkage",
    )
    assert dumps_newick(tree) == (
        "(((A:1,D:1)Inner2:2.5,E:3.5)Inner3:2,(B:1,C:1)Inner1:4.5)Inner4;"
    )
    assert report.method == "complete-linkage"
    assert report.method_policy.method == "complete-linkage"


def test_build_tree_from_imported_distance_matrix_rejects_asymmetric_complete_linkage_input() -> (
    None
):
    with pytest.raises(InvalidDistanceMatrixError) as error:
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_asymmetric.tsv"),
            method="complete-linkage",
        )
    assert error.value.code == "invalid_distance_matrix_error"
    assert "asymmetric directional entries" in str(error.value)

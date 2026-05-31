from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    build_tree_from_imported_distance_matrix,
    build_upgma_tree,
    load_imported_distance_matrix,
)
import bijux_phylogenetics.distance.average_linkage as average_linkage_module
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


def test_build_upgma_tree_matches_hand_computed_five_taxon_fixture() -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_upgma_five_taxon.tsv")
    )
    tree, report = build_upgma_tree(identifiers, lookup)
    assert tree.rooted is True
    assert dumps_newick(tree) == (
        "(((A:1,B:1)Inner1:2,C:3)Inner3:2,(D:2,E:2)Inner2:3)Inner4;"
    )
    assert tree.root_to_tip_pairs() == [
        ("A", 5.0),
        ("B", 5.0),
        ("C", 5.0),
        ("D", 5.0),
        ("E", 5.0),
    ]
    assert [
        (
            row.merge_index,
            row.left_cluster,
            row.right_cluster,
            row.left_cluster_size,
            row.right_cluster_size,
            row.pair_distance,
            row.merge_height,
            row.resulting_cluster,
            row.resulting_cluster_size,
        )
        for row in report.merge_history
    ] == [
        (1, "A", "B", 1, 1, 2.0, 1.0, "A|B", 2),
        (2, "D", "E", 1, 1, 4.0, 2.0, "D|E", 2),
        (3, "A|B", "C", 2, 1, 6.0, 3.0, "A|B|C", 3),
        (4, "A|B|C", "D|E", 3, 2, 10.0, 5.0, "A|B|C|D|E", 5),
    ]
    assert [
        (row.merge_index, row.cluster, row.height) for row in report.cluster_heights
    ] == [
        (1, "A|B", 1.0),
        (2, "D|E", 2.0),
        (3, "A|B|C", 3.0),
        (4, "A|B|C|D|E", 5.0),
    ]
    assert report.ultrametric_compatible is True
    assert report.assumption_warnings == []


def test_build_tree_from_imported_distance_matrix_surfaces_upgma_warning_for_nonultrametric_input() -> (
    None
):
    tree, report = build_tree_from_imported_distance_matrix(
        fixture("example_distance_matrix_nonultrametric.tsv"),
        method="upgma",
    )
    assert tree.rooted is True
    assert report.assumptions.ultrametric_compatible is False
    assert report.assumptions.warnings[-1] == (
        "pairwise distances are not ultrametric, so UPGMA's strict clock-like assumption is violated"
    )


def test_build_tree_from_imported_distance_matrix_rejects_asymmetric_upgma_input() -> (
    None
):
    try:
        build_tree_from_imported_distance_matrix(
            fixture("example_distance_matrix_asymmetric.tsv"),
            method="upgma",
        )
    except InvalidDistanceMatrixError as error:
        assert error.code == "invalid_distance_matrix_error"
        assert "asymmetric directional entries" in str(error)
    else:  # pragma: no cover - defensive assertion
        raise AssertionError("expected InvalidDistanceMatrixError")


def test_build_upgma_tree_routes_through_shared_agglomerative_engine(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_upgma_five_taxon.tsv")
    )

    def _boom(*_args: object, **_kwargs: object) -> object:
        raise AssertionError("upgma must call the shared agglomerative engine")

    monkeypatch.setattr(
        average_linkage_module,
        "build_agglomerative_clustering_tree",
        _boom,
    )
    with pytest.raises(AssertionError, match="shared agglomerative engine"):
        build_upgma_tree(identifiers, lookup)

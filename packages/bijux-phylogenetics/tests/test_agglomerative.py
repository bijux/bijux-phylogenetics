from __future__ import annotations

from pathlib import Path

from bijux_phylogenetics.distance import load_imported_distance_matrix
from bijux_phylogenetics.distance.agglomerative import (
    build_agglomerative_clustering_tree,
)
from bijux_phylogenetics.io.newick import dumps_newick

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


def test_agglomerative_tree_supports_taxon_count_updates_for_upgma_fixture() -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_upgma_five_taxon.tsv")
    )
    result = build_agglomerative_clustering_tree(
        identifiers,
        lookup,
        method_name="upgma",
        update_rule="taxon-count",
        assess_ultrametric_assumption=True,
    )
    assert dumps_newick(result.tree) == (
        "(((A:1,B:1)Inner1:2,C:3)Inner3:2,(D:2,E:2)Inner2:3)Inner4;"
    )
    assert [row.pair_distance for row in result.merge_history] == [2.0, 4.0, 6.0, 10.0]
    assert result.ultrametric_compatible is True
    assert result.assumption_warnings == []


def test_agglomerative_tree_supports_equal_updates_for_wpgma_fixture() -> None:
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_wpgma_uneven_cluster.tsv")
    )
    result = build_agglomerative_clustering_tree(
        identifiers,
        lookup,
        method_name="wpgma",
        update_rule="equal",
        assess_ultrametric_assumption=True,
    )
    assert dumps_newick(result.tree) == (
        "((((A:0.5,D:0.5)Inner1:2.25,E:2.75)Inner2:0.75,C:3.5)Inner3:0.5625,B:4.0625)Inner4;"
    )
    assert [row.pair_distance for row in result.merge_history] == [1.0, 5.5, 7.0, 8.125]
    assert result.ultrametric_compatible is False
    assert result.assumption_warnings == [
        "pairwise distances are not ultrametric, so wpgma's strict clock-like assumption is violated"
    ]


def test_agglomerative_tree_supports_minimum_updates_for_single_linkage_fixture() -> (
    None
):
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_single_linkage_chain.tsv")
    )
    result = build_agglomerative_clustering_tree(
        identifiers,
        lookup,
        method_name="single-linkage",
        update_rule="minimum",
        assess_ultrametric_assumption=False,
    )
    assert dumps_newick(result.tree) == (
        "((((A:0.5,B:0.5)Inner1:0.5,C:1)Inner2:0.5,D:1.5)Inner3:0.5,E:2)Inner4;"
    )
    assert [row.pair_distance for row in result.merge_history] == [1.0, 2.0, 3.0, 4.0]
    assert result.ultrametric_compatible is None
    assert result.assumption_warnings == []


def test_agglomerative_tree_supports_maximum_updates_for_complete_linkage_fixture() -> (
    None
):
    identifiers, lookup = _lookup_from_imported_matrix(
        fixture("example_distance_matrix_complete_linkage_compact_cluster.tsv")
    )
    result = build_agglomerative_clustering_tree(
        identifiers,
        lookup,
        method_name="complete-linkage",
        update_rule="maximum",
        assess_ultrametric_assumption=False,
    )
    assert dumps_newick(result.tree) == (
        "(((A:1,D:1)Inner2:2.5,E:3.5)Inner3:2,(B:1,C:1)Inner1:4.5)Inner4;"
    )
    assert [row.pair_distance for row in result.merge_history] == [2.0, 2.0, 7.0, 11.0]
    assert result.ultrametric_compatible is None
    assert result.assumption_warnings == []

from __future__ import annotations

from bijux_phylogenetics.clades import extract_tree_set_clades
from bijux_phylogenetics.shared_tree_set_fixtures import (
    get_shared_tree_set_fixture,
    list_shared_tree_set_fixtures,
)


def test_shared_tree_set_fixture_catalog_covers_governed_multiple_tree_case() -> None:
    fixtures = list_shared_tree_set_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {"multiple-trees", "branch-lengths", "shared-taxon-set"} <= feature_tags


def test_shared_tree_set_fixture_lookup_preserves_durable_ids() -> None:
    fixture = get_shared_tree_set_fixture("basic_newick_tree_set")

    assert fixture.relative_path == "trees/example_tree_set_left.nwk"
    assert fixture.tree_count == 3
    assert fixture.shared_taxa == ("A", "B", "C", "D")
    assert fixture.path.is_file()


def test_shared_tree_set_fixture_catalog_loads_tree_set_clades() -> None:
    fixture = get_shared_tree_set_fixture("basic_newick_tree_set")

    report = extract_tree_set_clades(fixture.path)

    assert report.tree_count == 3
    assert report.source_format == "newick"
    assert sorted({row.tree_index for row in report.rows if row.tree_index is not None}) == [
        1,
        2,
        3,
    ]
    assert sorted(
        {row.node_label for row in report.rows if row.tree_index == 1 and row.node_kind == "tip"}
    ) == ["A", "B", "C", "D"]

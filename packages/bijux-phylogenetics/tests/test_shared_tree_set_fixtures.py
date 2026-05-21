from __future__ import annotations

import pytest

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_tree_set_fixture,
    list_shared_tree_set_fixtures,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError
from bijux_phylogenetics.trees import extract_tree_set_clades


def test_shared_tree_set_fixture_catalog_covers_governed_multiple_tree_case() -> None:
    fixtures = list_shared_tree_set_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert {
        "multiple-trees",
        "branch-lengths",
        "shared-taxon-set",
        "topology-distance",
        "consensus",
        "prop-clades",
        "large-tip-count",
    } <= feature_tags
    assert max(len(fixture.shared_taxa) for fixture in fixtures) >= 100


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
    assert sorted(
        {row.tree_index for row in report.rows if row.tree_index is not None}
    ) == [
        1,
        2,
        3,
    ]
    assert sorted(
        {
            row.node_label
            for row in report.rows
            if row.tree_index == 1 and row.node_kind == "tip"
        }
    ) == ["A", "B", "C", "D"]


def test_shared_tree_set_fixture_catalog_preserves_large_topology_distance_pair() -> (
    None
):
    fixture = get_shared_tree_set_fixture("topology_distance_large_rooted_pair")

    report = extract_tree_set_clades(fixture.path)

    assert fixture.tree_count == 2
    assert len(fixture.shared_taxa) == 128
    assert report.tree_count == 2
    assert sorted(
        {row.tree_index for row in report.rows if row.tree_index is not None}
    ) == [1, 2]


def test_shared_tree_set_fixture_catalog_preserves_consensus_posterior_fixture() -> (
    None
):
    fixture = get_shared_tree_set_fixture("consensus_posterior_six_taxon_tree_set")

    report = extract_tree_set_clades(fixture.path)

    assert fixture.tree_count == 5
    assert fixture.shared_taxa == ("A", "B", "C", "D", "E", "F")
    assert report.tree_count == 5


def test_shared_tree_set_fixture_catalog_preserves_consensus_failure_fixture() -> None:
    fixture = get_shared_tree_set_fixture("consensus_mismatched_taxon_tree_set")

    assert fixture.shared_taxa == ()
    assert fixture.tree_count == 2
    with pytest.raises(
        InvalidAlignmentError,
        match="requires identical taxon sets",
    ):
        extract_tree_set_clades(fixture.path)


def test_shared_tree_set_fixture_catalog_preserves_prop_clades_posterior_fixture() -> (
    None
):
    fixture = get_shared_tree_set_fixture("prop_clades_posterior_six_taxon_tree_set")

    report = extract_tree_set_clades(fixture.path)

    assert fixture.tree_count == 5
    assert fixture.shared_taxa == ("A", "B", "C", "D", "E", "F")
    assert report.tree_count == 5


def test_shared_tree_set_fixture_catalog_preserves_prop_clades_failure_fixture() -> (
    None
):
    fixture = get_shared_tree_set_fixture("prop_clades_mismatched_taxon_tree_set")

    assert fixture.shared_taxa == ()
    assert fixture.tree_count == 2
    with pytest.raises(
        InvalidAlignmentError,
        match="requires identical taxon sets",
    ):
        extract_tree_set_clades(fixture.path)

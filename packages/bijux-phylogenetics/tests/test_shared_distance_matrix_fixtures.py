from __future__ import annotations

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_distance_matrix_fixture,
    list_shared_distance_matrix_fixtures,
)


def test_shared_distance_matrix_fixture_catalog_lists_neighbor_joining_corpus() -> None:
    fixtures = list_shared_distance_matrix_fixtures()
    assert [fixture.fixture_id for fixture in fixtures] == [
        "analytical_three_taxon",
        "ultrametric_four_taxon",
        "nonultrametric_four_taxon",
    ]
    assert fixtures[0].taxa == ("A", "B", "C")
    assert fixtures[1].feature_tags == (
        "four-taxon",
        "ultrametric",
        "neighbor-joining-tie-case",
    )


def test_get_shared_distance_matrix_fixture_resolves_existing_file() -> None:
    fixture = get_shared_distance_matrix_fixture("nonultrametric_four_taxon")
    assert fixture.path.exists()
    assert fixture.pair_count == 16
    assert "branch-length-case" in fixture.feature_tags

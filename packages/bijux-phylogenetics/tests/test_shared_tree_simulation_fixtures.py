from __future__ import annotations

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_tree_simulation_fixture,
    list_shared_tree_simulation_fixtures,
)


def test_shared_tree_simulation_fixture_catalog_lists_governed_parity_cases() -> None:
    fixtures = list_shared_tree_simulation_fixtures()

    assert [fixture.fixture_id for fixture in fixtures] == [
        "rtree_rooted_six_taxon_uniform_64",
        "rtree_rooted_twelve_taxon_uniform_128",
        "rcoal_rooted_six_taxon_64",
        "rcoal_rooted_twelve_taxon_128",
    ]
    assert fixtures[0].reference_function == "ape::rtree"
    assert fixtures[1].replicate_count == 128
    assert fixtures[2].population_size == 1.0
    assert "distributional-parity" in fixtures[3].feature_tags


def test_get_shared_tree_simulation_fixture_resolves_catalog_parameters() -> None:
    fixture = get_shared_tree_simulation_fixture("rcoal_rooted_twelve_taxon_128")

    assert fixture.path.exists()
    assert fixture.simulation_model == "coalescent"
    assert fixture.tip_count == 12
    assert fixture.seed == 29

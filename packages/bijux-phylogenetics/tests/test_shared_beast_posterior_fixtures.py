from __future__ import annotations

from bijux_phylogenetics.datasets.shared_fixtures import (
    get_shared_beast_posterior_fixture,
    list_shared_beast_posterior_fixtures,
)


def test_shared_beast_posterior_fixture_catalog_covers_governed_real_artifact_case() -> (
    None
):
    fixtures = list_shared_beast_posterior_fixtures()
    feature_tags = {tag for fixture in fixtures for tag in fixture.feature_tags}

    assert len(fixtures) >= 1
    assert {
        "real-artifact",
        "analysis-xml",
        "posterior-log",
        "posterior-trees",
        "consensus",
        "mcc",
        "burnin-calibrated",
        "ess-anchored",
    } <= feature_tags


def test_shared_beast_posterior_fixture_lookup_preserves_reference_bundle_paths() -> (
    None
):
    fixture = get_shared_beast_posterior_fixture("strict_yule_real_posterior")

    assert (
        fixture.analysis_xml_relative_path
        == "metadata/beast2_strict_yule_posterior.xml"
    )
    assert (
        fixture.posterior_log_relative_path
        == "metadata/beast2_strict_yule_posterior.log"
    )
    assert (
        fixture.posterior_trees_relative_path
        == "metadata/beast2_strict_yule_posterior.trees"
    )
    assert fixture.beast_version == "2.7"
    assert fixture.recommended_burnin_fraction == 0.1
    assert fixture.posterior_row_count == 101
    assert fixture.posterior_tree_count == 101
    assert fixture.shared_taxa == ("A", "B", "C", "D")
    assert fixture.analysis_xml_path.is_file()
    assert fixture.posterior_log_path.is_file()
    assert fixture.posterior_trees_path.is_file()
    assert fixture.consensus_tree_path.is_file()
    assert fixture.mcc_tree_path.is_file()
    assert fixture.reference_json_path.is_file()


def test_shared_beast_posterior_fixture_reference_loads_burnin_consensus_and_ess_expectations() -> (
    None
):
    reference = get_shared_beast_posterior_fixture(
        "strict_yule_real_posterior"
    ).load_reference()

    assert sorted(reference.burnin_reference) == [0.0, 0.1, 0.25, 0.5]
    assert reference.burnin_reference[0.1].kept_row_count == 91
    assert reference.burnin_reference[0.1].kept_tree_count == 91
    assert reference.consensus_reference.annotated_node_count == 2
    assert reference.consensus_reference.minimum_posterior_probability > 0.5
    assert reference.mcc_reference.selected_tree_index == 13
    assert sorted(reference.parameter_reference) == [
        "birthRate",
        "clockRate",
        "likelihood",
        "posterior",
        "prior",
        "tree.height",
    ]
    assert reference.parameter_reference["posterior"].effective_sample_size > 4.0

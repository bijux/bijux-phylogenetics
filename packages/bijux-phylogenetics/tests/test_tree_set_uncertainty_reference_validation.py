from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_tree_set_uncertainty_reference_fixtures,
)


def test_validate_tree_set_uncertainty_reference_fixtures_governs_publication_audits() -> (
    None
):
    report = validate_tree_set_uncertainty_reference_fixtures()

    assert report.goal_id == 231
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["multi_topology_tree_set_uncertainty_visible"].observed[
            "publication_ready"
        ]
        is True
    )
    assert (
        observed["single_topology_tree_set_keeps_empty_instability_panel"].observed[
            "publication_ready"
        ]
        is True
    )
    assert (
        observed["single_topology_tree_set_keeps_empty_instability_panel"].observed[
            "unstable_taxon_count"
        ]
        == 0
    )

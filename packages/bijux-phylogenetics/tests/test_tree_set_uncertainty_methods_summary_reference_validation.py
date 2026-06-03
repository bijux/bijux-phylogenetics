from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_tree_set_uncertainty_methods_summary_reference_fixtures,
)


def test_validate_tree_set_uncertainty_methods_summary_reference_fixtures_governs_uncertainty_counts() -> (
    None
):
    report = validate_tree_set_uncertainty_methods_summary_reference_fixtures()

    assert report.goal_id == 253
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed[
            "multi_topology_methods_summary_keeps_support_and_instability"
        ].observed["warning_count"]
        == 4
    )
    assert (
        observed[
            "multi_topology_methods_summary_keeps_support_and_instability"
        ].observed["multimodal"]
        is True
    )
    assert (
        observed["single_topology_methods_summary_keeps_clean_lane_explicit"].observed[
            "warning_count"
        ]
        == 0
    )
    assert (
        observed["single_topology_methods_summary_keeps_clean_lane_explicit"].observed[
            "topology_cluster_count"
        ]
        == 1
    )

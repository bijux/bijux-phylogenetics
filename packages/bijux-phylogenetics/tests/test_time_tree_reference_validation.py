from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_time_tree_reference_fixtures,
)


def test_validate_time_tree_reference_fixtures_governs_visible_uncertainty() -> None:
    report = validate_time_tree_reference_fixtures()

    assert report.goal_id == 228
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["visible_uncertainty_time_tree"].observed["publication_ready"] is True
    )
    assert (
        observed["visible_uncertainty_time_tree"].observed["readiness_decision"]
        == "ready"
    )
    assert (
        observed["invalid_tip_dates_block_publication"].observed["publication_ready"]
        is False
    )
    assert (
        observed["invalid_tip_dates_block_publication"].observed[
            "has_tip_date_limitation"
        ]
        is True
    )

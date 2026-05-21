from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_biogeography_figure_reference_fixtures,
)


def test_validate_biogeography_figure_reference_fixtures_governs_publication_audits() -> (
    None
):
    report = validate_biogeography_figure_reference_fixtures()

    assert report.goal_id == 230
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["visible_state_probabilities_and_transitions"].observed[
            "publication_ready"
        ]
        is True
    )
    assert (
        observed["missing_centroid_blocks_publication"].observed["publication_ready"]
        is False
    )
    assert (
        observed["missing_centroid_blocks_publication"].observed[
            "has_exclusion_limitation"
        ]
        is True
    )

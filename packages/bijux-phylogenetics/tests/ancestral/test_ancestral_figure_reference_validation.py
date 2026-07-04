from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_ancestral_figure_reference_fixtures,
)


def test_validate_ancestral_figure_reference_fixtures_governs_uncertainty_visibility() -> (
    None
):
    report = validate_ancestral_figure_reference_fixtures()

    assert report.goal_id == 229
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["continuous_intervals_visible"].observed["publication_ready"] is True
    )
    assert (
        observed["continuous_intervals_visible"].observed["uncertainty_visible"] is True
    )
    assert (
        observed["discrete_probabilities_interpretable"].observed["publication_ready"]
        is True
    )
    assert (
        observed["discrete_probabilities_interpretable"].observed[
            "rendered_internal_pie_count_matches"
        ]
        is True
    )

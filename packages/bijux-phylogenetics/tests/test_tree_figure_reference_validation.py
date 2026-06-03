from __future__ import annotations

from bijux_phylogenetics.validation.reference import validate_figure_reference_fixtures


def test_validate_figure_reference_fixtures_governs_publication_audits() -> None:
    report = validate_figure_reference_fixtures()

    assert report.goal_id == 95
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert observed["validated_support_figure"].observed["legend_complete"] is True
    assert observed["validated_support_figure"].observed["caption_ready"] is True
    assert observed["validated_support_figure"].observed["legible"] is True
    assert (
        observed["withheld_invalid_support_figure"].observed["legend_complete"] is True
    )

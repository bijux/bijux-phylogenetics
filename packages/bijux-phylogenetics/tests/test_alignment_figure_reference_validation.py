from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_alignment_figure_reference_fixtures,
)


def test_validate_alignment_figure_reference_fixtures_governs_publication_audits() -> (
    None
):
    report = validate_alignment_figure_reference_fixtures()

    assert report.goal_id == 232
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["clean_alignment_quality_figures_ready"].observed["publication_ready"]
        is True
    )
    assert (
        observed["missingness_alignment_quality_figures_blocked"].observed[
            "publication_ready"
        ]
        is False
    )
    assert (
        observed["missingness_alignment_quality_figures_blocked"].observed[
            "suspicious_alignment"
        ]
        is True
    )

from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_diversification_figure_reference_fixtures,
)


def test_validate_diversification_figure_reference_fixtures_governs_publication_audits() -> (
    None
):
    report = validate_diversification_figure_reference_fixtures()

    assert report.goal_id == 233
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["sampling_complete_diversification_figures_ready"].observed[
            "publication_ready"
        ]
        is True
    )
    assert (
        observed["sampling_complete_diversification_figures_ready"].observed[
            "better_model"
        ]
        == "yule"
    )
    assert (
        observed["sampling_complete_diversification_figures_ready"].observed[
            "methods_summary_present"
        ]
        is True
    )
    assert (
        observed["incomplete_sampling_blocks_diversification_readiness"].observed[
            "publication_ready"
        ]
        is False
    )
    assert (
        observed["incomplete_sampling_blocks_diversification_readiness"].observed[
            "sampling_metadata_complete"
        ]
        is False
    )
    assert (
        observed["incomplete_sampling_blocks_diversification_readiness"].observed[
            "methods_summary_present"
        ]
        is True
    )

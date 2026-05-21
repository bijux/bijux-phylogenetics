from __future__ import annotations

import pytest

from bijux_phylogenetics.validation.reference import (
    validate_comparative_model_figure_reference_fixtures,
)


@pytest.mark.slow
def test_validate_comparative_model_figure_reference_fixtures_governs_publication_audits() -> (
    None
):
    report = validate_comparative_model_figure_reference_fixtures()

    assert report.goal_id == 234
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["signal_strong_twenty_four_taxa_model_figures_ready"].observed[
            "publication_ready"
        ]
        is True
    )
    assert (
        observed["signal_strong_twenty_four_taxa_model_figures_ready"].observed[
            "selected_model"
        ]
        == "brownian"
    )
    assert (
        observed[
            "signal_strong_one_hundred_twenty_eight_taxa_support_ambiguous"
        ].observed["publication_ready"]
        is False
    )
    assert (
        observed[
            "signal_strong_one_hundred_twenty_eight_taxa_support_ambiguous"
        ].observed["support_distinct"]
        is False
    )

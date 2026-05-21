from __future__ import annotations

from bijux_phylogenetics.validation.reference import (
    validate_trait_tree_reference_fixtures,
)


def test_validate_trait_tree_reference_fixtures_governs_publication_audits() -> None:
    report = validate_trait_tree_reference_fixtures()

    assert report.goal_id == 227
    assert report.passed is True
    observed = {fixture.name: fixture for fixture in report.fixtures}
    assert (
        observed["complete_annotated_trait_tree"].observed["publication_ready"] is True
    )
    assert (
        observed["incomplete_metadata_strip_blocks_publication"].observed[
            "publication_ready"
        ]
        is False
    )
    assert observed["incomplete_metadata_strip_blocks_publication"].observed[
        "location_missing_taxa"
    ] == ["C"]

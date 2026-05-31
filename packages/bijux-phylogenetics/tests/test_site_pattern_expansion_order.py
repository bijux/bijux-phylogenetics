from __future__ import annotations

import pytest

from bijux_phylogenetics.phylo.alignment.models import AlignmentRecord
from bijux_phylogenetics.phylo.likelihood import (
    compress_alignment_site_patterns_from_records,
)
from bijux_phylogenetics.phylo.likelihood.patterns import (
    AlignmentSitePattern,
    CompressedAlignmentSitePatterns,
)
from bijux_phylogenetics.phylo.likelihood.sites import (
    expanded_site_log_likelihood_rows_from_patterns,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError


def test_expanded_site_log_likelihood_rows_follow_alignment_site_order() -> None:
    records = [
        AlignmentRecord(identifier="A", sequence="ACAC"),
        AlignmentRecord(identifier="B", sequence="GTGT"),
    ]
    compressed_patterns = compress_alignment_site_patterns_from_records(records)

    rows, total = expanded_site_log_likelihood_rows_from_patterns(
        compressed_patterns,
        site_log_likelihood=lambda states: 0.0 if states[0] == "A" else -1.0,
    )

    assert [row.site_position for row in rows] == [1, 2, 3, 4]
    assert [row.pattern_id for row in rows] == [
        "pattern-1",
        "pattern-2",
        "pattern-1",
        "pattern-2",
    ]
    assert total == -2.0


def test_expanded_site_log_likelihood_rows_reject_duplicate_site_positions() -> None:
    compressed_patterns = CompressedAlignmentSitePatterns(
        source_path=None,
        taxon_order=["A", "B"],
        alignment_length=3,
        pattern_count=2,
        patterns=[
            AlignmentSitePattern(
                pattern_id="pattern-1",
                states=("A", "A"),
                weight=2,
                site_positions=[1, 3],
            ),
            AlignmentSitePattern(
                pattern_id="pattern-2",
                states=("C", "C"),
                weight=1,
                site_positions=[3],
            ),
        ],
    )

    with pytest.raises(
        InvalidAlignmentError,
        match="cover the alignment exactly once",
    ):
        expanded_site_log_likelihood_rows_from_patterns(
            compressed_patterns,
            site_log_likelihood=lambda _states: 0.0,
        )

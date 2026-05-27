from __future__ import annotations

from pathlib import Path

import pytest

from bijux_phylogenetics.distance import (
    build_distance_tree,
    diagnose_distance_saturation,
)
from bijux_phylogenetics.runtime.errors import InvalidAlignmentError

FIXTURES = Path(__file__).parent / "fixtures"
FIXTURE_GROUPS = ("trees", "alignments", "metadata", "expected")


def fixture(name: str) -> Path:
    direct = FIXTURES / name
    if direct.exists():
        return direct
    for group in FIXTURE_GROUPS:
        candidate = FIXTURES / group / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(name)


def test_diagnose_distance_saturation_reports_jc69_pair_level_warnings() -> None:
    report = diagnose_distance_saturation(
        fixture("example_alignment_distance_saturated.fasta"),
        model="jc69",
    )

    assert report.blocks_tree_inference is True
    assert report.blocking_warning_count == 2
    assert [
        (
            row.left_identifier,
            row.right_identifier,
            row.warning_kind,
            row.blocks_tree_inference,
            row.reason,
        )
        for row in report.warning_rows
    ] == [
        (
            "A",
            "B",
            "undefined-corrected-distance",
            True,
            "p-distance exceeds the Jukes-Cantor correction range, so the corrected distance is undefined",
        ),
        (
            "B",
            "C",
            "infinite-corrected-distance",
            True,
            "p-distance is at the Jukes-Cantor correction limit, so the corrected distance tends to infinity",
        ),
    ]


def test_diagnose_distance_saturation_reports_k80_pair_level_warnings() -> None:
    report = diagnose_distance_saturation(
        fixture("example_alignment_distance.fasta"),
        model="k80",
    )

    assert report.blocks_tree_inference is True
    assert report.blocking_warning_count == 4
    assert [
        (row.left_identifier, row.right_identifier, row.warning_kind)
        for row in report.warning_rows
    ] == [
        ("A", "C", "infinite-corrected-distance"),
        ("A", "D", "infinite-corrected-distance"),
        ("B", "C", "infinite-corrected-distance"),
        ("B", "D", "infinite-corrected-distance"),
    ]


def test_diagnose_distance_saturation_reports_tn93_pair_level_warnings() -> None:
    report = diagnose_distance_saturation(
        fixture("example_alignment_distance.fasta"),
        model="tn93",
    )

    assert report.blocks_tree_inference is True
    assert report.blocking_warning_count == 6
    assert report.warnings[0] == (
        "alignment-wide resolved base composition omits at least one nucleotide, so TN93 assumptions break"
    )
    assert [
        (row.left_identifier, row.right_identifier, row.warning_kind)
        for row in report.warning_rows
    ] == [
        ("A", "B", "undefined-corrected-distance"),
        ("A", "C", "undefined-corrected-distance"),
        ("A", "D", "undefined-corrected-distance"),
        ("B", "C", "undefined-corrected-distance"),
        ("B", "D", "undefined-corrected-distance"),
        ("C", "D", "undefined-corrected-distance"),
    ]


def test_build_distance_tree_blocks_before_tree_inference_on_impossible_distances() -> (
    None
):
    with pytest.raises(
        InvalidAlignmentError,
        match="blocked before tree inference",
    ) as error:
        build_distance_tree(
            fixture("example_alignment_distance.fasta"),
            method="neighbor-joining",
            model="k80",
        )

    assert "A/C (infinite-corrected-distance)" in error.value.message

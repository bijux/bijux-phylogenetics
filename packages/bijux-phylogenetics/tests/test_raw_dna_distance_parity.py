from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.distance import compute_pairwise_genetic_distance_matrix
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


def test_raw_distance_alias_matches_existing_p_distance_surface() -> None:
    raw_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="raw",
    )
    p_distance_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="p-distance",
    )

    assert raw_report.model == "p-distance"
    assert raw_report.identifiers == p_distance_report.identifiers
    assert raw_report.pairs == p_distance_report.pairs


def test_raw_distance_alias_supports_complete_deletion_on_gapped_alignment() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        model="raw",
        gap_handling="complete-deletion",
    )

    assert [
        (
            pair.left_identifier,
            pair.right_identifier,
            pair.distance,
            pair.comparable_sites,
        )
        for pair in report.pairs
    ] == [
        ("A", "A", 0.0, 4),
        ("A", "B", 0.0, 4),
        ("A", "C", 0.25, 4),
        ("A", "D", 0.0, 4),
        ("B", "B", 0.0, 4),
        ("B", "C", 0.25, 4),
        ("B", "D", 0.0, 4),
        ("C", "C", 0.0, 4),
        ("C", "D", 0.25, 4),
        ("D", "D", 0.0, 4),
    ]


def test_raw_distance_alias_tracks_ignored_ambiguity_sites_explicitly() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="raw",
        gap_handling="pairwise-deletion",
        ambiguity_policy="ignore",
    )

    pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )

    assert pair.distance == 0.0
    assert pair.comparable_sites == 4
    assert pair.ambiguity_sites == 1
    assert pair.skipped_sites == 2


def test_raw_distance_alias_rejects_unequal_length_alignment_input() -> None:
    with pytest.raises(InvalidAlignmentError):
        compute_pairwise_genetic_distance_matrix(
            fixture("example_alignment_invalid_lengths.fasta"),
            model="raw",
        )


def test_cli_alignment_distance_matrix_accepts_raw_model_alias(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance.tsv"

    exit_code = main(
        [
            "alignment",
            "distance-matrix",
            str(fixture("example_alignment_distance.fasta")),
            "--model",
            "raw",
            "--out",
            str(output_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "p-distance"
    assert output_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tdistance\tcomparable_sites",
        "A\tA\t0\t8",
        "A\tB\t0.125\t8",
        "A\tC\t0.5\t8",
    ]

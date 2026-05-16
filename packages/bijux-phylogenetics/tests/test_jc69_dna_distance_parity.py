from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.distance import compute_pairwise_genetic_distance_matrix

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


def test_jc69_alias_matches_existing_jukes_cantor_surface() -> None:
    jc69_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="jc69",
    )
    jukes_cantor_report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance.fasta"),
        model="jukes-cantor",
    )

    assert jc69_report.model == "jukes-cantor"
    assert jc69_report.identifiers == jukes_cantor_report.identifiers
    assert jc69_report.pairs == jukes_cantor_report.pairs


def test_jc69_alias_supports_complete_deletion_on_gapped_alignment() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_gaps.fasta"),
        model="jc69",
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
        ("A", "C", 0.304098831081123, 4),
        ("A", "D", 0.0, 4),
        ("B", "B", 0.0, 4),
        ("B", "C", 0.304098831081123, 4),
        ("B", "D", 0.0, 4),
        ("C", "C", 0.0, 4),
        ("C", "D", 0.304098831081123, 4),
        ("D", "D", 0.0, 4),
    ]


def test_jc69_alias_tracks_ignored_ambiguity_sites_explicitly() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_ambiguity.fasta"),
        model="jc69",
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


def test_jc69_alias_distinguishes_undefined_and_infinite_saturation() -> None:
    report = compute_pairwise_genetic_distance_matrix(
        fixture("example_alignment_distance_saturated.fasta"),
        model="jc69",
    )

    undefined_pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "A" and row.right_identifier == "B"
    )
    infinite_pair = next(
        row
        for row in report.pairs
        if row.left_identifier == "B" and row.right_identifier == "C"
    )

    assert undefined_pair.distance is None
    assert undefined_pair.saturated is True
    assert (
        undefined_pair.saturation_reason
        == "p-distance exceeds the Jukes-Cantor correction range, so the corrected distance is undefined"
    )
    assert infinite_pair.distance is None
    assert infinite_pair.saturated is True
    assert (
        infinite_pair.saturation_reason
        == "p-distance is at the Jukes-Cantor correction limit, so the corrected distance tends to infinity"
    )


def test_cli_alignment_distance_matrix_accepts_jc69_model_alias(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "distance.tsv"

    exit_code = main(
        [
            "alignment",
            "distance-matrix",
            str(fixture("example_alignment_distance.fasta")),
            "--model",
            "jc69",
            "--out",
            str(output_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)

    assert exit_code == 0
    assert payload["metrics"]["model"] == "jukes-cantor"
    assert output_path.read_text(encoding="utf-8").splitlines()[:4] == [
        "left_identifier\tright_identifier\tdistance\tcomparable_sites",
        "A\tA\t0\t8",
        "A\tB\t0.136741167595466\t8",
        "A\tC\t0.823959216501082\t8",
    ]

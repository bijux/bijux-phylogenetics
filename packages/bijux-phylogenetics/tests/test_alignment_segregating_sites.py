from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.io.fasta import (
    compute_alignment_segregating_site_report,
    write_alignment_segregating_site_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).resolve().parent / "fixtures" / "alignments" / name


def test_compute_alignment_segregating_site_report_matches_clean_fixture() -> None:
    report = compute_alignment_segregating_site_report(
        fixture("example_alignment_distance.fasta")
    )

    assert report.sequence_count == 4
    assert report.alignment_length == 8
    assert report.segregating_site_positions == [1, 2, 3, 4, 8]
    assert report.ambiguity_policy == (
        "ambiguity states segregate only when they are surely incompatible with another observed state"
    )
    assert report.gap_policy == (
        "internal gap characters can create segregating sites against known or incompatible ambiguous states"
    )
    assert report.missing_data_policy == (
        "explicit missing characters do not create segregating sites"
    )
    assert report.trailing_gap_policy == (
        "leading and trailing gap runs are normalized to N before ape-style segregating-site detection"
    )
    assert [row.position for row in report.rows] == [1, 2, 3, 4, 8]
    assert report.rows[0].effective_states == "A|A|T|T"
    assert report.warnings == []


def test_compute_alignment_segregating_site_report_handles_gap_and_ambiguity_policy() -> (
    None
):
    gap_report = compute_alignment_segregating_site_report(
        fixture("example_alignment_distance_gaps.fasta")
    )
    ambiguity_report = compute_alignment_segregating_site_report(
        fixture("example_alignment_ambiguity.fasta")
    )

    assert gap_report.segregating_site_positions == [2, 5, 6]
    assert [row.effective_states for row in gap_report.rows] == [
        "C|C|T|-",
        "A|T|C|N",
        "A|A|G|A",
    ]
    assert ambiguity_report.segregating_site_positions == [5]
    assert ambiguity_report.rows[0].original_states == "N|R|-"
    assert ambiguity_report.rows[0].effective_states == "N|R|-"


def test_compute_alignment_segregating_site_report_handles_invariant_and_missing_cases() -> (
    None
):
    invariant_report = compute_alignment_segregating_site_report(
        fixture("example_alignment_invariant.fasta")
    )
    missing_report = compute_alignment_segregating_site_report(
        fixture("example_alignment_missingness.fasta")
    )
    all_gap_report = compute_alignment_segregating_site_report(
        fixture("example_alignment_all_gap_missing.fasta")
    )

    assert invariant_report.segregating_site_positions == []
    assert missing_report.segregating_site_positions == []
    assert all_gap_report.segregating_site_positions == []
    assert all_gap_report.warnings == [
        "alignment contains no canonical A/C/G/T residues, so ape-style segregating-site detection can only reflect ambiguity, gap, and missing states"
    ]


def test_compute_alignment_segregating_site_report_matches_one_variable_fixture() -> (
    None
):
    report = compute_alignment_segregating_site_report(
        fixture("example_alignment_one_variable_site.fasta")
    )

    assert report.segregating_site_positions == [4]
    assert report.rows[0].position == 4
    assert report.rows[0].effective_states == "A|A|T"


def test_write_alignment_segregating_site_table_writes_review_rows(
    tmp_path: Path,
) -> None:
    report = compute_alignment_segregating_site_report(
        fixture("example_alignment_distance_gaps.fasta")
    )
    output_path = tmp_path / "segregating-sites.tsv"

    write_alignment_segregating_site_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert (
        lines[0]
        == "position\toriginal_states\teffective_states\tknown_state_count\tambiguity_state_count\tgap_count\tmissing_count"
    )
    assert "2\tC|C|T|-\tC|C|T|-\t3\t0\t1\t0" in lines


def test_cli_alignment_segregating_sites_writes_site_table_and_json(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "segregating-sites.tsv"

    exit_code = main(
        [
            "alignment",
            "segregating-sites",
            str(fixture("example_alignment_distance_gaps.fasta")),
            "--site-table-out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["segregating_site_count"] == 3
    assert payload["data"]["segregating_site_positions"] == [2, 5, 6]
    assert payload["data"]["rows"][0]["position"] == 2
    assert output_path.exists()


def test_cli_alignment_segregating_sites_warns_on_all_gap_alignment(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "segregating-sites.tsv"

    exit_code = main(
        [
            "alignment",
            "segregating-sites",
            str(fixture("example_alignment_all_gap_missing.fasta")),
            "--site-table-out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["segregating_site_count"] == 0
    assert payload["warnings"] == [
        "alignment contains no canonical A/C/G/T residues, so ape-style segregating-site detection can only reflect ambiguity, gap, and missing states"
    ]
    assert output_path.exists()

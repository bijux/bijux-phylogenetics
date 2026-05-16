from __future__ import annotations

import json
from pathlib import Path

import pytest

from bijux_phylogenetics.command_line import main
from bijux_phylogenetics.io.fasta import (
    compute_alignment_base_frequency_report,
    write_alignment_base_frequency_table,
)


def fixture(name: str) -> Path:
    return Path(__file__).resolve().parent / "fixtures" / "alignments" / name


def _frequency_by_state(report) -> dict[str, float]:
    return {row.state: row.frequency for row in report.alignment_rows}


def test_compute_alignment_base_frequency_report_matches_lowercase_fixture() -> None:
    report = compute_alignment_base_frequency_report(
        fixture("example_alignment_lowercase.fasta")
    )

    assert report.sequence_count == 3
    assert report.alignment_length == 6
    assert report.ambiguity_policy == "count ambiguity codes as literal states"
    assert report.gap_policy == "count gap characters as literal states"
    assert (
        report.missing_data_policy
        == "count explicit missing characters as literal states"
    )
    assert _frequency_by_state(report) == {
        "a": pytest.approx(6 / 18),
        "c": pytest.approx(2 / 18),
        "g": pytest.approx(3 / 18),
        "t": pytest.approx(5 / 18),
        "r": 0.0,
        "m": 0.0,
        "w": 0.0,
        "s": 0.0,
        "k": 0.0,
        "y": 0.0,
        "v": 0.0,
        "h": 0.0,
        "d": 0.0,
        "b": 0.0,
        "n": pytest.approx(2 / 18),
        "-": 0.0,
        "?": 0.0,
    }
    lower_c_rows = {
        row.state: row.frequency
        for row in report.per_sequence_rows
        if row.identifier == "lower_c"
    }
    assert lower_c_rows["a"] == pytest.approx(1 / 6)
    assert lower_c_rows["t"] == pytest.approx(2 / 6)
    assert lower_c_rows["g"] == pytest.approx(1 / 6)
    assert lower_c_rows["n"] == pytest.approx(2 / 6)
    assert report.composition_outliers == []
    assert report.warnings == []


def test_compute_alignment_base_frequency_report_retains_ambiguity_gap_and_missing_states() -> (
    None
):
    report = compute_alignment_base_frequency_report(
        fixture("example_alignment_ambiguity.fasta")
    )

    assert _frequency_by_state(report)["r"] == pytest.approx(1 / 18)
    assert _frequency_by_state(report)["n"] == pytest.approx(1 / 18)
    assert _frequency_by_state(report)["-"] == pytest.approx(1 / 18)
    assert _frequency_by_state(report)["?"] == pytest.approx(3 / 18)


def test_compute_alignment_base_frequency_report_warns_on_all_gap_missing_alignment() -> (
    None
):
    report = compute_alignment_base_frequency_report(
        fixture("example_alignment_all_gap_missing.fasta")
    )

    assert _frequency_by_state(report)["-"] == pytest.approx(8 / 18)
    assert _frequency_by_state(report)["?"] == pytest.approx(10 / 18)
    assert all(
        _frequency_by_state(report)[state] == 0.0 for state in ["a", "c", "g", "t", "n"]
    )
    assert report.composition_outliers == []
    assert report.warnings == [
        "alignment contains no canonical A/C/G/T residues, so ape-style base frequencies reflect only ambiguity, gap, and missing states"
    ]


def test_write_alignment_base_frequency_table_writes_combined_scopes(
    tmp_path: Path,
) -> None:
    report = compute_alignment_base_frequency_report(
        fixture("example_alignment_lowercase.fasta")
    )
    output_path = tmp_path / "base-frequency.tsv"

    write_alignment_base_frequency_table(output_path, report)

    lines = output_path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "scope\tidentifier\tstate\tcount\tfrequency"
    assert "alignment\t\ta\t6\t0.333333333333333" in lines
    assert "sequence\tlower_c\tn\t2\t0.333333333333333" in lines


def test_cli_alignment_composition_writes_base_frequency_table_and_json(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "base-frequency.tsv"

    exit_code = main(
        [
            "alignment",
            "composition",
            str(fixture("example_alignment_lowercase.fasta")),
            "--base-frequency-out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["base_frequency_state_count"] == 17
    assert payload["metrics"]["composition_outlier_count"] == 0
    assert payload["data"]["alignment_state_frequencies"][0] == {
        "scope": "alignment",
        "identifier": None,
        "state": "a",
        "count": 6,
        "frequency": 0.333333333333333,
    }
    assert output_path.exists()


def test_cli_alignment_composition_flags_outliers_and_all_gap_warning(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "base-frequency.tsv"

    exit_code = main(
        [
            "alignment",
            "composition",
            str(fixture("example_alignment_all_gap_missing.fasta")),
            "--base-frequency-out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["warnings"] == [
        "alignment contains no canonical A/C/G/T residues, so ape-style base frequencies reflect only ambiguity, gap, and missing states"
    ]
    assert payload["metrics"]["composition_outlier_count"] == 0
    assert output_path.exists()

    exit_code = main(
        [
            "alignment",
            "composition",
            str(fixture("example_alignment_gc_outlier.fasta")),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["composition_outlier_count"] == 1
    assert payload["data"]["composition_outliers"][0]["identifier"] == "C"

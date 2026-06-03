from __future__ import annotations

import csv
import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

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


def test_cli_report_supplementary_diversification_table_writes_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-diversification.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-diversification-table",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions.tsv")),
            "--clade-model",
            "birth-death",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 3
    assert payload["metrics"]["clade_model"] == "birth-death"
    assert payload["metrics"]["high_clade_count"] == 1
    assert payload["metrics"]["low_clade_count"] == 1
    assert payload["metrics"]["sampling_metadata_complete"] is True
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert rows[0]["sampling_fraction"] == "0.75"


def test_cli_report_supplementary_diversification_table_surfaces_sampling_warnings(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-diversification-incomplete.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-diversification-table",
            "--tree",
            str(fixture("example_tree.nwk")),
            "--metadata",
            str(fixture("example_sampling_fractions_incomplete.tsv")),
            "--clade-model",
            "yule",
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 3
    assert payload["metrics"]["clade_model"] == "yule"
    assert payload["metrics"]["sampling_metadata_complete"] is False
    assert payload["metrics"]["warning_count"] >= 2
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3
    assert "B:missing-sampling-fraction:<missing>" in rows[0]["sampling_invalid_rows"]

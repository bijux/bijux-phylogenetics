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


def test_cli_report_supplementary_alignment_table_writes_table(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-alignment.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-alignment-table",
            "--alignment",
            str(fixture("example_taxon_workflow_alignment.fasta")),
            "--filtered-alignment",
            str(fixture("example_taxon_workflow_filtered_alignment.fasta")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert output_path.exists()
    assert payload["metrics"]["row_count"] == 3
    assert payload["metrics"]["retained_sequence_count"] == 2
    assert payload["metrics"]["removed_sequence_count"] == 1
    assert payload["metrics"]["filtered_only_sequence_count"] == 0
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["sequence_id"] == "A"
    assert "original_missing_fraction" in rows[0]

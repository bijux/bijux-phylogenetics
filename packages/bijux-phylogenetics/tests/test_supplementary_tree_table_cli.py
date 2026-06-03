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


def test_cli_report_supplementary_tree_table_writes_table(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-tree.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-tree-table",
            "--tree",
            str(fixture("example_tree_support_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert output_path.exists()
    assert payload["metrics"]["row_count"] == 1
    assert payload["metrics"]["tip_count"] == 4
    assert payload["metrics"]["supported_branch_count"] == 3
    assert payload["metrics"]["polytomy_count"] == 0
    assert payload["metrics"]["warning_count"] == 0
    assert payload["metrics"]["ultrametric"] is True
    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows[0]["tree_source"].endswith("example_tree_support_left.nwk")

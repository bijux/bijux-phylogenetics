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


def test_cli_report_supplementary_clade_support_table_writes_frequency_metrics(
    tmp_path: Path, capsys
) -> None:
    output_path = tmp_path / "supplementary-clade-support.tsv"

    exit_code = main(
        [
            "report",
            "supplementary-clade-support-table",
            "--tree",
            str(fixture("example_tree_support_left.nwk")),
            "--comparison-tree-set",
            str(fixture("example_tree_set_left.nwk")),
            "--out",
            str(output_path),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["row_count"] == 3
    assert payload["metrics"]["supported_clade_count"] == 3
    assert payload["metrics"]["frequency_scored_clade_count"] == 3
    assert payload["metrics"]["frequency_partial_support_count"] == 2
    assert payload["metrics"]["frequency_absent_clade_count"] == 0
    assert payload["metrics"]["frequency_unscored_clade_count"] == 0

    with output_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert any(row["descendant_taxa"] == "A|B" for row in rows)

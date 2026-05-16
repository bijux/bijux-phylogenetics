from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_cli_compare_support_reports_support_aware_conflicts_and_table_output(
    capsys, tmp_path: Path
) -> None:
    table_path = tmp_path / "support.tsv"

    exit_code = main(
        [
            "compare",
            "support",
            str(fixture("example_tree_support_conflict_left.nwk")),
            str(fixture("example_tree_support_conflict_right.nwk")),
            "--out",
            str(table_path),
            "--json",
        ]
    )
    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["shared_clades"] == 1
    assert payload["metrics"]["high_support_conflicts"] == 1
    assert payload["metrics"]["low_support_disagreements"] == 1
    assert payload["metrics"]["moderate_support_disagreements"] == 0
    assert payload["outputs"] == [str(table_path)]
    assert [
        row["conflict_classification"] for row in payload["data"]["conflicting_clades"]
    ] == ["low_support_disagreement", "high_support_conflict"]
    assert table_path.read_text(encoding="utf-8").startswith(
        "split_id\trow_kind\tcomparison_status\tleft_present\tright_present\t"
    )

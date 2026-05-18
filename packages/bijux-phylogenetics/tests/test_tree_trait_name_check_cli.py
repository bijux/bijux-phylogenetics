from __future__ import annotations

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


def test_cli_traits_name_check_json_output_and_table(tmp_path: Path, capsys) -> None:
    output_path = tmp_path / "trait-name-mismatches.tsv"
    exit_code = main(
        [
            "traits",
            "name-check",
            str(fixture("example_tree.nwk")),
            str(fixture("example_traits.tsv")),
            "--out",
            str(output_path),
            "--json",
        ]
    )

    captured = capsys.readouterr()
    payload = json.loads(captured.out)
    assert exit_code == 0
    assert payload["metrics"]["tree_not_data_count"] == 1
    assert payload["metrics"]["data_not_tree_count"] == 1
    assert payload["metrics"]["compatible"] is False
    assert payload["data"]["tree_not_data"] == ["D"]
    assert payload["data"]["data_not_tree"] == ["E"]
    assert output_path.exists()

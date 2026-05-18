from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main


def fixture(name: str) -> Path:
    return Path(__file__).parent / "fixtures" / "trees" / name


def test_speciational_cli_reports_sigma_squared(tmp_path: Path, capsys) -> None:
    output = tmp_path / "speciational.tsv"
    exit_code = main(
        [
            "simulate",
            "traits-speciational",
            str(fixture("example_tree_internal_long_branch.nwk")),
            "--root-state",
            "1.5",
            "--sigma-squared",
            "0.25",
            "--seed",
            "11",
            "--out",
            str(output),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["sigma_squared"] == 0.25
    assert payload["data"]["model"] == "speciational"
    assert output.read_text(encoding="utf-8").splitlines()[0] == "taxon\tvalue"

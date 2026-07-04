from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_consistency_index_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "consistency-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "consistency-index",
            str(fixture("fitch_tree.nwk")),
            str(fixture("consistency_index_matrix.tsv")),
            "--method",
            "fitch",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "parsimony-consistency-index"
    assert payload["metrics"]["method"] == "fitch"
    assert payload["metrics"]["minimum_possible_steps_total"] == 3.0
    assert payload["metrics"]["observed_steps_total"] == 4.0
    assert payload["metrics"]["consistency_index"] == 0.75
    assert (out_dir / "character_indices.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_consistency_index_cli_reports_constant_only_policy(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "consistency-index",
            str(fixture("fitch_tree.nwk")),
            str(fixture("consistency_index_constant_matrix.tsv")),
            "--method",
            "fitch",
            "--out-dir",
            str(tmp_path / "consistency-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["consistency_index"] is None
    assert payload["metrics"]["undefined_reason"] == "no_variable_characters"

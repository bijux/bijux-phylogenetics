from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_rescaled_consistency_index_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "rc-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "rescaled-consistency-index",
            str(fixture("fitch_tree.nwk")),
            str(fixture("rescaled_consistency_index_matrix.tsv")),
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
    assert payload["metrics"]["algorithm"] == "parsimony-rescaled-consistency-index"
    assert payload["metrics"]["method"] == "fitch"
    assert payload["metrics"]["ci"] == 0.75
    assert payload["metrics"]["ri"] == 0.5
    assert payload["metrics"]["rc"] == 0.375
    assert (out_dir / "character_indices.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_rescaled_consistency_index_cli_reports_undefined_policy(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "rescaled-consistency-index",
            str(fixture("fitch_tree.nwk")),
            str(fixture("rescaled_consistency_index_constant_matrix.tsv")),
            "--method",
            "fitch",
            "--out-dir",
            str(tmp_path / "rc-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["rc"] is None
    assert (
        payload["metrics"]["undefined_reason"]
        == "no_variable_characters|no_defined_retention_characters"
    )

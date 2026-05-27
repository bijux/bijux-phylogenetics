from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_acctran_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "acctran-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "acctran",
            str(fixture("acctran_tree_5_taxa.nwk")),
            str(fixture("acctran_ambiguous_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "acctran"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["character_count"] == 1
    assert payload["metrics"]["total_steps"] == 2
    assert (out_dir / "steps.tsv").is_file()
    assert (out_dir / "resolved_states.tsv").is_file()
    assert (out_dir / "branch_changes.tsv").is_file()
    assert (out_dir / "run.json").is_file()

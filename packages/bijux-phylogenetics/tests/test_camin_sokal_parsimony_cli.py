from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_camin_sokal_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "camin-sokal-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "camin-sokal",
            str(fixture("camin_sokal_tree_5_taxa.nwk")),
            str(fixture("camin_sokal_binary_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "camin-sokal"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["character_count"] == 3
    assert payload["metrics"]["root_state"] == "0"
    assert payload["metrics"]["total_gains"] == 3
    assert (out_dir / "steps.tsv").is_file()
    assert (out_dir / "branch_changes.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_camin_sokal_cli_reports_multistate_binarization_errors(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "camin-sokal",
            str(fixture("camin_sokal_tree_5_taxa.nwk")),
            str(fixture("camin_sokal_multistate_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "camin-sokal-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_matrix_multistate_not_binarized"
    assert (
        payload["errors"][0]["details"]["character_id"] == "char01_needs_binarization"
    )

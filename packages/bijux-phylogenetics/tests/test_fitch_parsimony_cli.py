from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_fitch_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "fitch-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "fitch",
            str(fixture("fitch_tree.nwk")),
            str(fixture("fitch_binary_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "unordered-fitch"
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["character_count"] == 2
    assert payload["metrics"]["total_steps"] == 2
    assert (out_dir / "steps.tsv").is_file()
    assert (out_dir / "node_state_sets.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_fitch_cli_reports_structured_missing_taxa_errors(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "fitch",
            str(fixture("fitch_tree.nwk")),
            str(fixture("fitch_missing_taxon_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "fitch-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_matrix_missing_taxa"
    assert payload["errors"][0]["details"]["missing_taxa"] == ["D"]


def test_phylo_parsimony_fitch_cli_honors_explicit_taxon_column(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "fitch-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "fitch",
            str(fixture("fitch_tree.nwk")),
            str(fixture("fitch_binary_species_matrix.tsv")),
            "--taxon-column",
            "species",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["character_count"] == 1
    assert payload["metrics"]["total_steps"] == 1

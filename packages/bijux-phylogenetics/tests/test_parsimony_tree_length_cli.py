from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_tree_length_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "tree-length-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "tree-length",
            str(fixture("fitch_tree.nwk")),
            str(fixture("fitch_binary_matrix.tsv")),
            "--method",
            "fitch",
            "--character-weights",
            str(fixture("fitch_character_weights.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "parsimony-tree-length"
    assert payload["metrics"]["method"] == "fitch"
    assert payload["metrics"]["raw_total_score"] == 2.0
    assert payload["metrics"]["total_score"] == 3.5
    assert (out_dir / "character_scores.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_tree_length_cli_reports_invalid_weight_errors(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "tree-length",
            str(fixture("fitch_tree.nwk")),
            str(fixture("fitch_binary_matrix.tsv")),
            "--method",
            "fitch",
            "--character-weights",
            str(fixture("parsimony_invalid_character_weight.tsv")),
            "--out-dir",
            str(tmp_path / "tree-length-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_character_weight_invalid_value"


def test_phylo_parsimony_tree_length_cli_requires_sankoff_cost_matrix(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "tree-length",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            "--method",
            "sankoff",
            "--out-dir",
            str(tmp_path / "tree-length-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_tree_length_cost_matrix_required"


def test_phylo_parsimony_tree_length_cli_rejects_asymmetric_sankoff_cost_matrix(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "tree-length",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            "--method",
            "sankoff",
            "--cost-matrix",
            str(fixture("sankoff_asymmetric_cost_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "tree-length-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_cost_matrix_asymmetric"


def test_phylo_parsimony_tree_length_cli_allows_asymmetric_sankoff_cost_matrix(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "tree-length",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            "--method",
            "sankoff",
            "--cost-matrix",
            str(fixture("sankoff_asymmetric_cost_matrix.tsv")),
            "--allow-asymmetric-costs",
            "--out-dir",
            str(tmp_path / "tree-length-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["method"] == "sankoff"

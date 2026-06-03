from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_sankoff_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "sankoff-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "sankoff",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            str(fixture("sankoff_cost_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "sankoff"
    assert payload["metrics"]["taxon_count"] == 5
    assert payload["metrics"]["character_count"] == 1
    assert payload["metrics"]["total_cost"] == 3.0
    assert payload["metrics"]["validation_warning_count"] == 0
    assert (out_dir / "steps.tsv").is_file()
    assert (out_dir / "node_costs.tsv").is_file()
    assert (out_dir / "selected_states.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_sankoff_cli_reports_negative_cost_matrix_errors(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "sankoff",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            str(fixture("sankoff_negative_cost_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "sankoff-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_cost_matrix_negative_cost"
    assert payload["errors"][0]["details"]["row_state"] == "red"
    assert payload["errors"][0]["details"]["column_state"] == "blue"


def test_phylo_parsimony_sankoff_cli_surfaces_validation_warnings(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "sankoff",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            str(fixture("sankoff_unused_state_cost_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "sankoff-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["validation_warning_count"] == 1
    assert payload["warnings"] == [
        "sankoff cost matrix includes states that are not observed in the current character matrix"
    ]


def test_phylo_parsimony_sankoff_cli_allows_asymmetry_explicitly(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "sankoff",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            str(fixture("sankoff_asymmetric_cost_matrix.tsv")),
            "--allow-asymmetric-costs",
            "--out-dir",
            str(tmp_path / "sankoff-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"


def test_phylo_parsimony_sankoff_cli_rejects_asymmetry_by_default(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "sankoff",
            str(fixture("sankoff_tree_5_taxa.nwk")),
            str(fixture("sankoff_character_matrix.tsv")),
            str(fixture("sankoff_asymmetric_cost_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "sankoff-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_cost_matrix_asymmetric"

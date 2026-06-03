from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_wagner_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "wagner-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "wagner",
            str(fixture("fitch_tree.nwk")),
            str(fixture("wagner_ordinal_matrix.tsv")),
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["algorithm"] == "ordered-wagner"
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["character_count"] == 2
    assert payload["metrics"]["total_cost"] == 6
    assert (out_dir / "steps.tsv").is_file()
    assert (out_dir / "node_costs.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_wagner_cli_reports_missing_state_order_errors(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "wagner",
            str(fixture("fitch_tree.nwk")),
            str(fixture("wagner_named_state_matrix.tsv")),
            "--out-dir",
            str(tmp_path / "wagner-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert payload["errors"][0]["code"] == "parsimony_state_order_required"
    assert payload["errors"][0]["details"]["character_id"] == "char01_gradient"


def test_phylo_parsimony_wagner_cli_accepts_explicit_state_order(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "wagner",
            str(fixture("fitch_tree.nwk")),
            str(fixture("wagner_named_state_matrix.tsv")),
            "--state-order",
            "low,medium,high,very_high",
            "--out-dir",
            str(tmp_path / "wagner-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["metrics"]["total_cost"] == 6

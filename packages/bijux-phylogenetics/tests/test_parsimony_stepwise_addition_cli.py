from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_stepwise_addition_cli_writes_governed_artifacts(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "stepwise-addition-cli"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "stepwise-addition",
            str(fixture("nni_search_matrix.tsv")),
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
    assert payload["metrics"]["algorithm"] == "greedy-stepwise-addition-tree"
    assert payload["metrics"]["method"] == "fitch"
    assert payload["metrics"]["character_count"] == 2
    assert payload["metrics"]["taxon_count"] == 4
    assert payload["metrics"]["final_score"] == 2.0
    assert payload["metrics"]["insertion_step_count"] == 2
    assert (out_dir / "tree.nwk").is_file()
    assert (out_dir / "trace.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_stepwise_addition_cli_reports_taxon_order_mismatch(
    tmp_path: Path,
    capsys,
) -> None:
    exit_code = main(
        [
            "phylo",
            "parsimony",
            "stepwise-addition",
            str(fixture("nni_search_matrix.tsv")),
            "--method",
            "fitch",
            "--insertion-order",
            "A,B,C,X",
            "--out-dir",
            str(tmp_path / "stepwise-addition-cli"),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 2
    assert payload["status"] == "error"
    assert (
        payload["errors"][0]["code"]
        == "parsimony_stepwise_addition_insertion_order_taxa_mismatch"
    )

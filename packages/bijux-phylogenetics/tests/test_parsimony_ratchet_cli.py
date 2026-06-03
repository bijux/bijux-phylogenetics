from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_ratchet_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "ratchet"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "ratchet",
            str(fixture("ratchet_search_start_tree_5_taxa.nwk")),
            str(fixture("ratchet_search_matrix.tsv")),
            "--method",
            "fitch",
            "--cycle-count",
            "3",
            "--seed",
            "1",
            "--perturbed-character-count",
            "1",
            "--perturbation-factor",
            "2.0",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"] == {
        "algorithm": "parsimony-ratchet",
        "method": "fitch",
        "taxon_count": 5,
        "character_count": 3,
        "cycle_count": 3,
        "random_seed": 1,
        "perturbed_character_count": 1,
        "perturbation_factor": 2.0,
        "start_score": 5.0,
        "final_score": 4.0,
        "best_score": 4.0,
        "best_tree_history_count": 2,
    }
    assert (out_dir / "start_tree.nwk").is_file()
    assert (out_dir / "final_tree.nwk").is_file()
    assert (out_dir / "best_tree.nwk").is_file()
    assert (out_dir / "cycle_history.tsv").is_file()
    assert (out_dir / "best_tree_history.tsv").is_file()
    assert (out_dir / "run.json").is_file()

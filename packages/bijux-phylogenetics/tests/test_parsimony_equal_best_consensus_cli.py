from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_equal_best_consensus_cli_writes_governed_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "equal-best-consensus"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "equal-best-consensus",
            str(fixture("bootstrap_matrix.tsv")),
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
    assert payload["metrics"] == {
        "algorithm": "parsimony-equal-best-consensus",
        "method": "fitch",
        "taxon_count": 4,
        "character_count": 4,
        "candidate_tree_count": 15,
        "best_score": 5.0,
        "equal_best_tree_count": 5,
        "retained_equal_best_tree_count": 5,
        "retained_all_equal_best_trees": True,
    }
    assert (out_dir / "equal_best_trees.nwk").is_file()
    assert (out_dir / "equal_best_scores.tsv").is_file()
    assert (out_dir / "strict_consensus_tree.nwk").is_file()
    assert (out_dir / "majority_consensus_tree.nwk").is_file()
    assert (out_dir / "clade_frequencies.tsv").is_file()
    assert (out_dir / "run.json").is_file()


def test_phylo_parsimony_equal_best_consensus_cli_omits_partial_consensus_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    out_dir = tmp_path / "equal-best-consensus-truncated"

    exit_code = main(
        [
            "phylo",
            "parsimony",
            "equal-best-consensus",
            str(fixture("bootstrap_matrix.tsv")),
            "--method",
            "fitch",
            "--max-retained-equal-best-trees",
            "3",
            "--out-dir",
            str(out_dir),
            "--json",
        ]
    )
    payload = json.loads(capsys.readouterr().out)

    assert exit_code == 0
    assert payload["status"] == "ok"
    assert payload["metrics"]["equal_best_tree_count"] == 5
    assert payload["metrics"]["retained_equal_best_tree_count"] == 3
    assert payload["metrics"]["retained_all_equal_best_trees"] is False
    assert (out_dir / "equal_best_trees.nwk").is_file()
    assert (out_dir / "equal_best_scores.tsv").is_file()
    assert not (out_dir / "strict_consensus_tree.nwk").exists()
    assert not (out_dir / "majority_consensus_tree.nwk").exists()
    assert not (out_dir / "clade_frequencies.tsv").exists()
    assert (out_dir / "run.json").is_file()

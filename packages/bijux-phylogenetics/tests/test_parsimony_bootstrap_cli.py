from __future__ import annotations

import json
from pathlib import Path

from bijux_phylogenetics.command_line import main

FIXTURES = Path(__file__).parent / "fixtures" / "parsimony"


def fixture(name: str) -> Path:
    return FIXTURES / name


def test_phylo_parsimony_bootstrap_cli_writes_governed_seeded_outputs(
    tmp_path: Path,
    capsys,
) -> None:
    left_out_dir = tmp_path / "bootstrap-left"
    right_out_dir = tmp_path / "bootstrap-right"
    command = [
        "phylo",
        "parsimony",
        "bootstrap",
        str(fixture("bootstrap_matrix.tsv")),
        "--method",
        "fitch",
        "--replicate-count",
        "20",
        "--seed",
        "1",
        "--out-dir",
    ]

    left_exit_code = main([*command, str(left_out_dir), "--json"])
    left_payload = json.loads(capsys.readouterr().out)
    right_exit_code = main([*command, str(right_out_dir), "--json"])
    right_payload = json.loads(capsys.readouterr().out)

    assert left_exit_code == 0
    assert right_exit_code == 0
    assert left_payload["status"] == "ok"
    assert right_payload["status"] == "ok"
    assert left_payload["metrics"] == {
        "algorithm": "parsimony-bootstrap",
        "method": "fitch",
        "taxon_count": 4,
        "character_count": 4,
        "replicate_count": 20,
        "candidate_tree_count": 15,
        "reference_score": 5.0,
        "support_row_count": 2,
    }
    assert right_payload["metrics"] == left_payload["metrics"]
    assert (left_out_dir / "reference_tree.nwk").is_file()
    assert (left_out_dir / "replicate_trees.nwk").is_file()
    assert (left_out_dir / "replicate_scores.tsv").is_file()
    assert (left_out_dir / "replicate_draws.tsv").is_file()
    assert (left_out_dir / "clade_support.tsv").is_file()
    assert (left_out_dir / "consensus_tree.nwk").is_file()
    assert (left_out_dir / "clade_frequencies.tsv").is_file()
    assert (left_out_dir / "run.json").is_file()
    assert (left_out_dir / "replicate_trees.nwk").read_text(encoding="utf-8") == (
        right_out_dir / "replicate_trees.nwk"
    ).read_text(encoding="utf-8")
    assert (left_out_dir / "clade_support.tsv").read_text(encoding="utf-8") == (
        right_out_dir / "clade_support.tsv"
    ).read_text(encoding="utf-8")
